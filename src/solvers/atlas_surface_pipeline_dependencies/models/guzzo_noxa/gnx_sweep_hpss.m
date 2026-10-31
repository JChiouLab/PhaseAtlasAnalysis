function S = gnx_sweep_hpss(P, ell, kind, opts)
%GNX_SWEEP_HPSS  Half-peak search on a segment [0, ell] (packet §4h).
%
%   S = GNX_SWEEP_HPSS(P, ell, kind)
%   S = GNX_SWEEP_HPSS(P, ell, kind, opts)
%   P: gnx_model('parameters'); ell: segment length (Lh = L/2 for HPSS,
%   L/4 for MPSS quarters); kind: 'hpss' | 'quarter' (validity kind).
%   Segment masses n = (rho + dn/2)*ell (same densities on every segment, so
%   the grid cell is unchanged; run_stage3_noxa.m:218-220 with rho = 1).
%   Myxo's search (run_stage3_noxa.m) with our gates:
%     anchor bins centre-first 3 x 3 (:69-83), 15 constructed guesses
%     (3 amplitudes x 5 widths, :257-275, 292-302), anchor column then rows
%     with break on the first failure (:99-131), constructed-seed pass on
%     every unsolved cell (:141-152), up to 3 neighbour rescue passes
%     (:158-182). Warm starts evaluate the neighbour on the uniform
%     121-point mesh (try_from, :281-290). No anchor is an empty result.
%   A solve is accepted iff its gnx_validity record passes; a converged
%   peak-at-ell return gets one orientation re-solve from
%   gnx_guess('reverse', sol) (RC-7) and only a converged peak-at-0 result
%   is accepted.
%   Bounds (RC-8): at most opts.max_calls bvp4c calls per sweep (default
%   20000), then the sweep stops and keeps its accepted cells
%   (S.cap_reached); S.n_nmax counts calls whose mesh reached NMax;
%   opts.ckpt_path (optional): SOLS and records saved after every phase.
%   opts.solve_opts is passed to gnx_solve (tests inject bvp_fun there).
%   S.attempt_log (E2, phase B): one row per attempt (phase, cell, guess,
%   first solve, orientation re-solve, validity outcome), for the §5.3
%   attribution; checkpointed with SOLS. Recording only; no decision reads it.

if nargin < 4 || isempty(opts), opts = struct(); end
D = gnx_model('descriptor');
max_calls = getf(opts, 'max_calls', D.sweep.max_calls);
ckpt_path = getf(opts, 'ckpt_path', '');
solve_opts = getf(opts, 'solve_opts', struct());
verbose = getf(opts, 'verbose', true);

res = P.res;
dns = linspace(-P.dn_half_range, P.dn_half_range, res);
xg = linspace(0, ell, D.numerics.n_init);

S = struct('kind', kind, 'ell', ell, 'dn', dns, 'xg', xg, ...
    'ok_grid', false(res, res), 'SOLS', {cell(res, res)}, ...
    'validity', {cell(res, res)}, 'last_reason', strings(res, res), ...
    'n_attempts', zeros(res, res), 'n_calls', 0, 'n_nmax', 0, ...
    'cap_reached', false, 'reason_counts', struct(), 'anchor', [], ...
    'phase_log', {{}}, 'max_calls', max_calls, 'n_flip_resolves', 0, ...
    'cur_phase', 'anchor', 'cur_guess', '', 'attempt_log', []);
alog('reset');
t0 = tic;

% ---- anchor ----
q = round(res * [0.5 0.25 0.75]);
[bi, bj] = meshgrid(q, q);
bins = [bi(:), bj(:)];
anchor = [];
anchor_sol = [];
for b = 1:size(bins, 1)
    i = bins(b, 1); j = bins(b, 2);
    [sol, rec, S] = try_seeds(S, P, i, j, ell, xg, max_calls, solve_opts);
    if ~isempty(sol)
        anchor = [i, j]; anchor_sol = sol;
        S = keep(S, i, j, sol, rec);
        break
    end
    if S.cap_reached, break; end
end
march = ~isempty(anchor);
if march
    S.anchor = anchor;
else
    anchor = [1, 1];
end

% ---- anchor column ----
S.cur_phase = 'anchor_column';
j = anchor(2);
if march && ~S.cap_reached, dirs = [-1, 1]; else, dirs = []; end
for dir = dirs
    prev = anchor_sol;
    if dir < 0, rng_i = anchor(1)-1:-1:1; else, rng_i = anchor(1)+1:res; end
    for i = rng_i
        S.cur_guess = sprintf('w%d,%d', i - dir, j);
        [sol, rec, S] = try_from(S, P, i, j, ell, xg, prev, max_calls, solve_opts);
        if isempty(sol), break; end
        S = keep(S, i, j, sol, rec);
        prev = sol;
    end
end
S = phase_done(S, 'anchor_column', t0, ckpt_path, verbose);

% ---- rows ----
S.cur_phase = 'rows';
for i = 1:res
    if S.cap_reached, break; end
    if ~S.ok_grid(i, j), continue; end
    for dir = [-1, 1]
        prev = S.SOLS{i, j};
        if dir < 0, rng_j = j-1:-1:1; else, rng_j = j+1:res; end
        for jj = rng_j
            S.cur_guess = sprintf('w%d,%d', i, jj - dir);
            [sol, rec, S] = try_from(S, P, i, jj, ell, xg, prev, max_calls, solve_opts);
            if isempty(sol), break; end
            S = keep(S, i, jj, sol, rec);
            prev = sol;
        end
    end
end
S = phase_done(S, 'rows', t0, ckpt_path, verbose);

% ---- constructed-seed pass ----
S.cur_phase = 'constructed_seeds';
for i = 1:res
    for jj = 1:res
        if S.cap_reached, break; end
        if S.ok_grid(i, jj), continue; end
        [sol, rec, S] = try_seeds(S, P, i, jj, ell, xg, max_calls, solve_opts);
        if isempty(sol), continue; end
        S = keep(S, i, jj, sol, rec);
    end
end
S = phase_done(S, 'constructed_seeds', t0, ckpt_path, verbose);

% ---- rescue passes ----
for pass = 1:D.sweep.rescue_passes
    if S.cap_reached, break; end
    S.cur_phase = sprintf('rescue_%d', pass);
    n_new = 0;
    for i = 1:res
        for jj = 1:res
            if S.ok_grid(i, jj), continue; end
            nb = [i-1 jj; i+1 jj; i jj-1; i jj+1];
            nb = nb(all(nb >= 1, 2) & all(nb <= res, 2), :);
            for k = 1:size(nb, 1)
                if S.cap_reached, break; end
                if ~S.ok_grid(nb(k,1), nb(k,2)), continue; end
                prev = S.SOLS{nb(k,1), nb(k,2)};
                S.cur_guess = sprintf('w%d,%d', nb(k,1), nb(k,2));
                [sol, rec, S] = try_from(S, P, i, jj, ell, xg, prev, max_calls, solve_opts);
                if ~isempty(sol)
                    S = keep(S, i, jj, sol, rec);
                    n_new = n_new + 1; break
                end
            end
        end
    end
    S = phase_done(S, sprintf('rescue_%d', pass), t0, ckpt_path, verbose);
    if n_new == 0, break; end
end
S.wall_s = toc(t0);
S.attempt_log = alog('get');
end

% =========================================================================
function [n_R, n_B] = masses(P, S, i, j, ell)
n_R = (P.rho_R + S.dn(j) / 2) * ell;
n_B = (P.rho_B + S.dn(i) / 2) * ell;
end

function [sol, rec, S] = try_seeds(S, P, i, j, ell, xg, max_calls, solve_opts)
D = gnx_model('descriptor');
sol = []; rec = [];
[nR, nB] = masses(P, S, i, j, ell);
for amp_frac = D.sweep.amp_fracs
    for w = D.sweep.widths
        if S.n_calls >= max_calls, S.cap_reached = true; return; end
        S.cur_guess = sprintf('c%.1f/%.2f', amp_frac, w);
        guess = gnx_guess('constructed', xg, nR, nB, amp_frac, w, P);
        [s, r, S] = attempt(S, P, i, j, ell, guess, nR, nB, max_calls, solve_opts);
        if ~isempty(s), sol = s; rec = r; return; end
    end
end
end

function [sol, rec, S] = try_from(S, P, i, j, ell, xg, prev, max_calls, solve_opts)
sol = []; rec = [];
if S.n_calls >= max_calls, S.cap_reached = true; return; end
[nR, nB] = masses(P, S, i, j, ell);
y = deval(prev, xg);
guess = gnx_guess('profile', xg, y([1 3 5 7], :));
[sol, rec, S] = attempt(S, P, i, j, ell, guess, nR, nB, max_calls, solve_opts);
end

function [sol, rec, S] = attempt(S, P, i, j, ell, guess, nR, nB, max_calls, solve_opts)
% One constructed/warm guess: solve, orient (one re-solve), validate.
sol = []; rec = [];
ctx = struct('par', P, 'ell', ell, 'nR', nR, 'nB', nB);
[s, sr] = gnx_solve(guess, P, nR, nB, solve_opts);
S = count(S, sr);
S.n_attempts(i, j) = S.n_attempts(i, j) + 1;
row = row0();
row.i = i; row.j = j; row.phase = S.cur_phase; row.guess = S.cur_guess;
row.call = sr.call_index; row.conv = sr.converged; row.reason = sr.reason;
row.maxres = sr.maxres; row.n_x = sr.n_x; row.guard = sr.guard_hits; row.warn = sr.warn_id;
if ~sr.converged
    S = reason(S, i, j, sr.reason);
    row.outcome = sr.reason; alog('add', row);
    return;
end
if s.y(5, 1) < s.y(5, end)
    row.flip = true;
    if S.n_calls >= max_calls
        S.cap_reached = true; S = reason(S, i, j, 'cap');
        row.outcome = 'cap'; alog('add', row); return;
    end
    [s2, sr2] = gnx_solve(gnx_guess('reverse', s), P, nR, nB, solve_opts);
    S = count(S, sr2);
    S.n_flip_resolves = S.n_flip_resolves + 1;
    row.flip_conv = sr2.converged; row.flip_reason = sr2.reason; row.flip_maxres = sr2.maxres;
    if ~sr2.converged
        S = reason(S, i, j, ['orientation_resolve_' sr2.reason]);
        row.outcome = ['orientation_resolve_' sr2.reason]; alog('add', row);
        return;
    end
    s = s2; sr = sr2;
end
r = gnx_validity(S.kind, s, sr, ctx);
row.valid = r.valid;
if isfield(r.checks, 'positivity'), row.pos_min = r.checks.positivity.value; end
if isfield(r.checks, 'morphology')
    row.morph = char(r.checks.morphology.value);
    row.far_gap_rel = r.checks.morphology.far_trough_gap_rel;
end
if r.valid
    sol = s; rec = r;
    row.outcome = 'accepted';
else
    S = reason(S, i, j, first_failed(r));
    row.outcome = first_failed(r);
end
alog('add', row);
end

function row = row0()
row = struct('i', NaN, 'j', NaN, 'phase', '', 'guess', '', 'call', NaN, 'conv', false, ...
    'reason', '', 'maxres', NaN, 'n_x', NaN, 'guard', NaN, 'warn', '', 'flip', false, ...
    'flip_conv', false, 'flip_reason', '', 'flip_maxres', NaN, 'valid', false, ...
    'pos_min', NaN, 'morph', '', 'far_gap_rel', NaN, 'outcome', '');
end

function out = alog(action, row)
% Attempt log in a persistent buffer (no per-attempt copy of a growing array).
persistent B n
out = [];
switch action
    case 'reset'
        B = repmat(row0(), 1024, 1); n = 0;
    case 'add'
        if isempty(n), B = repmat(row0(), 1024, 1); n = 0; end
        n = n + 1;
        if n > numel(B), B(2 * numel(B)) = row0(); end
        B(n) = row;
    case 'get'
        if isempty(n), out = repmat(row0(), 0, 1); else, out = B(1:n); end
end
end

function S = count(S, sr)
S.n_calls = S.n_calls + 1;
if sr.nmax_reached, S.n_nmax = S.n_nmax + 1; end
end

function S = reason(S, i, j, why)
if isempty(why), why = 'unknown'; end
why = matlab.lang.makeValidName(why);
S.last_reason(i, j) = string(why);
if isfield(S.reason_counts, why)
    S.reason_counts.(why) = S.reason_counts.(why) + 1;
else
    S.reason_counts.(why) = 1;
end
end

function why = first_failed(r)
names = gnx_validity('required', r.kind);
why = 'unknown';
for k = 1:numel(names)
    if ~isfield(r.checks, names{k}) || ~r.checks.(names{k}).pass
        why = names{k};
        if strcmp(why, 'morphology') && isfield(r.checks, 'morphology')
            why = ['morphology_' char(r.checks.morphology.value)];
        end
        return;
    end
end
end

function S = keep(S, i, j, sol, rec)
S.ok_grid(i, j) = true;
S.SOLS{i, j} = sol;
S.validity{i, j} = rec;
S.last_reason(i, j) = "accepted";
end

function S = phase_done(S, name, t0, ckpt_path, verbose)
S.attempt_log = alog('get');
S.phase_log{end+1} = sprintf('%s: ok %d, calls %d, nmax %d, cap %d, %.1f s', ...
    name, nnz(S.ok_grid), S.n_calls, S.n_nmax, S.cap_reached, toc(t0));
if verbose, fprintf('[gnx_sweep_hpss %s ell=%g] %s\n', S.kind, S.ell, S.phase_log{end}); end
if ~isempty(ckpt_path)
    ok_grid = S.ok_grid; SOLS = S.SOLS; validity = S.validity;
    last_reason = S.last_reason; n_calls = S.n_calls; n_nmax = S.n_nmax;
    phase_log = S.phase_log; attempt_log = S.attempt_log; n_attempts = S.n_attempts;
    save(ckpt_path, 'ok_grid', 'SOLS', 'validity', 'last_reason', 'n_calls', ...
        'n_nmax', 'phase_log', 'attempt_log', 'n_attempts', '-v7.3');
end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
