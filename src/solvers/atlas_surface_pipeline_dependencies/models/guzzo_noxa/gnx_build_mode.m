function R = gnx_build_mode(P, n, out_dir, opts, varargin)
%GNX_BUILD_MODE  Higher-mode half-cell library, mode n: n copies of a reflected half-peak piece
%   (M. xanthus counterpart of build_mode_branch; the paper uses only no-peak, half-peak and
%   full-peak states, so modes n >= 3 are an extension).
%
%   R = GNX_BUILD_MODE(P, n, out_dir)
%   R = GNX_BUILD_MODE(P, n, out_dir, opts)
%       P = gnx_model('parameters', ...); n >= 2 integer (production n >= 3;
%       n = 2 is the full-peak state and is allowed for tests); out_dir = the
%       library folder, which must already hold branch_{HSS,HPSS,MPSS}_latest
%       .mat of the same P (the grid, Lh, res and model_version are copied
%       from branch_HSS_latest.mat and checked bit-equal to the HPSS and MPSS
%       siblings; never taken from the sweep).
%       1. The piece: gnx_sweep_hpss(P, P.L/(2n), 'quarter', so), unchanged,
%          so = opts.sweep with so.ckpt_path = <out_dir>/mode<n>_work/
%          checkpoint_piece_<run_id>.mat. Piece masses (rho + dn/2) * l are
%          1/n of the half masses at the same densities (see gnx_sweep_hpss).
%       2. Every found piece is stitched T-start (a trough at the pole) and
%          tested (gnx_build_mode('record', ...)); a cell is in the mask iff
%          every check passes; fail_reason = the first failing check or
%          'quarter_not_found'.
%       3. Writes branch_MODE<n>_<run_id>.mat and copies it to
%          branch_MODE<n>_latest.mat in out_dir. run_status 'empty' when the
%          mask is all false (the file is still written), else 'complete'.
%   opts: sweep (struct for gnx_sweep_hpss: max_calls, solve_opts (tests
%       inject bvp_fun there), verbose), run_id (default a timestamp),
%       n2_branch ('legacy' (default) | 'general'; which record runs at n = 2: 'legacy' is the
%       full-peak record of gnx_validity.m, 'general' the code for any n).
%   R: path, latest, n_valid, run_status, n_calls, ell, sweep_file.
%
%   rec = GNX_BUILD_MODE('record', n, qrec, solq, ctx)       general n
%   rec = GNX_BUILD_MODE('record_legacy', qrec, solq, ctx)   n = 2, the
%       mpss record of gnx_validity.m, copied unchanged (kind 'mpss').
%       ctx: par, ell_q, nR_half, nB_half (as gnx_build_mpss).
%       General n (kind 'mode'), checks in the n = 2 order:
%         quarter_valid  the piece's 'quarter' record passes
%         seam_slope     |Bm'(0)| of the piece (its peak end) <= mpss_slope
%         trough_slope   n >= 3 only: |Bm'(l)| (its trough end, an interior
%                        join for n >= 3) <= mpss_slope
%         symmetry       even n: max|ys - fliplr(ys)| of the T-start stitch
%                        (reflection about the half's centre, as at n = 2);
%                        odd n: max|ys_T - fliplr(ys_P)|, the T-start stitch
%                        against the reflected P-start stitch of the same
%                        piece (a T-start stitch of odd n is the reflection
%                        of the P-start one); both test the stitching
%                        arithmetic; <= mpss_sym
%         mass_doubling  dense mass of the stitch, the piece's Hermite
%                        interpolant at the T fold of n (n_dense - 1) + 1
%                        uniform points (n = 2: 2 n_dense - 1),
%                        relative to the half masses (rho + dn/2) Lh;
%                        <= mass_rel; records n_points
%         peak_fraction  the piece's argmax position (last maximal node)
%                        x*/l < 1 - 2 mpss_peak_lo. At n = 2 the
%                        stitch's first maximum is the mirrored copy's node
%                        at l - x*, so pf = (l - x*)/(2l) = 1/2 - x*/(2l) <= 1/2
%                        < mpss_peak_hi always, and pf > lo <=> x*/l < 1 - 2 lo
%         contrast       (max - min)/max of Bm on the stitch >= contrast
%       ('quarter_valid' and 'mass_doubling' keep the names used for n = 2: the
%       piece and the n-fold mass.)
%   c = GNX_BUILD_MODE('stitch', xq, yq, n, start)   c = {x, y}; 'T' | 'P' (identical
%       to sag_seed_shape 'stitch'); n = 2 'T' = the peak-to-peak stitch of the full-peak state
%   xm = GNX_BUILD_MODE('fold', x, l, start)   piece coordinate (identical to
%       sag_seed_shape 'fold')
%   names = GNX_BUILD_MODE('required', n)
%   ok = GNX_BUILD_MODE('verdict', rec)
%   path = GNX_BUILD_MODE('write', out_dir, n, P, data, meta)   the writer
%       (data: existence_mask, SOL, SOL_piece, validity, fail_reason,
%       acceptance, piece (ok_grid, last_reason, validity, phase_log,
%       cap_reached), n_calls, n_nmax, source_sweep_file, record_branch;
%       meta: run_id plus provenance)

if ischar(P)
    a = {}; if nargin >= 2, a{1} = n; end
    if nargin >= 3, a{2} = out_dir; end
    if nargin >= 4, a{3} = opts; end
    a = [a, varargin];
    switch P
        case 'record',        R = mode_record(a{:});
        case 'record_legacy', R = legacy_record(a{:});
        case 'stitch',        [R, x2] = stitch(a{:}); R = {R, x2};
        case 'fold',          R = fold(a{:});
        case 'required',      R = required_mode(a{:});
        case 'verdict',       R = verdict_any(a{:});
        case 'write',         R = write_mode(a{:});
        otherwise, error('GNX:MODE', 'Unknown action %s.', P);
    end
    return
end
if nargin < 4 || isempty(opts), opts = struct(); end
if ~(isscalar(n) && n == round(n) && n >= 2)
    error('GNX:MODE', 'n must be an integer >= 2.');
end
run_id = getf(opts, 'run_id', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
n2_branch = getf(opts, 'n2_branch', 'legacy');
so = getf(opts, 'sweep', struct());
sib = siblings(out_dir, P);
Lh = P.L / 2;
ell = P.L / (2 * n);
work = fullfile(out_dir, sprintf('mode%d_work', n));
if exist(work, 'dir') ~= 7, mkdir(work); end
so.ckpt_path = fullfile(work, sprintf('checkpoint_piece_%s.mat', run_id));
Q = gnx_sweep_hpss(P, ell, 'quarter', so);
if ~isequal(Q.dn, sib.dn_R) || ~isequal(Q.dn, sib.dn_B)
    error('GNX:MODE', 'The sweep''s axis differs from the sibling library grid.');
end
sweep_file = fullfile(work, sprintf('sweep_piece_%s.mat', run_id));
save(sweep_file, '-struct', 'Q', '-v7.3');
res = P.res;
data = struct();
data.existence_mask = false(res, res);
data.SOL = cell(res, res);
data.SOL_piece = cell(res, res);
data.validity = cell(res, res);
data.fail_reason = strings(res, res);
acc = struct('seam_slope', nan(res), 'trough_slope', nan(res), 'symmetry', nan(res), ...
    'mass_rel', nan(res, res, 2), 'peak_position', nan(res), 'contrast', nan(res), ...
    'n_mass_points', NaN);
Dd = gnx_model('descriptor');
acc.n_mass_points = n * (Dd.numerics.n_dense - 1) + 1;
use_legacy = n == 2 && strcmp(n2_branch, 'legacy');
for i = 1:res
    for j = 1:res
        if ~Q.ok_grid(i, j)
            data.fail_reason(i, j) = "quarter_not_found";
            continue;
        end
        ctx = struct('par', P, 'ell_q', ell, ...
            'nR_half', (P.rho_R + Q.dn(j) / 2) * Lh, ...
            'nB_half', (P.rho_B + Q.dn(i) / 2) * Lh);
        if use_legacy
            rec = legacy_record(Q.validity{i, j}, Q.SOLS{i, j}, ctx);
        else
            rec = mode_record(n, Q.validity{i, j}, Q.SOLS{i, j}, ctx);
        end
        data.validity{i, j} = rec;
        acc = note_acceptance(acc, rec, i, j);
        if rec.valid
            data.existence_mask(i, j) = true;
            data.SOL{i, j} = struct('x', rec.x, 'y', rec.y);
            data.SOL_piece{i, j} = Q.SOLS{i, j};
        else
            data.fail_reason(i, j) = string(first_failed(rec, n));
        end
    end
end
data.acceptance = acc;
data.piece = struct('ok_grid', Q.ok_grid, 'last_reason', Q.last_reason, ...
    'validity', {Q.validity}, 'phase_log', {Q.phase_log}, 'cap_reached', Q.cap_reached);
data.n_calls = Q.n_calls;
data.n_nmax = Q.n_nmax;
data.source_sweep_file = sweep_file;
if use_legacy, data.record_branch = 'legacy_n2'; else, data.record_branch = 'general'; end
meta = struct('run_id', run_id, 'matlab_version', version, 'opts_sweep', rmfield_safe(so, 'solve_opts'));
path_ = write_mode(out_dir, n, P, data, meta);
R = struct('path', path_, 'latest', fullfile(out_dir, sprintf('branch_MODE%d_latest.mat', n)), ...
    'n_valid', nnz(data.existence_mask), 'run_status', run_status_of(data.existence_mask), ...
    'n_calls', Q.n_calls, 'ell', ell, 'sweep_file', sweep_file);
fprintf('[gnx_build_mode n=%d ell=%.17g] %d/%d cells (%d pieces), %d calls, %s\n', n, ell, ...
    R.n_valid, res^2, nnz(Q.ok_grid), Q.n_calls, R.run_status);
end

% =========================================================================
function a = rmfield_safe(s, f)
a = s;
if isstruct(s) && isfield(s, f), a = rmfield(s, f); end
end

function s = run_status_of(mask)
if any(mask(:)), s = 'complete'; else, s = 'empty'; end
end

function sib = siblings(out_dir, P)
% grid, Lh, res and model_version of branch_HSS_latest.mat, checked
% bit-equal to the HPSS and MPSS siblings; par equal to P in all three
flds = {'dn_R', 'dn_B', 'R_DIFF', 'B_DIFF', 'Lh', 'res', 'model_version', 'par'};
sib = [];
for ty = {'HSS', 'HPSS', 'MPSS'}
    f = fullfile(out_dir, sprintf('branch_%s_latest.mat', ty{1}));
    if exist(f, 'file') ~= 2
        error('GNX:MODE', 'Sibling library file missing: %s.', f);
    end
    S = load(f, flds{:});
    if ~isequal(S.par, P)
        error('GNX:MODE', '%s was built at other parameters than P.', f);
    end
    if isempty(sib)
        sib = S;
    elseif ~isequal(rmfield(S, 'par'), rmfield(sib, 'par'))
        error('GNX:MODE', '%s: grid, Lh, res or model_version differ from branch_HSS_latest.mat.', f);
    end
end
if ~isequal(sib.Lh, P.L / 2) || ~isequal(sib.res, P.res)
    error('GNX:MODE', 'Sibling Lh/res differ from P.');
end
end

function acc = note_acceptance(acc, rec, i, j)
c = rec.checks;
if isfield(c, 'seam_slope'), acc.seam_slope(i, j) = c.seam_slope.value; end
if isfield(c, 'trough_slope'), acc.trough_slope(i, j) = c.trough_slope.value; end
if isfield(c, 'symmetry'), acc.symmetry(i, j) = c.symmetry.value; end
if isfield(c, 'mass_doubling'), acc.mass_rel(i, j, :) = c.mass_doubling.value; end
if isfield(c, 'peak_fraction'), acc.peak_position(i, j) = c.peak_fraction.value; end
if isfield(c, 'contrast'), acc.contrast(i, j) = c.contrast.value; end
end

function why = first_failed(rec, n)
if strcmp(rec.kind, 'mpss'), names = gnx_validity('required', 'mpss'); else, names = required_mode(n); end
why = 'unknown';
for k = 1:numel(names)
    if ~isfield(rec.checks, names{k}) || ~rec.checks.(names{k}).pass
        why = names{k}; return;
    end
end
end

% -------------------------------------------------------------------------
function names = required_mode(n)
if n >= 3
    names = {'quarter_valid', 'seam_slope', 'trough_slope', 'symmetry', 'mass_doubling', ...
             'peak_fraction', 'contrast'};
else
    names = {'quarter_valid', 'seam_slope', 'symmetry', 'mass_doubling', ...
             'peak_fraction', 'contrast'};
end
end

function ok = verdict_any(rec)
if isstruct(rec) && isfield(rec, 'kind') && strcmp(rec.kind, 'mode')
    ok = false;
    if ~isfield(rec, 'checks') || ~isfield(rec, 'n'), return; end
    names = required_mode(rec.n);
    for k = 1:numel(names)
        if ~isfield(rec.checks, names{k}), return; end
        c = rec.checks.(names{k});
        if ~isstruct(c) || ~isfield(c, 'pass') || ~isequal(c.pass, true), return; end
    end
    ok = true;
else
    ok = gnx_validity('verdict', rec);
end
end

function ok = verdict(rec)
% the n = 2 record's local verdict (gnx_validity's own)
ok = gnx_validity('verdict', rec);
end

function c = chk(pass, value, tol, varargin)
% copied from gnx_validity.m
c = struct('pass', logical(pass), 'value', value, 'tol', tol);
for k = 1:2:numel(varargin), c.(varargin{k}) = varargin{k+1}; end
end

% -------------------------------------------------------------------------
function rec = legacy_record(qrec, solq, ctx)
rec = mpss_record(qrec, solq, ctx);
rec.valid = verdict(rec);                 % as in gnx_validity.m
end

function rec = mpss_record(qrec, solq, ctx)
% Copy of mpss_record in gnx_validity.m (unchanged).
D = gnx_model('descriptor'); tol = D.tol; nm = D.numerics;
rec = struct('kind', 'mpss', 'checks', struct(), 'diag', struct(), 'valid', false, ...
    'x', [], 'y', []);
qok = isstruct(qrec) && isfield(qrec, 'kind') && strcmp(qrec.kind, 'quarter') && verdict(qrec);
rec.checks.quarter_valid = chk(qok, double(qok), NaN);
if ~qok, return; end
xq = solq.x(:).'; ellq = xq(end);
yq = solq.y([1 3 5 7], :);
% peak-to-peak reflection (same convention as build_mpss_branch.m): the
% quarter's peak (x = 0) goes to the seam x = ell_q of the stitched half.
x_left = ellq - fliplr(xq);  x_right = ellq + xq;
xs = [x_left(1:end-1), x_right];
ys = [fliplr(yq(:, 2:end)), yq];
rec.x = xs; rec.y = ys;
seam = abs(solq.y(6, 1));
rec.checks.seam_slope = chk(seam <= tol.mpss_slope, seam, tol.mpss_slope);
sym = max(abs(ys - fliplr(ys)), [], 'all');
rec.checks.symmetry = chk(sym <= tol.mpss_sym, sym, tol.mpss_sym);
% dense mass of the stitched profile from the quarter's Hermite interpolant
xd = linspace(0, 2 * ellq, 2 * nm.n_dense - 1);
xm = abs(xd - ellq);
xm(xm > ellq) = ellq;
Yd = deval(solq, xm);
mR = trapz(xd, Yd(1, :) + Yd(3, :));
mB = trapz(xd, Yd(5, :) + Yd(7, :));
eR = abs(mR - ctx.nR_half) / ctx.nR_half; eB = abs(mB - ctx.nB_half) / ctx.nB_half;
rec.checks.mass_doubling = chk(max(eR, eB) <= tol.mass_rel, [eR eB], tol.mass_rel, ...
    'dense_mass', [mR mB]);
[bmax, ipk] = max(ys(3, :));
pf = xs(ipk) / xs(end);
rec.checks.peak_fraction = chk(pf > tol.mpss_peak_lo && pf < tol.mpss_peak_hi, pf, ...
    [tol.mpss_peak_lo tol.mpss_peak_hi]);
ctr = (bmax - min(ys(3, :))) / max(bmax, 1e-12);
rec.checks.contrast = chk(ctr >= tol.contrast, ctr, tol.contrast);
rec.diag.outer_trough = ys(:, end).';
end

% -------------------------------------------------------------------------
function rec = mode_record(n, qrec, solq, ctx)
% the n-general record (see the header)
D = gnx_model('descriptor'); tol = D.tol; nm = D.numerics;
rec = struct('kind', 'mode', 'n', n, 'checks', struct(), 'diag', struct(), 'valid', false, ...
    'x', [], 'y', []);
qok = isstruct(qrec) && isfield(qrec, 'kind') && strcmp(qrec.kind, 'quarter') && verdict(qrec);
rec.checks.quarter_valid = chk(qok, double(qok), NaN);
if ~qok, return; end
xq = solq.x(:).'; ell = xq(end);
yq = solq.y([1 3 5 7], :);
[xs, ys] = stitch(xq, yq, n, 'T');          % T-start: a trough at the pole
rec.x = xs; rec.y = ys;
seam = abs(solq.y(6, 1));                   % the piece's peak end
rec.checks.seam_slope = chk(seam <= tol.mpss_slope, seam, tol.mpss_slope);
if n >= 3                                   % the trough end is an interior join
    tro = abs(solq.y(6, end));
    rec.checks.trough_slope = chk(tro <= tol.mpss_slope, tro, tol.mpss_slope);
end
if mod(n, 2) == 0
    sym = max(abs(ys - fliplr(ys)), [], 'all');
    how = 'T-start stitch against its reflection about the half centre';
else
    [~, ysP] = stitch(xq, yq, n, 'P');
    sym = max(abs(ys - fliplr(ysP)), [], 'all');
    how = 'T-start stitch against the reflected P-start stitch of the same piece';
end
rec.checks.symmetry = chk(sym <= tol.mpss_sym, sym, tol.mpss_sym, 'comparison', how);
npts = n * (nm.n_dense - 1) + 1;
xd = linspace(0, n * ell, npts);
xm = fold(xd, ell, 'T');
xm(xm > ell) = ell;
Yd = deval(solq, xm);
mR = trapz(xd, Yd(1, :) + Yd(3, :));
mB = trapz(xd, Yd(5, :) + Yd(7, :));
eR = abs(mR - ctx.nR_half) / ctx.nR_half; eB = abs(mB - ctx.nB_half) / ctx.nB_half;
rec.checks.mass_doubling = chk(max(eR, eB) <= tol.mass_rel, [eR eB], tol.mass_rel, ...
    'dense_mass', [mR mB], 'n_points', npts);
bq = yq(3, :);
kq = find(bq == max(bq), 1, 'last');
pos = xq(kq) / ell;
rec.checks.peak_fraction = chk(pos < 1 - 2 * tol.mpss_peak_lo, pos, 1 - 2 * tol.mpss_peak_lo, ...
    'definition', 'piece argmax position x*/l (last maximal node) < 1 - 2 mpss_peak_lo');
bmax = max(ys(3, :));
ctr = (bmax - min(ys(3, :))) / max(bmax, 1e-12);
rec.checks.contrast = chk(ctr >= tol.contrast, ctr, tol.contrast);
rec.diag.seam_end = ys(:, end).';
rec.valid = verdict_any(rec);
end

% -------------------------------------------------------------------------
function [x, y] = stitch(xq, yq, n, start)
% identical to sap_seed_shape / sag_seed_shape 'stitch'
xq = xq(:).';
l = xq(end);
x = zeros(1, 0); y = zeros(size(yq, 1), 0);
for k = 0:n-1
    if strcmp(start, 'T'), mirrored = mod(k, 2) == 0; else, mirrored = mod(k, 2) == 1; end
    if mirrored
        xk = (k + 1) * l - fliplr(xq); yk = fliplr(yq);
    else
        xk = k * l + xq; yk = yq;
    end
    if k < n - 1
        xk = xk(1:end-1); yk = yk(:, 1:end-1);
    end
    x = [x, xk]; %#ok<AGROW>
    y = [y, yk]; %#ok<AGROW>
end
end

function xm = fold(x, l, start)
% identical to sag_seed_shape 'fold'
r = mod(x, 2 * l);
if strcmp(start, 'T')
    xm = abs(r - l);
else
    xm = min(r, 2 * l - r);
end
end

% -------------------------------------------------------------------------
function out_path = write_mode(out_dir, n, P, data, meta)
D = gnx_model('descriptor');
here = fileparts(mfilename('fullpath'));
sib = siblings(out_dir, P);
res = P.res;
Lh = P.L / 2;
ell = P.L / (2 * n);
L = struct();
L.model_id = D.id;
L.model_version = sib.model_version;           % the siblings' (the reader requires one value)
L.builder_model_version = gnx_model('version', {here});   % this folder now (provenance)
L.schema = D.schema;
L.fields = D.fields;
L.state_layout = D.state_layout;
L.branch_type = sprintf('MODE%d', n);
L.dn_R = sib.dn_R;                             % copied from branch_HSS_latest.mat
L.dn_B = sib.dn_B;
L.R_DIFF = sib.R_DIFF;
L.B_DIFF = sib.B_DIFF;
L.Lh = sib.Lh;
L.res = sib.res;
L.par = P;
L.existence_mask = logical(data.existence_mask);
L.SOL = data.SOL;
L.validity = data.validity;
L.flux_values = nan(res, res, 4);
L.flux_fields = D.fields;
L.flux_source_location = strings(res, res);
L.mean_cytosol = nan(res, res, 2);
L.mean_cytosol_fields = {'Rc', 'Bc'};
L.is_stable = false(res, res);
L.is_stable_computed = false(res, res);
L.stability_source = 'not_computed';
L.run_id = meta.run_id;
L.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
L.provenance = meta;
L.orientation = 'trough_at_0';                 % the T-start stitch (drawn n-)
L.domain = [0 Lh; 0 ell];
L.storage_kind = 'piece_bvp4c+stitched_xy';
L.construction = 'reflected_piece_HPSS_n';
L.mode_n = n;
L.piece_length = ell;
L.SOL_piece = data.SOL_piece;                  % bvp4c struct, peak at x = 0
L.fail_reason = data.fail_reason;
L.piece_ok_grid = data.piece.ok_grid;
L.piece_last_reason = data.piece.last_reason;
L.piece_validity = data.piece.validity;
L.phase_log = data.piece.phase_log;
L.n_calls = data.n_calls; L.n_nmax = data.n_nmax;
L.cap_reached = data.piece.cap_reached;
L.run_status = run_status_of(L.existence_mask);
L.acceptance = data.acceptance;
L.record_branch = data.record_branch;
L.source_sweep_file = data.source_sweep_file;
L.source_sweep_sha256 = sha256_file(data.source_sweep_file);
L.builder_sha256 = sha256_file(fullfile(here, 'gnx_build_mode.m'));
for i = 1:res
    for j = 1:res
        if ~L.existence_mask(i, j), continue; end
        sq = L.SOL_piece{i, j};
        check_struct(sq, (P.rho_R + L.dn_R(j) / 2) * ell, (P.rho_B + L.dn_B(i) / 2) * ell, D, i, j);
        L.flux_values(i, j, :) = L.SOL{i, j}.y(:, end);
        L.flux_source_location(i, j) = "seam_x_Lh";
        v = data.piece.validity{i, j};
        if isstruct(v) && isfield(v, 'diag') && isfield(v.diag, 'cyt_mean_121')
            L.mean_cytosol(i, j, :) = v.diag.cyt_mean_121;
        end
    end
end
ex = L.existence_mask;
nonempty = ~cellfun(@isempty, L.SOL);
fv = reshape(L.flux_values, res * res, 4);
finite_flux = reshape(all(isfinite(fv), 2), res, res);
if ~isequal(ex, nonempty) || ~isequal(ex, finite_flux) || ~isequal(ex, ~cellfun(@isempty, L.SOL_piece))
    error('GNX:MODE', 'MODE%d: existence, SOL, SOL_piece and flux disagree.', n);
end
for k = find(ex)'
    if ~verdict_any(L.validity{k})
        error('GNX:MODE', 'MODE%d: stored cell %d has a failing record.', n, k);
    end
end
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
out_path = fullfile(out_dir, sprintf('branch_MODE%d_%s.mat', n, meta.run_id));
save(out_path, '-struct', 'L');
copyfile(out_path, fullfile(out_dir, sprintf('branch_MODE%d_latest.mat', n)));
end

function check_struct(s, nR, nB, D, i, j)
% as in gnx_write_library.m, on the piece
if ~(isstruct(s) && isfield(s, 'y') && size(s.y, 1) == 12)
    error('GNX:MODE', 'cell (%d,%d): stored piece is not a 12-row BVP struct.', i, j);
end
a0 = s.y([11 12], 1);
aend = s.y([11 12], end);
rel = abs(aend - [nR; nB]) ./ [nR; nB];
if any(a0 ~= 0) || any(aend == 0) || any(rel > D.tol.acc_end_rel)
    error('GNX:MODE', ['cell (%d,%d): piece accumulators start at [%g %g] and end at ' ...
        '[%g %g] (targets %g %g); refusing a zeroed or off-mass struct.'], ...
        i, j, a0(1), a0(2), aend(1), aend(2), nR, nB);
end
end

function h = sha256_file(p)
fid = fopen(p, 'r');
if fid < 0, error('GNX:MODE', 'Cannot open %s.', p); end
c = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes(:));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2)', 1, []));
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
