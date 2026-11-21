function rows = sa_t42_am8(E)
%SA_T42_AM8  The starting-profile mesh rules
%   continuation.guess_mesh / continuation.guess_nodes (sa_guess_mesh).
%   On the stub search with the wall rule (continuation edges, substeps, R1
%   diagonal edges and R2 secant calls all occur), every guess the corrector
%   receives is logged and checked bit for bit against the rule:
%     default 'parent' (parent iterate on its own mesh;
%       secant 2*deval(P1, P1.x) - deval(P2, P1.x)); guess_source unchanged;
%     'uniform' N: x = linspace(P.x(1), P.x(end), N), y = deval(P, x);
%     'subsample' N: at most N of the parent's own nodes, y = P.y at them;
%       the parent mesh itself when it has <= N nodes;
%   determinism (two runs, identical guesses and records); bounded guess
%   size; resume equivalence (kill mid-call + resume == uninterrupted) with
%   each new rule on; the leaf is scientific (resume with a changed rule is
%   refused); multipoint iterates refused under the new rules. Zero solves.
%   Also: myxo_restart under 'subsample'
%   and 'uniform' (each level's first guess = the rule on the parent, each
%   within-level guess = the rule on the last valid substep, bit for bit);
%   the R3 deep rung occurs and its guesses follow the rule; boundary sizes
%   n in {N-1, N, N+1, 2N-1, 2N} of sa_guess_mesh('mesh'); 'uniform'
%   enlarges a parent with n < N, 'subsample' never does; guess_nodes is
%   inert under 'parent' but scientific (a resume with it changed is
%   refused); the report-only guess node count (attempts.csv / attempts.mat
%   column guess_nodes) equals the guess each call received, and replays
%   identically on resume.
T = 'T42 AM-8 guess mesh';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't42');
fail = @(dd, b, ex) sa_m3_wall(dd, ex);
box = [3 1];
M = sa_model('stub_two_pool');

% ---- settings ------------------------------------------------------------------
S0 = sa_settings('resolve', sa_m3_stub('settings', fullfile(d, 'x')), M);
ok = strcmp(S0.continuation.guess_mesh, 'parent') && S0.continuation.guess_nodes == 1000;
bad = {struct('guess_mesh', 'cold'), struct('guess_nodes', 2), struct('guess_nodes', 10.5)};
nref = 0;
for k = 1:numel(bad)
    s = sa_m3_stub('settings', fullfile(d, 'x')); s.continuation = bad{k};
    try, sa_settings('resolve', s, M); catch ME, nref = nref + strcmp(ME.identifier, 'SA:SETTINGS:E_FIELD'); end
end
rows = sa_check(rows, T, 'defaults guess_mesh = ''parent'', guess_nodes = 1000; unknown rule / N < 3 / non-integer N refused', ...
    ok && nref == 3, sprintf('%d of 3 refused', nref));

% ---- default 'parent' rule, bit for bit ---------------------------------------
[r0, st0] = sa_m3_run(fullfile(d, 'parent'), box, 'fail', fail, 'log', {});
[okp, info] = check_guesses(r0, st0('log'), 'parent', Inf);
rows = sa_check(rows, T, ['default: every non-cold guess = v0.3 (struct(x = P.x, y = P.y); secant ' ...
    '2*deval(P1, P1.x) - deval(P2, P1.x)), bitwise; guess_source without a tag'], okp, info);

% ---- 'uniform' and 'subsample' -----------------------------------------------------
N = 50;
ou = struct('continuation', struct('guess_mesh', 'uniform', 'guess_nodes', N));
os = struct('continuation', struct('guess_mesh', 'subsample', 'guess_nodes', N));
[ru, stu] = sa_m3_run(fullfile(d, 'uniform'), box, 'fail', fail, 'log', {}, 'settings', ou);
[oku, info] = check_guesses(ru, stu('log'), 'uniform', N);
rows = sa_check(rows, T, sprintf(['uniform N = %d: every non-cold guess has exactly N nodes, ' ...
    'x = linspace over the parent range, y = deval(parent, x) (secant on that mesh), bitwise; tag in guess_source'], N), ...
    oku, info);
[rs, sts] = sa_m3_run(fullfile(d, 'subsample'), box, 'fail', fail, 'log', {}, 'settings', os);
[oks, info] = check_guesses(rs, sts('log'), 'subsample', N);
rows = sa_check(rows, T, sprintf(['subsample N = %d (parents 201 nodes): every non-cold guess = the ' ...
    'parent''s nodes round(linspace(1, n, N)), y = parent y there, bitwise; <= N nodes'], N), oks, info);
ob = struct('continuation', struct('guess_mesh', 'subsample', 'guess_nodes', 500));
[rb, stb] = sa_m3_run(fullfile(d, 'subsample500'), box, 'fail', fail, 'log', {}, 'settings', ob);
[okb, info] = check_guesses(rb, stb('log'), 'parent_mesh_tagged', 500);
rows = sa_check(rows, T, 'subsample N = 500 >= parent nodes (201): guess = the parent mesh and values, bitwise', ...
    okb, info);
% ---- boundary sizes of sa_guess_mesh('mesh') ---------------------------
okz = true; det = {};
for n = [N - 1, N, N + 1, 2 * N - 1, 2 * N]
    P = boundary_sol(n);
    ss = S0; ss.continuation.guess_mesh = 'subsample'; ss.continuation.guess_nodes = N;
    x = sa_guess_mesh('mesh', P, ss);
    g = sa_guess_mesh('from', P, ss);
    okn = numel(x) == min(n, N) && all(diff(x) > 0) && x(1) == P.x(1) && x(end) == P.x(end) && ...
        all(ismember(x, P.x)) && isequal(g.x, x) && size(g.y, 2) == numel(x);
    su = S0; su.continuation.guess_mesh = 'uniform'; su.continuation.guess_nodes = N;
    xu = sa_guess_mesh('mesh', P, su);
    okn = okn && numel(xu) == N && all(diff(xu) > 0) && xu(1) == P.x(1) && xu(end) == P.x(end);
    okz = okz && okn;
    det{end+1} = sprintf('n=%d: subsample %d, uniform %d', n, numel(x), numel(xu)); %#ok<AGROW>
end
rows = sa_check(rows, T, sprintf(['boundary sizes (N = %d): subsample gives exactly min(n, N) strictly ' ...
    'increasing parent nodes with both ends kept for n in {N-1, N, N+1, 2N-1, 2N}; uniform gives N ' ...
    '(ENLARGES n = N-1 to N), subsample never enlarges'], N), okz, strjoin(det, '; '));

% ---- report-only guess node count ------------------------------------------------
okg = true; det = {};
for c = {{r0, st0, 'parent'}, {rs, sts, 'subsample'}, {ru, stu, 'uniform'}}
    rr = c{1}{1}; G = c{1}{2}('log');
    want = cellfun(@(g) numel(g.x), G);
    Z = load(fullfile(sa_layout('state', rr.run_dir), 'attempts.mat'), 'guess_nodes');
    gcsv = csv_guess_nodes(fullfile(rr.run_dir, 'attempts.csv'));
    okg = okg && isequal(rr.guess_nodes, want) && isequal(Z.guess_nodes, want) && isequal(gcsv, want);
    det{end+1} = sprintf('%s: %d rows, guess nodes %d..%d', c{1}{3}, numel(want), min(want), max(want)); %#ok<AGROW>
end
rows = sa_check(rows, T, ['R-8.2: attempts.csv / attempts.mat column guess_nodes (report-only, not a ' ...
    'row field) = node count of the guess each call received, every call incl. cold'], okg, strjoin(det, '; '));

same = isequal({r0.attempts.outcome}, {ru.attempts.outcome}, {rs.attempts.outcome}) && ...
    isequal(strip_tag({ru.attempts.unit_key}), {r0.attempts.unit_key}) && ...
    isequal(strip_tag({rs.attempts.unit_key}), {r0.attempts.unit_key});
rows = sa_check(rows, T, ['stub (guess-mesh blind): the three rules give the same work units (up to the ' ...
    'guess_source tag) and outcomes'], same, sprintf('%d calls each', r0.calls));

% ---- determinism --------------------------------------------------------------------------
[rs2, sts2] = sa_m3_run(fullfile(d, 'subsample_2'), box, 'fail', fail, 'log', {}, 'settings', os);
[ru2, stu2] = sa_m3_run(fullfile(d, 'uniform_2'), box, 'fail', fail, 'log', {}, 'settings', ou);
det = isequal(sts('log'), sts2('log')) && isequal(stu('log'), stu2('log')) && ...
    isequal({rs.attempts.unit_key}, {rs2.attempts.unit_key}) && isequal({ru.attempts.unit_key}, {ru2.attempts.unit_key});
rows = sa_check(rows, T, 'determinism: a second run of each rule receives bit-identical guesses and makes the same units', ...
    det, sprintf('%d guesses', numel(sts('log'))));

% ---- resume equivalence with each new rule on -------------------------------------------
ok = true; det = {};
for rule = {'uniform', 'subsample'}
    o = struct('continuation', struct('guess_mesh', rule{1}, 'guess_nodes', N));
    U = fullfile(d, ['U_' rule{1}]);
    rU = sa_m3_run(U, box, 'fail', fail, 'settings', o);
    n = rU.calls;
    for k = unique([3, round(n / 2), n - 1])
        K = fullfile(d, sprintf('K_%s_%03d', rule{1}, k));
        try
            sa_m3_run(K, box, 'fail', fail, 'kill_at', k, 'settings', o);
        catch ME
            assert(strcmp(ME.identifier, 'SA:TEST:KILL'), ME.message);
        end
        oo = o; oo.resume = true;
        r = sa_m3_run(K, box, 'fail', fail, 'settings', oo);
        D = sa_m3_equiv(U, K);
        gk = r.guess_nodes; gu = rU.guess_nodes;
        okn = numel(gk) == numel(gu) + 1 && isequal(gk([1:k-1, k+1:end]), gu) && gk(k) == gk(k + 1);
        ok = ok && D.pass && D.n_interrupted == 1 && isequal(r.replay_interrupted, k) && okn;
        det{end+1} = sprintf('%s k=%d %d', rule{1}, k, D.pass); %#ok<AGROW>
    end
end
rows = sa_check(rows, T, ['resume equivalence: kill mid-call (3 points per rule) + resume == uninterrupted, ' ...
    'uniform and subsample (and the replayed guess_nodes column equals the uninterrupted one)'], ...
    ok, strjoin(det, ' '));

% ---- the leaf is scientific: resume with another rule is refused ----------------------
C = fullfile(d, 'copy_uniform');
copyfile(fullfile(d, 'U_uniform'), C);
o = os; o.resume = true;
code = '';
try
    sa_m3_run(C, box, 'fail', fail, 'settings', o);
catch ME
    code = ME.identifier;
end
rows = sa_check(rows, T, 'resume of a uniform run with guess_mesh = subsample is refused (scientific leaf)', ...
    strcmp(code, 'SA:RESUME:E_SETTINGS'), code);

% ---- guess_nodes is inert under 'parent' but scientific --------------------------
C = fullfile(d, 'copy_parent');
copyfile(fullfile(d, 'parent'), C);
o = struct('continuation', struct('guess_nodes', 500), 'resume', true);
code = ''; msg = '';
try
    sa_m3_run(C, box, 'fail', fail, 'settings', o);
catch ME
    code = ME.identifier; msg = ME.message;
end
rows = sa_check(rows, T, ['resume of a parent run with only guess_nodes changed (1000 -> 500, inert under ' ...
    '''parent'') is refused (scientific leaf)'], strcmp(code, 'SA:RESUME:E_SETTINGS') && ...
    contains(msg, 'guess_nodes'), code);

% ---- myxo_restart levels and the R3 deep rung under the new rules -----------------
for rule = {'subsample', 'uniform'}
    o = struct('continuation', struct('step_rule', 'myxo_restart', 'guess_mesh', rule{1}, 'guess_nodes', N), ...
        'rescue', struct('rungs', {{'diagonal', 'deep'}}));
    [rm, stm] = sa_m3_run(fullfile(d, ['myxo_' rule{1}]), box, 'fail', fail, 'log', {}, 'settings', o);
    [okm, info] = check_guesses(rm, stm('log'), rule{1}, N, ...
        {'continuation', 'substep', 'rescue_diagonal', 'rescue_deep'});
    [okl, linfo] = check_myxo_levels(rm);
    rows = sa_check(rows, T, sprintf(['myxo_restart, %s N = %d, rungs diagonal + deep: each level''s first ' ...
        'guess = the rule on the edge parent, each within-level guess = the rule on the last valid ' ...
        'substep (bitwise); R3 deep calls occur and follow the rule'], rule{1}, N), okm && okl, ...
        [info '; ' linfo]);
end

% ---- multipoint refusal -----------------------------------------------------------------------
th = sa_m3_stub('params');
P = sa_m3_stub('sol', th, [2 2 2 2], 'A');
P.x(100) = P.x(101);
nmp = 0;
for rule = {'uniform', 'subsample'}
    s = S0; s.continuation.guess_mesh = rule{1}; s.continuation.guess_nodes = N;
    try, sa_guess_mesh('from', P, s); catch ME, nmp = nmp + strcmp(ME.identifier, 'SA:GUESSMESH:E_MULTIPOINT'); end
end
g = sa_guess_mesh('from', P, S0);
rows = sa_check(rows, T, 'a multipoint iterate (repeated x) is refused under uniform and subsample; parent passes it through', ...
    nmp == 2 && isequal(g, struct('x', P.x, 'y', P.y)), sprintf('%d of 2 refused', nmp));
rmpath(E.stub);
end

% =========================================================================

function [ok, info] = check_guesses(res, G, rule, N, need)
if nargin < 5, need = {'continuation', 'substep', 'rescue_diagonal', 'rescue_secant'}; end
A = res.attempts;
ok = numel(G) == numel(A);
bad = {}; nk = struct('continuation', 0, 'substep', 0, 'rescue_diagonal', 0, 'rescue_secant', 0, ...
    'rescue_deep', 0);
maxn = 0;
for k = 1:min(numel(G), numel(A))
    a = A(k);
    if strcmp(a.kind, 'cold'), continue; end
    g = G{k};
    if strcmp(a.route, 'rescue_secant')
        t = regexp(a.guess_source, '^secant 2\*(A\d+) - (A\d+)', 'tokens', 'once');
        P1 = load_sol(res, t{1}); P2 = load_sol(res, t{2});
        x1 = expected_mesh(P1, rule, N);
        e = struct('x', x1, 'y', 2 * deval(P1, x1) - deval(P2, x1));
        nk.rescue_secant = nk.rescue_secant + 1;
    else
        P = load_sol(res, a.parent_attempt);
        [x, idx] = expected_mesh(P, rule, N);
        if isempty(idx), y = deval(P, x); else, y = P.y(:, idx); end
        e = struct('x', x, 'y', y);
        if a.fraction < 1, nk.substep = nk.substep + 1;
        elseif strcmp(a.route, 'continuation'), nk.continuation = nk.continuation + 1;
        elseif strcmp(a.route, 'rescue_diagonal'), nk.rescue_diagonal = nk.rescue_diagonal + 1;
        elseif strcmp(a.route, 'rescue_deep'), nk.rescue_deep = nk.rescue_deep + 1;
        end
    end
    want_tag = '';
    if ~strcmp(rule, 'parent')
        r = rule; if strcmp(r, 'parent_mesh_tagged'), r = 'subsample'; end
        want_tag = sprintf(' [%s %d]', r, N);
    end
    if ~isequal(g, e) || ~endsWith(a.guess_source, want_tag) || ...
            (isempty(want_tag) && contains(a.guess_source, '['))
        bad{end+1} = a.attempt_id; %#ok<AGROW>
    end
    maxn = max(maxn, numel(g.x));
end
ok = ok && isempty(bad) && all(cellfun(@(f) nk.(f) > 0, need));
if any(strcmp(rule, {'uniform', 'subsample'})), ok = ok && maxn <= N; end
info = sprintf('%d calls: continuation %d, substeps %d, R1 %d, R2 %d, R3 %d; max guess nodes %d; mismatches %s', ...
    numel(A), nk.continuation, nk.substep, nk.rescue_diagonal, nk.rescue_secant, nk.rescue_deep, maxn, ...
    strjoin(bad, ','));
end

function [ok, info] = check_myxo_levels(res)
% Every myxo_restart edge: the first call and every call after a refused one
% (a level restart) are built from the edge parent; every call after a valid
% substep (within a level) from that substep. Together with check_guesses
% (guess == rule applied to the attempt named parent_attempt) this fixes the rule.
A = res.attempts; ids = {A.attempt_id};
E = res.Q.edges(strcmp({res.Q.edges.rule}, 'myxo_restart'));
ok = ~isempty(E); nlev = 0; nin = 0; bad = {};
for ei = 1:numel(E)
    e = E(ei);
    prev = '';
    for j = 1:numel(e.attempts)
        a = A(strcmp(ids, e.attempts{j}));
        if isempty(prev) || ~A(strcmp(ids, prev)).valid
            want = e.parent_attempt;
            if ~isempty(prev), nlev = nlev + 1; end
        else
            want = prev; nin = nin + 1;
        end
        if ~strcmp(a.parent_attempt, want), bad{end+1} = a.attempt_id; end %#ok<AGROW>
        prev = a.attempt_id;
    end
end
ok = ok && isempty(bad) && nlev > 0 && nin > 0;
info = sprintf('myxo_restart edges %d, level restarts %d, within-level guesses %d, mismatches %s', ...
    numel(E), nlev, nin, strjoin(bad, ','));
end

function P = boundary_sol(n)
% a two-point iterate with a non-uniform, strictly increasing mesh of n nodes
x = linspace(0, 1.5, n) .^ 2 / 1.5;
y = [sin(x); cos(x); x];
P = struct('solver', 'bvp4c', 'x', x, 'y', y, 'yp', [cos(x); -sin(x); ones(size(x))], ...
    'stats', struct('maxres', 1e-9, 'nmeshpoints', n));
end

function g = csv_guess_nodes(f)
L = strsplit(strtrim(fileread(f)), newline);
h = strsplit(L{1}, ',');
% a quoted plain-words reason column follows guess_nodes
assert(strcmp(h{end}, 'reason') && strcmp(h{end-1}, 'guess_nodes'), 'attempts.csv: guess_nodes is not the column before reason');
L(2:end) = regexprep(L(2:end), ',"(?:[^"]|"")*"$', '');
g = cellfun(@(r) str2double(r(find(r == ',', 1, 'last') + 1:end)), L(2:end));
end

function [x, idx] = expected_mesh(P, rule, N)
n = numel(P.x);
switch rule
    case {'parent', 'parent_mesh_tagged'}
        x = P.x; idx = 1:n;
    case 'uniform'
        x = linspace(P.x(1), P.x(end), N); idx = [];
    case 'subsample'
        if n <= N, idx = 1:n; else, idx = round(linspace(1, n, N)); end
        x = P.x(idx);
end
end

function s = load_sol(res, id)
Z = load(fullfile(sa_layout('state', res.run_dir), 'attempts', [id '.mat']), 'sol');
s = Z.sol;
end

function c = strip_tag(c)
c = regexprep(c, ' \[(uniform|subsample) \d+\]', '');
end
