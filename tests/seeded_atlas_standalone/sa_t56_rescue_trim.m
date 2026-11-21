function rows = sa_t56_rescue_trim(E)
%SA_T56_RESCUE_TRIM  Rescue
%   trimmed to secant. (1) The mode defaults are secant only, no ring
%   rescue; a stub run under them makes only secant rescue edges, at most
%   2 x n_pools = 4 calls per node, none into the ring. (2) a node
%   whose ladder makes no call (no lineage pair at q-e, q-2e) records
%   rescue_result 'no_basis', is NOT marked rescue_done and makes 0 rescue
%   calls; once a lineage pair exists (a later closure pass, after another
%   node was rescued and joined the lineage), secant fires exactly once there
%   and rescues it. (3) With the explicit full ladder, sa_t27_rescue, sa_t48_am10
%   and sa_t49_am11_rescue_ring pass unchanged. Zero solves (stub corrector).
T = 'T56 rescue trim';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't56');
% ---- (1) defaults and a secant-only run ---------------------------------------------------
Mf = sa_model('original_polarity'); Mn = sa_model('guzzo_noxa');
ok = true;
for mm = {Mf, Mn}
    for m = sa_modes('modes')
        u = struct('model', mm{1}.id, 'run_dir', 'x', 'mode', m{1});
        if any(strcmp(m{1}, {'production', 'focus'})), u.budget = struct('max_calls', 10, 'max_wall_s', 10); end
        if strcmp(m{1}, 'focus'), u.import = struct('source_dir', 'x'); end   % focus mode needs a source
        S = sa_settings('resolve', u, mm{1});
        ok = ok && isequal(S.rescue.rungs, {'secant'}) && ~S.rescue.ring && S.routes.rescue;
    end
end
trim = struct('rescue', struct('rungs', {{'secant'}}, 'ring', false));
r1 = sa_m3_run(fullfile(d, 'trim'), [], 'fail', @(dd, b, ex) walls(dd, ex, 'ringwall'), 'settings', trim);
Er = r1.Q.edges(startsWith({r1.Q.edges.kind}, 'rescue'));
rc = r1.Q.rescue_calls(r1.Q.rescue_calls > 0);
rows = sa_check(rows, T, ['defaults (both models, every mode): rescue on, rungs {secant}, ring false; a stub run under ' ...
    'them: only rescue_secant edges, <= 2 x n_pools = 4 rescue calls per node, no rescue edge into the ring'], ...
    ok && ~isempty(Er) && all(strcmp({Er.kind}, 'rescue_secant')) && all(rc <= 4) && ...
    ~any(r1.Q.ring([Er.q])) && ~any(r1.Q.rescue_done & r1.Q.ring), ...
    sprintf('%d rescue edges (%s), max calls per node %d', numel(Er), strjoin(unique({Er.kind}), ','), max([0; rc])));
% ---- (2) no_basis ------------------------------------------------------------------------------
% walls at k:-1,0 and k:-2,0 refuse every non-exact guess. Pass 1: k:-2,0 has
% cardinal roots (k:-3,0, k:-2,+-1) but no lineage pair in any direction ->
% no_basis, 0 calls; k:-1,0 is rescued by secant from k:0,0 and k:1,0 and
% joins the anchor lineage. Pass 2: k:-2,0 now has the pair (k:-1,0, k:0,0)
% -> secant fires once and rescues it.
o = trim; o.continuation = struct('box', [3 1]);
r2 = sa_m3_run(fullfile(d, 'rc52'), [3 1], 'fail', @(dd, b, ex) walls(dd, ex, 'pair'), 'settings', o);
q2 = find(all(r2.L.index == [-2 0], 2)); q1 = find(all(r2.L.index == [-1 0], 2));
E2 = r2.Q.edges;
into2 = E2([E2.q] == q2 & startsWith({E2.kind}, 'rescue'));
rows = sa_check(rows, T, ['the no-basis node is retried in the closure pass: secant fires exactly once at ' ...
    'k:-2,0 (search pass 2) and rescues it; k:-1,0 rescued in pass 1'], ...
    numel(into2) == 1 && strcmp(into2.kind, 'rescue_secant') && into2.calls == 1 && into2.pass == 2 && ...
    strcmp(into2.status, 'completed') && strcmp(r2.Q.rescue_result{q2}, 'rescued_secant') && ...
    r2.Q.rescue_done(q2) && strcmp(r2.Q.rescue_result{q1}, 'rescued_secant') && r2.Q.rescue_calls(q2) == 1, ...
    sprintf('k:-2,0 %s (%d calls, pass %s); k:-1,0 %s', r2.Q.rescue_result{q2}, r2.Q.rescue_calls(q2), ...
    mat2str([into2.pass]), r2.Q.rescue_result{q1}));
o1 = o; o1.search = struct('max_passes', 1);
r3 = sa_m3_run(fullfile(d, 'rc52_one_pass'), [3 1], 'fail', @(dd, b, ex) walls(dd, ex, 'pair'), 'settings', o1);
E3 = r3.Q.edges;
rows = sa_check(rows, T, ['with one search pass, k:-2,0 (a cold-free lineage-pair-free ladder) makes 0 calls, ' ...
    'records no_basis and is not rescue_done; summary.txt counts it'], ...
    strcmp(r3.Q.rescue_result{q2}, 'no_basis') && ~r3.Q.rescue_done(q2) && r3.Q.rescue_calls(q2) == 0 && ...
    ~any([E3.q] == q2 & startsWith({E3.kind}, 'rescue')) && ...
    contains(fileread(fullfile(d, 'rc52_one_pass', 'summary.txt')), 'Rescue no_basis: 1 node(s)'), ...
    sprintf('%s, done %d, calls %d', r3.Q.rescue_result{q2}, r3.Q.rescue_done(q2), r3.Q.rescue_calls(q2)));
% under the full ladder the same geometry still rescues k:-2,0 in pass 1 (deep), so the
% the earlier behaviour is unchanged wherever a ladder makes a call
r4 = sa_m3_run(fullfile(d, 'rc52_m4'), [3 1], 'fail', @(dd, b, ex) walls(dd, ex, 'pair'), 'settings', ...
    struct('continuation', struct('box', [3 1])));
rows = sa_check(rows, T, 'older ladder on the same geometry: every ladder makes calls; no node records no_basis', ...
    ~any(strcmp(r4.Q.rescue_result, 'no_basis')) && all(r4.Q.rescue_calls(r4.Q.rescue_done) > 0), ...
    strjoin(r4.Q.rescue_result(r4.Q.rescue_done).', ','));
rmpath(E.stub);
% ---- (3) the full ladder fixtures --------------------------------------------------------------
r27 = sa_t27_rescue(E); r48 = sa_t48_am10(E); r49 = sa_t49_am11_rescue_ring(E);
rows = sa_check(rows, T, 'with the explicit older ladder, T27, T48 and T49 pass', ...
    all([r27.pass]) && all([r48.pass]) && all([r49.pass]), sprintf('T27 %d/%d, T48 %d/%d, T49 %d/%d', ...
    sum([r27.pass]), numel(r27), sum([r48.pass]), numel(r48), sum([r49.pass]), numel(r49)));
end

function c = walls(dd, ex, which)
c = 'S';
switch which
    case 'ringwall'
        % the ring wall of sa_t49_am11_rescue_ring plus the interior wall of sa_t27_rescue: refuse non-exact guesses
        if abs(dd(1) - 3.9) < 1e-9 && abs(dd(2)) < 0.7 && ~ex, c = 'R'; end
        if abs(dd(1) - 1.3) < 1e-9 && abs(dd(2)) < 0.7 && ~ex, c = 'R'; end
    case 'pair'
        if (abs(dd(1) + 0.65) < 1e-9 || abs(dd(1) + 1.3) < 1e-9) && abs(dd(2)) < 1e-9 && ~ex, c = 'R'; end
end
end
