function rows = sa_t60_reports(E)
%SA_T60_REPORTS  The report lines for the cold-only mesh cap, the secant-only
%   rescue and the per-root pass. Zero solves (stub corrector). Run A (smoke
%   mode, box [3 1]): pass on (a branch-B cold root at the anchor starts no
%   lineage), cold cap 10,000 with one planted NMax-bound cold return (cut),
%   secant-only rescue with walls at k:-1,0 and k:-2,0 and one search pass
%   (k:-2,0 stays no_basis). Run B: secant-only rescue, nothing refused (the
%   b_bound NOTE). summary.txt: the search line prints rescue.rungs,
%   rescue.ring and root_pass.enabled; the cut-by-cap line; the no_basis
%   count; the pass rounds. Pilot report (rebuilt, a pure function of the
%   folder): cut_by_cap reason class, the NMax split, the pass yield and
%   NOTE, pilot/root_pass.csv = Q.root_pass, the no_basis count, b_bound
%   conservative under secant-only; the thinning note (MULTIPLICITY
%   PARTIALLY TESTED, identity.csv multiplicity_tested 0).
T = 'T60 reports';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't60');
% ---- run A ----------------------------------------------------------------------------------
A = fullfile(d, 'A');
s = sa_m3_stub('settings', A, [3 1]);
s.families = {'X_X', 'Y_X'};
s.targets = struct('policy', 'explicit', 'list', [0 0; 1 0]);
s.continuation.starts = 'anchor';
s.rescue = struct('rungs', {{'secant'}}, 'ring', false);
s.search = struct('max_passes', 1);
s.root_pass = struct('enabled', true);
s.cold = struct('mesh_cap', 10000);
st = sa_m3_stub('state', 'fail', @(dd, b, ex) walls(dd, ex));
c0 = sa_m3_stub('corrector', st);
cnt = containers.Map('KeyType', 'char', 'ValueType', 'double'); cnt('n') = 0;
rA = seeded_atlas_run(sa_m3_stub('params'), struct(), s, @(th, g, t, o) wrap(th, g, t, o, c0, cnt));
sm = fileread(fullfile(A, 'summary.txt'));
nb = sum(strcmp(rA.Q.rescue_result, 'no_basis') & ~rA.Q.rescue_done);
rows = sa_check(rows, T, ['summary.txt: search line prints rescue.rungs {secant}, rescue.ring 0, root_pass.enabled 1; ' ...
    '"cold calls cut by the cap (NMax = mesh_reject = 10000): 1 of 4"; the no_basis count; the pass rounds'], ...
    contains(sm, 'route settings: rescue.rungs {secant}; rescue.ring 0; root_pass.enabled 1') && ...
    contains(sm, 'cold calls cut by the cap (NMax = mesh_reject = 10000): 1 of 4 -- a refusal of the capped search, not evidence that no root exists') && ...
    nb == 1 && contains(sm, 'Rescue no_basis: 1 node(s)') && contains(sm, 'Per-root pass (root_pass.enabled 1') && ...
    contains(sm, sprintf('pass 1 interior round %d:', numel(rA.Q.root_pass))) && sum([rA.attempts.cut_by_cap]) == 1, ...
    sprintf('rounds %d, cut %d, no_basis %d', numel(rA.Q.root_pass), sum([rA.attempts.cut_by_cap]), nb));
P = sa_report('pilot', A);
pr = fileread(fullfile(A, 'pilot_report.txt'));
rp = strsplit(strtrim(fileread(fullfile(sa_layout('state', A), 'pilot', 'root_pass.csv'))), newline);
vb = fileread(fullfile(sa_layout('state', A), 'pilot', 'validity_by_check.csv'));
rows = sa_check(rows, T, ['pilot report: cut_by_cap reason class (cold 1), the NMax split line, the pass yield (targets, ' ...
    'calls, new roots, calls per new root) and its NOTE, root_pass.csv one row per round, the no_basis count'], ...
    contains(pr, 'reason class cut_by_cap') && contains(vb, 'cut_by_cap,cold,0,1,1') && ...
    contains(pr, 'NMax-bound calls split: at the cold cap (nmax 10000)') && ...
    contains(pr, sprintf('per-root pass (root_pass.enabled 1): targets posed %d', P.yield.root_pass.targets)) && ...
    P.yield.root_pass.targets > 0 && P.yield.root_pass.new_roots > 0 && contains(pr, 'NOTE: the per-root pass is on') && ...
    numel(rp) - 1 == numel(rA.Q.root_pass) && P.yield.no_basis == 1 && ...
    contains(pr, 'nodes with a zero-call ladder (no_basis: no lineage pair at q-e, q-2e): 1'), ...
    sprintf('targets %d, new %d, calls/new %.3g, csv rows %d', P.yield.root_pass.targets, P.yield.root_pass.new_roots, ...
    P.yield.root_pass.calls_per_new_root, numel(rp) - 1));
% ---- run B: b_bound note ------------------------------------------------------------------------
B = fullfile(d, 'B');
rB = sa_m3_run(B, [1 1], 'settings', struct('rescue', struct('rungs', {{'secant'}}, 'ring', false)));
PB = sa_report('pilot', B);
pb = fileread(fullfile(B, 'pilot_report.txt'));
rows = sa_check(rows, T, ['b_bound under secant-only rescue with no refused node: charged as before and stated ' ...
    'conservative (secant-only per-node bound 2 x n_pools = 4)'], PB.projection.n_refused_measured == 0 && ...
    contains(pb, 'conservative under secant-only rescue (its per-node bound is 2 x n_pools = 4 calls, not max_calls_per_node = 16)') && ...
    PB.projection.b_bound == 2 * 4 + 1 + 16, sprintf('calls %d; b = %g', rB.calls, PB.projection.b));
% ---- thinning note ---------------------------------------------------------------------------------
Bt = fullfile(d, 'B_thin'); copyfile(B, Bt);
Z = load(fullfile(sa_layout('state', Bt), 'settings_effective.mat'));
Z.settings.targets.policy = 'all_seedable'; Z.settings.targets.list = zeros(0, 2);
Z.settings.targets.thin = struct('families', {{'X_X'}}, 'stride', 2, 'offset', [0 1]);
sa_write('mat', fullfile(sa_layout('state', Bt), 'settings_effective.mat'), Z);
Pt = sa_report('pilot', Bt);
pt = fileread(fullfile(Bt, 'pilot_report.txt'));
idc = fileread(fullfile(sa_layout('state', Bt), 'pilot', 'identity.csv'));
nt = sa_report('multiplicity_note', Z.settings);
keep = sa_targets('thin_keep', Z.settings, 'X_X', rB.L.index);
cen = load(fullfile(sa_layout('state', Bt), 'settings_effective.mat'), 'census'); cen = cen.census;
fi = find(ismember(cen.names, Z.targets.families));
want = 0;
for f = fi(:).'
    if strcmp(cen.names{f}, 'X_X'), want = want + sum(cen.seedable(:, f) & keep(:));
    else, want = want + sum(cen.seedable(:, f)); end
end
rows = sa_check(rows, T, ['thinning note: MULTIPLICITY PARTIALLY TESTED (families, stride, offset) in section 7, ' ...
    'identity.csv multiplicity_tested 0; the projection''s cold line uses the thinned count'], ...
    startsWith(nt, 'MULTIPLICITY PARTIALLY TESTED: cold targets thinned for X_X (stride 2, offset [0 1])') && ...
    contains(pt, 'MULTIPLICITY PARTIALLY TESTED') && contains(idc, 'multiplicity_tested,0') && ...
    ~Pt.identity.multiplicity_tested && Pt.projection.cold_all == want && want < sum(sum(cen.seedable(:, fi))), ...
    sprintf('cold %d (unthinned %d)', Pt.projection.cold_all, sum(sum(cen.seedable(:, fi)))));
rmpath(E.stub);
end

function [sol, info] = wrap(th, g, t, o, c0, cnt)
% call 2 (cold X_X at k:1,0) returns an NMax-bound iterate under the cap
n = cnt('n') + 1; cnt('n') = n;
[sol, info] = c0(th, g, t, o);
if n == 2 && ~isempty(sol)
    sol.stats.nmeshpoints = 12000; sol.stats.maxres = 2e-6;
    info.warning_id = 'MATLAB:bvp4c:RelTolNotMet';
end
end

function c = walls(dd, ex)
c = 'S';
if (abs(dd(1) + 0.65) < 1e-9 || abs(dd(1) + 1.3) < 1e-9) && abs(dd(2)) < 1e-9 && ~ex, c = 'R'; end
end
