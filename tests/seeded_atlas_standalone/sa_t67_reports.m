function rows = sa_t67_reports(E)
%SA_T67_REPORTS  The report
%   lines of the mesh trace. Zero solves (stub correctors).
%   - summary.txt: "mesh proposals recorded for n of m calls";
%   - the cut-by-cap line says "count exact" when every cold call's
%     proposals are recorded and "a lower bound" when not;
%   - the pilot report's NMax split uses the recorded P and says for how
%     many calls;
%   - sa_report('cutset') on the run folder: the csv, the histogram and the
%     cap-loss curve, consistent with the run's own cut_by_cap labels.
T = 'T67 reports';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't67');
% run A: every call notes (one NMax-bound cold return above the cap); run B: no call notes
plan = {'cut', 'ok', 'ok', 'ok'};
A = fullfile(d, 'A');
rA = run_planted(A, plan, true);
B = fullfile(d, 'B');
run_planted(B, plan, false);
sa = fileread(fullfile(A, 'summary.txt')); sb = fileread(fullfile(B, 'summary.txt'));
rows = sa_check(rows, T, ['summary.txt: "mesh proposals recorded for 4 of 4 calls" and the cut line "count exact" when ' ...
    'every call notes; "0 of 4" and "a lower bound" when none does'], ...
    contains(sa, 'mesh proposals recorded for 4 of 4 calls') && ...
    contains(sa, 'cold calls cut by the cap (NMax = mesh_reject = 10000): 1 of 4 -- a refusal of the capped search, not evidence that no root exists; count exact') && ...
    contains(sb, 'mesh proposals recorded for 0 of 4 calls') && contains(sb, 'count a lower bound of the calls the cap changed') && ...
    sum([rA.attempts.cut_by_cap]) == 1, sprintf('%s || %s', line_of(sa, 'mesh proposals'), line_of(sb, 'cold calls cut')));
sa_report('pilot', A); sa_report('pilot', B);
pa = fileread(fullfile(A, 'pilot_report.txt')); pb = fileread(fullfile(B, 'pilot_report.txt'));
rows = sa_check(rows, T, ['pilot report: the NMax split uses the recorded P (1 at the cold cap; "from the recorded ' ...
    'proposals ... for 4 of 4 calls"); without traces the M5 line is unchanged (the stub classifier names no NMax-bound ' ...
    'call, so it counts 0 there)'], ...
    contains(pa, 'NMax-bound calls split: at the cold cap (nmax 10000) 1,') && ...
    contains(pa, 'NMax-bound from the recorded proposals (P > the call''s nmax) for 4 of 4 calls') && ...
    contains(pb, 'NMax-bound calls split: at the cold cap (nmax 10000) 0,') && ~contains(pb, 'from the recorded proposals'), ...
    line_of(pa, 'NMax-bound calls split'));
X = sa_report('cutset', A, 10000, {}, struct('out_dir', fullfile(d, 'cut')));
cs = fileread(fullfile(d, 'cut', 'cutset_10000.csv'));
sm = fileread(fullfile(d, 'cut', 'cutset_10000_summary.txt'));
lab = {rA.attempts([rA.attempts.cut_by_cap]).attempt_id};
rows = sa_check(rows, T, ['sa_report(''cutset'', run, 10000): EXACT, the cut set equals the run''s cut_by_cap labels; ' ...
    'cutset_10000.csv, the P histogram and the cap-loss curve written'], X.exact && isequal(X.cut_ids, lab) && ...
    contains(cs, 'A000001,cold,') && contains(sm, 'histogram of P') && contains(sm, 'c = 10000: cut 1'), X.lines{1});
rmpath(E.stub);
end

function res = run_planted(R, plan, note)
s = sa_m3_stub('settings', R, []);
s.routes = struct('continuation', false, 'rescue', false);
s.targets = struct('policy', 'explicit', 'list', [0 0; 1 0; 2 0; 3 0]);
s.cold = struct('mesh_cap', 10000);
s.budget = struct('max_calls', numel(plan));
st = containers.Map('KeyType', 'char', 'ValueType', 'any'); st('n') = 0;
res = seeded_atlas_run(sa_m3_stub('params'), struct(), s, @(th, g, t, o) planted(th, g, t, o, st, plan, note));
end

function [sol, info] = planted(th, g, t, o, st, plan, note)
n = st('n') + 1; st('n') = n;
info = struct('class', 'STUB', 'warning_id', '');
sol = sa_m3_stub('sol', th, t.n, 'A');
if note, sa_mesh_trace('initial', numel(g.x)); end
switch plan{n}
    case 'cut'                     % NMax-bound: proposal 10,801 above the cap, previous mesh returned
        if note, sa_mesh_trace('note', 201); sa_mesh_trace('note', 10801); end
        sol.stats.nmeshpoints = 10801;
        if 10801 > o.nmax, info.warning_id = 'MATLAB:bvp4c:RelTolNotMet'; sol.stats.maxres = 2e-6; end
    otherwise
        if note, sa_mesh_trace('note', 201); end
end
end

function s = line_of(txt, key)
L = strsplit(txt, newline);
k = find(contains(L, key), 1);
if isempty(k), s = ''; else, s = L{k}; end
end
