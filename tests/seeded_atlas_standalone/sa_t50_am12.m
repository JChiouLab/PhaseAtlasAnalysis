function rows = sa_t50_am12(E)
%SA_T50_AM12  Milestone-4 T50 (packet v0.3 AM-12): the report-only
%   edge_large_change on planted distances (d > factor x the median of the
%   completed 1-call interior continuation edges; a multi-call edge uses
%   max(d(parent, child), max d_sub / s_sub); fewer than 30 basis edges ->
%   not_computed; the factor is a report leaf) and the second-difference
%   diagnostic on a planted linear branch against a planted jump.
T = 'T50 AM-12 edge_large_change';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't50');
U = fullfile(d, 'U');
seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
r = sa_m3_run(U, [3 3], 'fail', @(dd, b, ex) first_refused(dd, seen));   % two multi-call edges
Q0 = load(fullfile(sa_layout('state', U), 'queue_final.mat'));
E0 = Q0.Q.edges;
basis = find(strcmp({E0.kind}, 'continuation') & strcmp({E0.phase}, 'interior') & strcmp({E0.status}, 'completed') & ...
    [E0.calls] == 1);
% ---- planted distances -------------------------------------------------------------
P = fullfile(d, 'planted'); copyfile(U, P);
Z = load(fullfile(sa_layout('state', P), 'queue_final.mat'));
Ed = Z.Q.edges;
for k = basis, Ed(k).d_parent_child = 0.05; Ed(k).d_sub = 0.05; Ed(k).s_sub = 1; end
for k = find(strcmp({Ed.status}, 'completed') & [Ed.calls] > 1)        % the natural multi-call edges
    Ed(k).d_parent_child = 0.05; Ed(k).d_sub = 0.05 * Ed(k).s_sub;
end
hi = basis([3 5]); lo = basis(7);
Ed(hi(1)).d_parent_child = 0.2; Ed(hi(2)).d_parent_child = 0.16; Ed(lo).d_parent_child = 0.14;
mc = basis(9);                                     % a planted multi-call edge (leaves the basis)
Ed(mc).calls = 3; Ed(mc).d_parent_child = 0.05; Ed(mc).d_sub = [0.01 0.09]; Ed(mc).s_sub = [0.5 0.5];
mc2 = basis(11);
Ed(mc2).calls = 3; Ed(mc2).d_parent_child = 0.05; Ed(mc2).d_sub = [0.02 0.07]; Ed(mc2).s_sub = [0.5 0.5];
Z.Q.edges = Ed;
sa_write('mat', fullfile(sa_layout('state', P), 'queue_final.mat'), Z);
LC = sa_pilot_report(P, 'large_change');
want = sort([hi, mc]);
rows = sa_check(rows, T, ['planted: median 0.05, factor 3: flagged exactly 0.2, 0.16 (not 0.14) and the multi-call edge with ' ...
    'max d_sub/s_sub = 0.18 (not the one with 0.14)'], strcmp(LC.status, 'computed') && abs(LC.median - 0.05) < 1e-15 && ...
    isequal(sort(LC.flagged_ids), sort([Ed(want).id])), sprintf('basis %d; flagged %s', LC.n_basis, mat2str(LC.flagged_ids)));
E2 = load(fullfile(sa_layout('state', P), 'settings_effective.mat'));
E2.settings.report.large_change_factor = 2.5;
sa_write('mat', fullfile(sa_layout('state', P), 'settings_effective.mat'), E2);
LC2 = sa_pilot_report(P, 'large_change');
rows = sa_check(rows, T, 'the factor is the report leaf report.large_change_factor (2.5 also flags 0.14)', ...
    isequal(sort(LC2.flagged_ids), sort([Ed(sort([want lo mc2])).id])), mat2str(LC2.flagged_ids));
r1 = sa_m3_run(fullfile(d, 'small'), [1 1]);
LCs = sa_pilot_report(fullfile(d, 'small'), 'large_change');
rows = sa_check(rows, T, 'fewer than 30 basis edges -> not_computed (box [1 1])', strcmp(LCs.status, 'not_computed') && ...
    isempty(LCs.flagged_ids) && LCs.n_basis < 30, sprintf('basis %d', LCs.n_basis));
% ---- recomputation from attempt files (milestone-3 folders) equals the stored d_sub -------------
Rm = fullfile(d, 'm3form'); copyfile(U, Rm);
Z = load(fullfile(sa_layout('state', Rm), 'queue_final.mat'));
Z.Q.edges = rmfield(Z.Q.edges, {'d_sub', 's_sub'});
sa_write('mat', fullfile(sa_layout('state', Rm), 'queue_final.mat'), Z);
LCa = sa_pilot_report(U, 'large_change'); LCb = sa_pilot_report(Rm, 'large_change');
nmc = sum([LCa.edges.calls] > 1);
rows = sa_check(rows, T, 'without d_sub (milestone-3 folders) the step distances are recomputed from the attempt files: same d per edge', ...
    nmc >= 2 && isequal([LCa.edges.d], [LCb.edges.d]) && isequal([LCa.edges.d_sub_norm], [LCb.edges.d_sub_norm]) && ...
    LCb.dmax_mismatch == 0, sprintf('%d multi-call edges', nmc));
% ---- second difference: linear branch vs jump ----------------------------------------------
force = @(dd, b) jump(dd);
rJ = sa_m3_run(fullfile(d, 'jump'), [3 0], 'force', force);
LCj = sa_pilot_report(fullfile(d, 'jump'), 'large_change');
Ej = LCj.edges;
Q = rJ.Q;
sd = [Ej.second_diff];
qn = arrayfun(@(t) Q.index(Q.edges(t.id).q, 1), Ej);
pn = arrayfun(@(t) Q.index(Q.edges(t.id).p, 1), Ej);
atj = qn == 2 & pn == 1;
smooth = isfinite(sd) & ~atj & ~(pn == 2 & qn == 3);
rows = sa_check(rows, T, 'second difference: ~0 on the planted linear branch, ~1 (> 0.3) at the planted jump', ...
    any(atj & isfinite(sd)) && all(sd(atj & isfinite(sd)) > 0.3) && any(smooth) && max(sd(smooth)) < 1e-10, ...
    sprintf('jump %s; smooth max %.3g', mat2str(sd(atj), 3), max([sd(smooth) 0])));
rmpath(E.stub);
end

function b = jump(dd)
b = '';
if dd(1) > 1.3 - 1e-9, b = 'B'; end
end

function c = first_refused(dd, seen)
% the first call AT the nodes [1.3 0.65] and [-1.3 -0.65] is refused, later ones pass
c = 'S';
if max(abs(abs(dd) - [1.3 0.65])) < 1e-9 && sign(dd(1)) == sign(dd(2))
    k = sprintf('%.6f,%.6f', dd);
    if ~isKey(seen, k), seen(k) = true; c = 'R'; end
end
end
