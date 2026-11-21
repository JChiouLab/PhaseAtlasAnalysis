function rows = sa_t38_39_budget_checkpoint(E)
%SA_T38_39_BUDGET_CHECKPOINT  Budget stops and checkpoints.
%   T38 a budget stop in the middle of an edge leaves its node 'unfinished';
%       resuming with a raised budget completes the run equal to the
%       uninterrupted one, the raised budget recorded in resume_1.mat.
%   T39 checkpoints are compared at their own attempt count (an older
%       queue.mat after a kill passes; a tampered queue.mat or shard is an
%       integrity stop); resume_<n>.mat holds budget, time, pid, attempt
%       count, override; edge_large_change flags exactly the edges of a
%       planted branch jump (stub at R = B = 0.4, lattice step 0.01: one
%       lattice step moves the iterate by ~5e-3 < 1e-2, a branch jump by
%       ~3e-2).
T = 'T38-T39 budget and checkpoints';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't38');
fail = @(dd, b, ex) sa_m3_wall(dd, ex);
box = [2 1];
U = fullfile(d, 'U');
rU = sa_m3_run(U, box, 'fail', fail);
A = rU.attempts;
ks = find(arrayfun(@(a) a.fraction < 1, A), 1);       % first substep (f = 1/2) of a wall edge
K = fullfile(d, 'K');
rK = sa_m3_run(K, box, 'fail', fail, 'settings', struct('budget', struct('max_calls', ks)));
q = find(strcmp(rK.L.key, regexprep(A(ks).target_key, '<-.*$', '')));
ctx = struct('Q', rK.Q, 'shards', rK.shards, 'L', rK.L, 'S', rK.settings);
N = sa_continuation('node_classes', ctx);
cl = N([N.node] == q).class;
rows = sa_check(rows, T, 'T38 budget stop mid-edge: stop budget_calls; the edge''s node is unfinished', ...
    strcmp(rK.stop_reason, 'budget_calls') && rK.calls == ks && strcmp(cl, 'unfinished') && ...
    any(strcmp({rK.Q.edges.status}, 'unfinished')), sprintf('stop %s at %d; node %s %s', rK.stop_reason, ...
    rK.calls, rK.L.key{q}, cl));
rR = sa_m3_run(K, box, 'fail', fail, 'settings', struct('resume', true, 'budget', struct('max_calls', 5000)));
D = sa_m3_equiv(U, K);
Z = load(fullfile(sa_layout('state', K), 'resume_1.mat'));
rows = sa_check(rows, T, 'T38 resume with a raised budget: completed, == uninterrupted; budget in resume_1.mat', ...
    strcmp(rR.stop_reason, 'completed') && D.pass && Z.budget.max_calls == 5000 && Z.budget_stored.max_calls == ks, ...
    strjoin(D.diffs, ' | '));
% ---- T39 checkpoints at their own count ------------------------------------
kc = 23;
C1 = fullfile(d, 'C1');
try, sa_m3_run(C1, box, 'fail', fail, 'kill_at', kc); catch, end
Zq = load(fullfile(sa_layout('state', C1), 'queue.mat'));
C2 = fullfile(d, 'C2'); copyfile(C1, C2);
C3 = fullfile(d, 'C3'); copyfile(C1, C3);
rC = sa_m3_run(C1, box, 'fail', fail, 'settings', struct('resume', true));
D = sa_m3_equiv(U, C1);
rows = sa_check(rows, T, 'T39 queue.mat stamped 20 (every 5 calls) after a kill at call 23: met at its own count, resume == uninterrupted', ...
    Zq.attempt_count == 20 && D.pass && strcmp(rC.stop_reason, 'completed'), strjoin(D.diffs, ' | '));
Zq.Q.slots_used = 99; save(fullfile(sa_layout('state', C2), 'queue.mat'), '-struct', 'Zq');
got = '';
try, sa_m3_run(C2, box, 'fail', fail, 'settings', struct('resume', true)); catch ME, got = [ME.identifier ' ' ME.message]; end
rows = sa_check(rows, T, 'T39 tampered queue.mat -> integrity stop at its own count', ...
    startsWith(got, 'SA:INTEGRITY') && contains(got, 'queue.mat (count 20)'), got);
sf = dir(fullfile(sa_layout('state', C3), 'roots', '*.mat')); Zs = load(fullfile(sf(1).folder, sf(1).name));
Zs.store.roots(1).representative_profile(1) = Zs.store.roots(1).representative_profile(1) + 1;
save(fullfile(sf(1).folder, sf(1).name), '-struct', 'Zs');
got = '';
try, sa_m3_run(C3, box, 'fail', fail, 'settings', struct('resume', true)); catch ME, got = [ME.identifier ' ' ME.message]; end
rows = sa_check(rows, T, 'T39 tampered shard -> integrity stop at the shard''s own count', ...
    startsWith(got, 'SA:INTEGRITY') && contains(got, 'shard'), got);
Z = load(fullfile(sa_layout('state', C1), 'resume_1.mat'));
rows = sa_check(rows, T, 'T39 resume_1.mat: budget, checkpoint, time, pid, attempt count, override, interrupted ids', ...
    all(isfield(Z, {'budget', 'checkpoint', 'time', 'pid', 'attempt_count', 'override', 'interrupted'})) && ...
    Z.attempt_count == kc && isequal(Z.interrupted, kc), sprintf('attempt_count %d', Z.attempt_count));
% ---- T39 edge_large_change on a planted jump ---------------------------------
o = struct('lattice', struct('step', 0.01), 'routes', struct('rescue', false), ...
    'params_override', struct('R', 0.4, 'B', 0.4, 'Lh', 1.5));
force = @(dd, b) jump(dd);
rL = sa_m3_run(fullfile(d, 'L'), [2 2], 'force', force, 'settings', o);
Ed = rL.Q.edges;
j = find(all(rL.L.index == [1 0], 2));
into = [Ed.q] == j & strcmp({Ed.status}, 'completed');
rows = sa_check(rows, T, 'T39 edge_large_change set exactly on the edges landing on the planted branch jump', ...
    any(into) && isequal([Ed.large_change], into) && all([Ed(into).d_parent_child] > 1e-2) && ...
    all([Ed(~into).d_parent_child] < 1e-2), sprintf('%d edges, %d flagged; d jump %s; others max %.3g', ...
    numel(Ed), sum([Ed.large_change]), mat2str([Ed(into).d_parent_child], 3), max([Ed(~into).d_parent_child])));
rmpath(E.stub);
end

function b = jump(dd)
b = '';
if abs(dd(1) - 0.01) < 1e-12 && dd(2) == 0, b = 'B'; end
end
