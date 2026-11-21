function rows = sa_t65_replay(E)
%SA_T65_REPLAY  Replay and
%   kill/resume with mesh traces. Zero solves (stub corrector that notes a
%   trace like the traced bvp4c). The other kill/resume tests (sa_t29_resume,
%   sa_t30_resume_refusals, sa_t38_39_budget_checkpoint, sa_t57_cold_cap) run
%   unchanged in the same suite.
%   - a stub run killed mid-edge (inside a continuation edge, at a substep)
%     resumes equivalent to its uninterrupted twin (sa_m3_equiv, every
%     trace compared without the per-process token);
%   - the interrupted record carries max_proposed_mesh NaN, its retry a
%     fresh trace;
%   - a pure replay of the finished folder (budget = its record count)
%     replays every record with zero corrector calls; the per-record replay
%     time is printed (report-only).
T = 'T65 replay';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't65');
box = [2 1];
U = fullfile(d, 'U');
rU = run_noting(U, box, Inf, false, []);
A = rU.attempts;
k = find(strcmp({A.kind}, 'continuation') & [A.fraction] < 1, 1);
if isempty(k), k = find(strcmp({A.kind}, 'continuation'), 1) + 1; end
K = fullfile(d, 'K');
try, run_noting(K, box, k, false, []); catch, end
rK = run_noting(K, box, Inf, true, []);
D = sa_m3_equiv(U, K);
AK = rK.attempts;
ii = find(strcmp({AK.outcome}, 'interrupted'));
rt = find(strcmp({AK.retry_of}, sprintf('A%06d', k)));
ok = D.pass && D.n_interrupted == 1 && isequal(ii, k) && D.n_trace_compared == numel(A) && numel(rt) == 1 && ...
    isnan(AK(ii).max_proposed_mesh) && isfinite(AK(rt).max_proposed_mesh) && strcmp(AK(rt).kind, 'continuation');
rows = sa_check(rows, T, sprintf(['a stub run killed mid-edge (call %d of %d: %s, fraction %g) resumes equivalent to its ' ...
    'twin with every trace equal; the interrupted record has P NaN and its retry (%s) a fresh trace'], k, numel(A), ...
    A(k).kind, A(k).fraction, sprintf('A%06d', rt)), ok, sprintf('equiv %d, interrupted %s, traces %d of %d; %s', ...
    D.pass, mat2str(ii), D.n_trace_compared, numel(A), strjoin(D.diffs, ' | ')));
% pure replay of the finished folder: zero corrector calls, per-record time printed
n = numel(AK);
cnt = containers.Map('KeyType', 'char', 'ValueType', 'double'); cnt('n') = 0;
t0 = tic;
rP = run_noting(K, box, Inf, true, cnt, n);
tr = toc(t0);
fprintf('T65: pure replay of %d records in %.2f s (%.4f s per record, stub model; report-only)\n', n, tr, tr / n);
rows = sa_check(rows, T, sprintf(['a pure replay of the finished folder (budget = its %d records) re-derives every record ' ...
    'with zero corrector calls (stop completed or budget_calls); per-record time printed (%.4f s per record, ' ...
    'report-only)'], n, tr / n), cnt('n') == 0 && rP.calls == n && any(strcmp(rP.stop_reason, {'budget_calls', 'completed'})), ...
    sprintf('corrector calls %d, stop %s', cnt('n'), rP.stop_reason));
rmpath(E.stub);
end

function res = run_noting(R, box, kill_at, resume, cnt, max_calls)
s = sa_m3_stub('settings', R, box);
if resume, s.resume = true; end
if nargin >= 6, s.budget.max_calls = max_calls; end
st = sa_m3_stub('state', 'kill_at', kill_at);
c0 = sa_m3_stub('corrector', st);
res = seeded_atlas_run(sa_m3_stub('params'), struct(), s, @(th, g, t, o) noting(c0, cnt, th, g, t, o));
end

function [sol, info] = noting(c0, cnt, th, g, t, o)
if ~isempty(cnt), cnt('n') = cnt('n') + 1; end
sa_mesh_trace('initial', numel(g.x));
[sol, info] = c0(th, g, t, o);
if ~isempty(sol), sa_mesh_trace('note', sol.stats.nmeshpoints); end
end
