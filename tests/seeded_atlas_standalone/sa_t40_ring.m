function rows = sa_t40_ring(E)
%SA_T40_RING  Milestone-3 T40 (packet 2h, K-20): with search.ring_last the
%   ring (lattice nodes outside the library's native range; |k| = 6 on the
%   stub's 13 x 13 lattice) runs after the interior: interior outcomes are
%   identical with budget.ring_calls 0, a cap, or Inf; ring-phase calls
%   target ring nodes only and no interior work follows them; reaching the
%   cap stops at budget_ring, and a resume with a raised cap equals the
%   uninterrupted run. With ring_last = false the ring joins one phase.
T = 'T40 ring';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't40');
rI = sa_m3_run(fullfile(d, 'inf'), []);
r0 = sa_m3_run(fullfile(d, 'zero'), [], 'settings', struct('budget', struct('ring_calls', 0)));
r10 = sa_m3_run(fullfile(d, 'ten'), [], 'settings', struct('budget', struct('ring_calls', 10)));
A = rI.attempts;
ph = {A.phase};
first_ring = find(strcmp(ph, 'ring'), 1);
nint = sum(~strcmp(ph, 'ring'));
D0 = sa_m3_equiv(fullfile(d, 'inf'), fullfile(d, 'zero'), struct('prefix_only', true));
D10 = sa_m3_equiv(fullfile(d, 'inf'), fullfile(d, 'ten'), struct('prefix_only', true));
rows = sa_check(rows, T, 'ring_calls 0 / 10 / Inf: identical interior outcomes (0 and 10 are prefixes of Inf)', ...
    D0.pass && D10.pass && r0.calls == nint && r10.calls == nint + 10 && ~isempty(first_ring) && ...
    first_ring == nint + 1, sprintf('interior calls %d; runs %d/%d/%d', nint, r0.calls, r10.calls, rI.calls));
rows = sa_check(rows, T, 'cap reached -> stop budget_ring (cap 0 and cap 10); uncapped completes', ...
    strcmp(r0.stop_reason, 'budget_ring') && strcmp(r10.stop_reason, 'budget_ring') && ...
    strcmp(rI.stop_reason, 'completed'), sprintf('%s / %s / %s', r0.stop_reason, r10.stop_reason, rI.stop_reason));
ring = A(strcmp(ph, 'ring'));
kk = cellfun(@(k) max(abs(sscanf(regexprep(k, '<-.*$', ''), 'k:%d,%d'))), {ring.target_key});
rows = sa_check(rows, T, 'ring-phase calls target ring nodes only; no interior call after the first ring call', ...
    all(kk == 6) && ~any(~strcmp(ph(first_ring:end), 'ring')), sprintf('%d ring calls', numel(ring)));
ctx = struct('Q', rI.Q, 'shards', rI.shards, 'L', rI.L, 'S', rI.settings);
N = sa_continuation('node_classes', ctx);
rows = sa_check(rows, T, 'all 48 ring nodes reached in the uncapped run (ours-only report)', ...
    sum([N.ring]) == 48 && all(strcmp({N([N.ring]).class}, 'reached')), '');
rR = sa_m3_run(fullfile(d, 'ten'), [], 'settings', struct('resume', true, 'budget', struct('ring_calls', Inf)));
D = sa_m3_equiv(fullfile(d, 'inf'), fullfile(d, 'ten'));
rows = sa_check(rows, T, 'budget_ring stop + resume with a raised cap == uninterrupted', ...
    D.pass && strcmp(rR.stop_reason, 'completed'), strjoin(D.diffs, ' | '));
rF = sa_m3_run(fullfile(d, 'mixed'), [], 'settings', struct('search', struct('ring_last', false)));
rows = sa_check(rows, T, 'ring_last = false: one phase ''all''; ring nodes reached in it', ...
    isequal(rF.Q.phases_run, {'all'}) && ~any(strcmp({rF.attempts.phase}, 'ring')) && rF.n_roots == 169, ...
    sprintf('%d roots', rF.n_roots));
rmpath(E.stub);
end
