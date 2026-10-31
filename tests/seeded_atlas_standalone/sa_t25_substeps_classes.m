function rows = sa_t25_substeps_classes(E)
%SA_T25_SUBSTEPS_CLASSES  Milestone-3 T25: substeps (non-lattice
%   coordinates) are route records only: never stored, never parents of other
%   edges, never coverage; node classes partition the box nodes.
T = 'T25 substeps and node classes';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't25');
res = sa_m3_run(fullfile(d, 'a'), [3 1], 'fail', @(dd, b, ex) sa_m3_wall(dd, ex));
A = res.attempts;
sub = A(arrayfun(@(a) a.fraction < 1, A));
rows = sa_check(rows, T, 'substep records exist and carry substep outcomes only', ...
    ~isempty(sub) && all(startsWith({sub.outcome}, 'substep_')), sprintf('%d substeps', numel(sub)));
rows = sa_check(rows, T, 'no substep record holds a root id', all(cellfun(@isempty, {sub.root_id})));
disc = {};
for k = res.shards.keys
    st = res.shards(k{1});
    for r = 1:numel(st.roots), disc = [disc, {st.roots(r).discoveries.attempt_id}]; end %#ok<AGROW>
end
isnode = arrayfun(@(id) any(strcmp({A(strcmp({A.outcome}, 'accepted')).attempt_id}, id)), disc);
fr = arrayfun(@(id) A(strcmp({A.attempt_id}, id)).fraction, disc);
rows = sa_check(rows, T, 'every stored discovery is an accepted node return (fraction 1)', ...
    all(isnode) && all(fr == 1) && numel(disc) == sum(strcmp({A.outcome}, 'accepted')), ...
    sprintf('%d discoveries', numel(disc)));
sk = {sub.attempt_id};
par = {};
for e = res.Q.edges, par{end+1} = e.parent_attempt; end %#ok<AGROW>
for l = res.Q.lin, par = [par, l.node_attempt(~cellfun(@isempty, l.node_attempt))]; end %#ok<AGROW>
rows = sa_check(rows, T, 'no substep is the parent iterate of an edge or a lineage node', ...
    ~any(ismember(par, sk)), sprintf('%d parent references', numel(par)));
ctx = struct('Q', res.Q, 'shards', res.shards, 'L', res.L, 'S', res.settings);
N = sa_continuation('node_classes', ctx);
cls = {'reached', 'refused', 'unattempted', 'policy_skipped', 'unfinished'};
one = all(cellfun(@(c) any(strcmp(c, cls)), {N.class}));
rows = sa_check(rows, T, 'node classes partition the box nodes (one class each, counts sum to the box)', ...
    one && numel(N) == sum(res.Q.box) && numel(unique([N.node])) == numel(N) && ...
    sum(cellfun(@(c) sum(strcmp({N.class}, c)), cls)) == numel(N), ...
    strjoin(cellfun(@(c) sprintf('%s %d', c, sum(strcmp({N.class}, c))), cls, 'UniformOutput', false), ', '));
reached_by_sub = arrayfun(@(n) strcmp(n.class, 'reached') && n.n_roots == 0, N);
rows = sa_check(rows, T, 'a node is reached only by a stored root (substeps never coverage)', ~any(reached_by_sub));
rmpath(E.stub);
end
