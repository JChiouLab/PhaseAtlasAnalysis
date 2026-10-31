function rows = sa_t28_lineages(E)
%SA_T28_LINEAGES  Milestone-3 T28 (RC-3, M-C): merge by membership;
%   call-free merges use no slot; cold-start preference; M-C propagation with
%   planted roots of two branches (sa_m3_mc_run).
T = 'T28 lineages and M-C';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't28');
res = sa_m3_mc_run(fullfile(d, 'a'));
Q = res.Q; L = res.L; A = res.attempts;
cold = A(strcmp({A.kind}, 'cold') & strcmp({A.outcome}, 'accepted'));
rows = sa_check(rows, T, 'six cold roots planted: X_X (branch A) and Y_X (branch B) at three nodes', ...
    numel(cold) == 6, sprintf('%s', strjoin(strcat({cold.family}, '@', {cold.target_key}), ' ')));
rows = sa_check(rows, T, 'two lineages: the anchor (X_X at k:0,0) and ONE non-anchor start; one slot used', ...
    numel(Q.lin) == 2 && strcmp(Q.lin(1).kind, 'anchor') && strcmp(Q.lin(2).kind, 'start') && Q.slots_used == 1, ...
    sprintf('%d lineages, slots %d', numel(Q.lin), Q.slots_used));
yx = cold(strcmp({cold.family}, 'Y_X'));
rows = sa_check(rows, T, 'cold-start preference: the non-anchor lineage starts at the B root nearest the anchor (k:0,0)', ...
    strcmp(Q.lin(2).start_root, yx(strcmp({yx.target_key}, 'k:0,0')).root_id), Q.lin(2).start_root);
own = @(r) sa_continuation('owner', Q, r);
xx = cold(strcmp({cold.family}, 'X_X'));
oA = arrayfun(@(a) own(a.root_id), xx); oB = arrayfun(@(a) own(a.root_id), yx);
rows = sa_check(rows, T, ['call-free merges: the far A cold roots joined lineage 1, the far B cold roots ' ...
    'joined lineage 2, consuming no slot'], all(oA == 1) && all(oB == 2), sprintf('A owners %s, B owners %s', ...
    mat2str(oA), mat2str(oB)));
m = find(all(L.index == [1 1], 2));
rows = sa_check(rows, T, 'merge by membership: lineage 2 lands on lineage 1''s root at k:1,1 and records the merge', ...
    isequal(Q.lin(2).merged_with, 1) && strcmp(Q.lin(2).node_root{m}, Q.lin(1).node_root{m}) && ...
    ~Q.lin(2).expandable(m), sprintf('merged_with %s', mat2str(Q.lin(2).merged_with)));
box = find(Q.box).';
nr = arrayfun(@(q) sa_continuation('n_roots', struct('shards', res.shards, 'L', L), q), box);
rows = sa_check(rows, T, 'M-C: every box node except k:1,1 holds two distinct roots (A and B branches)', ...
    sum(nr == 2) == numel(box) - 1 && nr(box == m) == 1, sprintf('roots per node %s', mat2str(nr)));
res0 = sa_m3_mc_run(fullfile(d, 'noslot'), struct('continuation', struct('max_lineages', 0)));
nr0 = arrayfun(@(q) sa_continuation('n_roots', struct('shards', res0.shards, 'L', L), q), box);
rows = sa_check(rows, T, 'max_lineages = 0: no non-anchor lineage; the 3 B roots recorded without a slot', ...
    numel(res0.Q.lin) == 1 && numel(res0.Q.no_slot) == 3 && sum(nr0 == 2) == 3, ...
    sprintf('no_slot %s', strjoin(res0.Q.no_slot, ',')));
rmpath(E.stub);
end
