function rows = sa_t45_targets(E)
%SA_T45_TARGETS  Milestone-4 T45 (packet v0.3 section 2b): central_nine centred
%   on the node nearest d = 0 (equal to the M2/M3 block when the anchor is 0;
%   ties by index distance to the anchor, then lattice order, on an offset
%   lattice); missing members listed, never backfilled; non-seedable members
%   policy-skipped; pilot_sample (7 x 7 over the native range, snapped,
%   de-duplicated, centre included); the cross region (lines through the
%   central node, ring included).
T = 'T45 targets';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
% ---- planted lattices (sa_lattice field layout) --------------------------------
L0 = lat([0 0], 0.15, [-13 -13], [13 13]);
c = sa_modes('central_node', L0);
[n9, miss] = sa_modes('central_nine', L0);
old = find(all(abs(L0.index) <= 1, 2));
rows = sa_check(rows, T, 'anchor at d = 0: central node = the anchor; central_nine = the M2/M3 block around index 0', ...
    isequal(L0.index(c, :), [0 0]) && isequal(n9, old) && isempty(miss));
L1 = lat([0.3 0.1], 0.2, [-8 -8], [8 8]);           % ties: axis 1 {-1,-2}, axis 2 {0,-1}
c1 = sa_modes('central_node', L1);
d1 = sqrt(sum(L1.coordinate.^2, 2));
cand = find(abs(d1 - min(d1)) < 1e-12);
rows = sa_check(rows, T, 'offset anchor: nearest d = 0 with distance ties broken by the smaller index distance to the anchor', ...
    numel(cand) == 4 && isequal(L1.index(c1, :), [-1 0]), sprintf('%d tied nodes; centre %s', numel(cand), ...
    mat2str(L1.index(c1, :))));
L2 = L0; L2.coordinate(:) = 1; L2.coordinate(find(all(L2.index == [0 1], 2)), :) = [0.1 0];   %#ok<FNDSB>
L2.coordinate(find(all(L2.index == [1 0], 2)), :) = [0 0.1]; %#ok<FNDSB>
c2 = sa_modes('central_node', L2);
rows = sa_check(rows, T, 'equal distance and equal index distance: lattice order decides (planted)', ...
    isequal(L2.index(c2, :), [0 1]), mat2str(L2.index(c2, :)));
L3 = lat([0.3 0], 0.15, [-2 -3], [3 3]);             % d = 0 lies outside: centre on the edge k1 = -2
[n93, miss3] = sa_modes('central_nine', L3);
kc = L3.index(sa_modes('central_node', L3), :);
rows = sa_check(rows, T, 'centre on the lattice edge: 6 members, the 3 off-lattice ones listed as missing', ...
    isequal(kc, [-2 0]) && numel(n93) == 6 && size(miss3, 1) == 3 && all(miss3(:, 1) == -3), mat2str(miss3));
% ---- sa_targets with a planted census ---------------------------------------------
C = census(L1, abs(L1.coordinate(:, 1)) <= 1.21 & abs(L1.coordinate(:, 2)) <= 1.21);
C.seedable(find(all(L1.index == [-1 1], 2)), 1) = false; %#ok<FNDSB>
S = struct('families', {{}}, 'targets', struct('policy', 'central_nine', 'list', [], 'per_family', struct()), ...
    'pilot', struct('sample_per_axis', 7));
Tt = sa_targets(L1, C, S);
ki = vertcat(Tt.list.index);
rows = sa_check(rows, T, 'sa_targets central_nine on the offset lattice: 8 seedable members around [-1 0]; the non-seedable one policy-skipped', ...
    numel(Tt.list) == 8 && all(max(abs(ki - [-1 0]), [], 2) <= 1) && numel(Tt.skipped) == 1 && ...
    isequal(Tt.skipped(1).index, [-1 1]), sprintf('%d targets, %d skipped', numel(Tt.list), numel(Tt.skipped)));
% ---- pilot_sample -----------------------------------------------------------------------
Cn = census(L0, all(abs(L0.index) <= 12, 2));            % noxa-like: native +-1.8 at step 0.15
Sp = struct('pilot', struct('sample_per_axis', 7));
ps = sa_modes('pilot_sample', L0, Cn, Sp);
kk = unique(L0.index(ps, 1)).';
rows = sa_check(rows, T, 'pilot_sample noxa-like: 7 x 7 = 49 nodes at k = -12:4:12 (exact snap), centre included', ...
    numel(ps) == 49 && isequal(kk, -12:4:12) && any(all(L0.index(ps, :) == 0, 2)), mat2str(kk));
Cs = census(L0, all(abs(L0.index) <= 2, 2));             % narrow range: snapping duplicates
ps2 = sa_modes('pilot_sample', L0, Cs, Sp);
rows = sa_check(rows, T, 'pilot_sample on a 5-node native range: de-duplicated (5 x 5 = 25), centre included', ...
    numel(ps2) == 25 && numel(unique(ps2)) == 25 && any(all(L0.index(ps2, :) == 0, 2)), sprintf('%d', numel(ps2)));
S.targets.policy = 'pilot_sample';
Cn1 = census(L1, all(abs(L1.index) <= 6, 2));
Tp = sa_targets(L1, Cn1, S);
ns = sa_modes('pilot_sample', L1, Cn1, Sp); [n9b, ~] = sa_modes('central_nine', L1);
rows = sa_check(rows, T, 'targets.policy pilot_sample = central_nine U pilot_sample (lattice order, seedable only)', ...
    isequal(sort([Tp.list.node]).', unique([ns; n9b])) && issorted([Tp.list.node]), ...
    sprintf('%d targets', numel(Tp.list)));
% ---- cross --------------------------------------------------------------------------------
x = sa_modes('cross', L0);
ring = ~Cn.in_range;
rows = sa_check(rows, T, 'cross: the two lattice lines through the central node, 27 + 27 - 1 = 53 nodes, 4 of them in the ring', ...
    sum(x) == 53 && all(L0.index(x, 1) == 0 | L0.index(x, 2) == 0) && sum(x & ring) == 4, ...
    sprintf('%d nodes, ring %d', sum(x), sum(x & ring)));
x1 = sa_modes('cross', L1);
rows = sa_check(rows, T, 'cross on the offset lattice runs through the central node [-1 0]', ...
    sum(x1) == 17 + 17 - 1 && all(L1.index(x1, 1) == -1 | L1.index(x1, 2) == 0));
end

function L = lat(anchor, step, kmin, kmax)
% the sa_lattice layout: rectangle, ascending lexicographic (axis 1 slowest)
[a, b] = ndgrid(kmin(1):kmax(1), kmin(2):kmax(2));
a = a.'; b = b.';
idx = [a(:) b(:)];
L = struct('step', step, 'anchor', anchor, 'margin', 1e-4, 'n_coordinates', 2, 'index_min', kmin, ...
    'index_max', kmax, 'index', idx, 'coordinate', anchor + step * idx, 'half_targets', zeros(size(idx, 1), 4), ...
    'n_nodes', size(idx, 1));
L.key = arrayfun(@(r) sprintf('k:%d,%d', idx(r, 1), idx(r, 2)), (1:size(idx, 1)).', 'UniformOutput', false);
L.rule = 'planted';
end

function C = census(L, inr)
n = L.n_nodes;
C = struct('names', {{'X_X'}}, 'counts', sum(inr), 'seedable', inr(:), 'in_range', inr(:), ...
    'cell', repmat([1 1], n, 1), 'mirror', repmat([1 1], n, 1), 'outside_reason', {repmat({'planted'}, n, 1)}, ...
    'family_types', {{'X', 'X'}});
end
