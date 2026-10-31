function T = sa_targets(L, C, S, varargin)
%SA_TARGETS  The declared target set: (family, lattice node) pairs.
%
%   T = SA_TARGETS(L, C, S) from the lattice L, the seedability census C and
%   the resolved settings S. Families: S.families, or every family with a
%   nonzero seedable count. Unknown family names are rejected. Policies:
%     central_nine    the node nearest the
%                     symmetric point d = 0 (sa_modes 'central_node': ties by
%                     index distance to the anchor, then lattice order) and
%                     its lattice neighbours (|k - k_c| <= 1 on every axis; the
%                     3x3 block when there are two axes). With the anchor at
%                     d = 0 this is the block around index 0. Members
%                     that are not lattice nodes are listed as skipped ('not a
%                     lattice node'), never backfilled.
%     pilot_sample    central_nine plus sa_modes('pilot_sample') (about
%                     pilot.sample_per_axis^2 nodes evenly spaced over the
%                     native range); lattice order
%     explicit        the rows of S.targets.list
%     all_seedable    every seedable node
%     per_seed        one target per library seed (production default). A
%                     seed is (family, library cell): the seedable nodes of
%                     the family that map to that cell (the mirror cell
%                     follows from the cell). Its target is the anchor node
%                     when the anchor is one of them (so the anchor node gets
%                     the same cold calls as under all_seedable and the anchor
%                     lineage starts as before); otherwise the node nearest
%                     the cell's own coordinate (Euclidean, in native
%                     spacings per axis, C.offset; distances within 1e-9
%                     spacings are ties, broken by lattice order). The other
%                     nodes of the seed go to T.skipped, reason 'per_seed:
%                     seed cell [r c] targeted at node <key>' (policy-skipped,
%                     never refused). Only the cold calls change: later
%                     phases visit their own node set (sa_continuation
%                     'targets'), not T.list
%     per_seed_centre per_seed's targets with the cold calls ordered
%                     nearest d = 0 first (ties: per_seed order); a run
%                     option (milestone 8), never a default
%     nearest_anchor  per family, the S.targets.per_family.(family) seedable
%                     nodes nearest the anchor (Euclidean index distance;
%                     ties by lattice order); families absent from per_family
%                     contribute nothing
%   Within each family only SEEDABLE nodes become targets; a requested node
%   that is not seedable is listed in T.skipped with its reason (it is
%   policy-skipped, not refused). Order: family order, then lattice order
%   (nearest_anchor: distance order).
%
%   T.list: struct array with family, left_type, right_type, node (lattice
%   row), index, coordinate, key, cell_left, cell_right.
%
%   Thinning lever (never a default): with
%   targets.thin.families non-empty (policy all_seedable only; refused with
%   per_seed), a pair
%   (f, node k) with f thinned is kept iff (k_c - offset_c) mod stride = 0 on
%   every axis (k = lattice index from the anchor; offset [] = zeros); the
%   anchor node is always kept, so thinning never changes the anchor's own
%   targets. Dropped pairs go to T.skipped with reason 'thinned: stride s
%   offset o' (policy-skipped, never refused). Thinned family names must be
%   families of the model's census.
%   keep = SA_TARGETS('thin_keep', S, family, index)   the rule on index rows
%   n = SA_TARGETS('seed_counts', C)   library seeds per census family (the
%       per_seed target counts; 1 x numel(C.names))

if ischar(L) && strcmp(L, 'thin_keep')
    T = thin_keep(C, S, varargin{1});
    return
end
if ischar(L) && strcmp(L, 'seed_counts')
    T = zeros(1, numel(C.names));
    for f = 1:numel(C.names)
        T(f) = size(unique(C.cell(C.seedable(:, f), :), 'rows'), 1);
    end
    return
end
fam = S.families;
if isempty(fam)
    fam = C.names(C.counts > 0);
end
for k = 1:numel(fam)
    if ~any(strcmp(fam{k}, C.names))
        error('SA:TARGETS:E_FAMILY', 'Family %s is not a family of this model (%s).', ...
            fam{k}, strjoin(C.names, ', '));
    end
end
thin = struct('families', {{}}, 'stride', 1, 'offset', []);
if isfield(S.targets, 'thin'), thin = S.targets.thin; end
for k = 1:numel(thin.families)
    if ~any(strcmp(thin.families{k}, C.names))
        error('SA:TARGETS:E_FAMILY', 'targets.thin.families names %s, which is not a family of this model (%s).', ...
            thin.families{k}, strjoin(C.names, ', '));
    end
end
if any(strcmp(S.targets.policy, {'per_seed', 'per_seed_centre'})) && ~isempty(thin.families)
    error('SA:TARGETS:E_THIN', ['targets.thin applies to targets.policy all_seedable only ' ...
        '(policy per_seed has one cold call per library seed already).']);
end
if strcmp(S.targets.policy, 'nearest_anchor')
    extra = setdiff(fieldnames(S.targets.per_family), fam);
    if ~isempty(extra)
        error('SA:TARGETS:E_FAMILY', ['targets.per_family names %s, which is not ' ...
            'among the selected families.'], strjoin(extra', ', '));
    end
end

T = struct('families', {fam}, 'list', [], 'skipped', []);
if any(strcmp(S.targets.policy, {'central_nine', 'pilot_sample'}))
    [c9, miss9] = sa_modes('central_nine', L);
    if strcmp(S.targets.policy, 'pilot_sample')
        c9 = unique([c9; sa_modes('pilot_sample', L, C, S)]);
    end
end
list = struct('family', {}, 'left_type', {}, 'right_type', {}, 'node', {}, ...
    'index', {}, 'coordinate', {}, 'key', {}, 'cell_left', {}, 'cell_right', {});
skipped = struct('family', {}, 'index', {}, 'reason', {});
for k = 1:numel(fam)
    f = find(strcmp(fam{k}, C.names));
    switch S.targets.policy
        case {'central_nine', 'pilot_sample'}
            nodes = c9;
            for r = 1:size(miss9, 1)
                skipped(end+1) = struct('family', fam{k}, 'index', miss9(r, :), ...
                    'reason', 'not a lattice node'); %#ok<AGROW>
            end
        case 'explicit'
            nodes = zeros(0, 1);
            for r = 1:size(S.targets.list, 1)
                hit = find(all(L.index == S.targets.list(r, :), 2));
                if isempty(hit)
                    skipped(end+1) = struct('family', fam{k}, ...
                        'index', S.targets.list(r, :), 'reason', 'not a lattice node'); %#ok<AGROW>
                else
                    nodes(end+1, 1) = hit; %#ok<AGROW>
                end
            end
        case 'all_seedable'
            nodes = find(C.seedable(:, f));
        case {'per_seed', 'per_seed_centre'}
            [nodes, sk] = per_seed(L, C, f, fam{k});
            skipped = [skipped, sk]; %#ok<AGROW>
        case 'nearest_anchor'
            nodes = zeros(0, 1);
            if isfield(S.targets.per_family, fam{k})
                cand = find(C.seedable(:, f));
                dist = sqrt(sum(L.index(cand, :).^2, 2));
                [~, ord] = sortrows([dist cand]);
                nodes = cand(ord(1:min(numel(cand), S.targets.per_family.(fam{k}))));
            end
    end
    if any(strcmp(fam{k}, thin.families))
        keep = thin_keep(S, fam{k}, L.index(nodes, :));
        for r = find(~keep(:)).'
            skipped(end+1) = struct('family', fam{k}, 'index', L.index(nodes(r), :), ...
                'reason', thin_reason(thin)); %#ok<AGROW>
        end
        nodes = nodes(keep);
    end
    for r = 1:numel(nodes)
        q = nodes(r);
        if ~C.seedable(q, f)
            why = 'not seedable';
            if ~C.in_range(q), why = ['outside native range: ' C.outside_reason{q}]; end
            skipped(end+1) = struct('family', fam{k}, 'index', L.index(q, :), ...
                'reason', why); %#ok<AGROW>
            continue
        end
        list(end+1) = struct('family', fam{k}, 'left_type', C.family_types{f, 1}, ...
            'right_type', C.family_types{f, 2}, 'node', q, 'index', L.index(q, :), ...
            'coordinate', L.coordinate(q, :), 'key', L.key{q}, ...
            'cell_left', C.cell(q, :), 'cell_right', C.mirror(q, :)); %#ok<AGROW>
    end
end
if strcmp(S.targets.policy, 'per_seed_centre') && ~isempty(list)
    % milestone 8 run option: per_seed's targets, the cold calls ordered
    % nearest d = 0 first (Euclidean in lattice coordinates); ties keep the
    % per_seed order (family, then lattice order)
    d = sqrt(sum(vertcat(list.coordinate).^2, 2));
    [~, o] = sortrows([d (1:numel(list)).']);
    list = list(o);
end
T.list = list;
T.skipped = skipped;
end

function [nodes, skipped] = per_seed(L, C, f, family)
% One node per (family, library cell): the anchor node when it is in the
% seed's group, else the nearest to the cell (ties within 1e-9 native
% spacings by lattice order). nodes: lattice order; skipped: the rest, in
% lattice order.
skipped = struct('family', {}, 'index', {}, 'reason', {});
cand = find(C.seedable(:, f));
nodes = zeros(0, 1);
if isempty(cand), return; end
[~, ~, g] = unique(C.cell(cand, :), 'rows');
d = sqrt(sum(C.offset(cand, :).^2, 2));
anchor = all(L.index(cand, :) == 0, 2);
pick = false(numel(cand), 1);
tgt = zeros(numel(cand), 1);                  % the chosen candidate of each seed
for s = 1:max(g)
    m = find(g == s);                          % ascending: lattice order
    a = m(anchor(m));
    if isempty(a)
        a = m(find(d(m) <= min(d(m)) + 1e-9, 1));
    end
    pick(a) = true;
    tgt(m) = a;
end
nodes = cand(pick);
r = find(~pick);
if isempty(r), return; end
idx = num2cell(L.index(cand(r), :), 2);
why = arrayfun(@(k) sprintf('per_seed: seed cell %s targeted at node %s', ...
    mat2str(C.cell(cand(k), :)), L.key{cand(tgt(k))}), r, 'UniformOutput', false);
skipped = struct('family', family, 'index', idx(:).', 'reason', why(:).');
end

function keep = thin_keep(S, family, index)
% True for every index row kept by targets.thin for this family (all rows
% when the family is not thinned); the anchor (index 0) is always kept.
keep = true(size(index, 1), 1);
if ~isfield(S.targets, 'thin') || ~any(strcmp(family, S.targets.thin.families)), return; end
st = S.targets.thin.stride;
o = S.targets.thin.offset;
if isempty(o), o = zeros(1, size(index, 2)); end
keep = all(mod(index - o, st) == 0, 2) | all(index == 0, 2);
end

function t = thin_reason(thin)
o = thin.offset;
if isempty(o), o = 0; end
t = sprintf('thinned: stride %d offset %s', thin.stride, mat2str(o));
end
