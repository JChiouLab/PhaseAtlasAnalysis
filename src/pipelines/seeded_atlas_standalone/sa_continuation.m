function varargout = sa_continuation(action, varargin)
%SA_CONTINUATION  Continuation route.
%
%   Q = SA_CONTINUATION('init', ctx)
%       Scheduler state: lattice geometry, box, ring, lineages (the anchor
%       lineage from the first accepted cold or import root at the anchor
%       node, in attempt order; milestone 10), root
%       ownership, edge table, rescue table. Q holds no clock and no handle.
%   ctx = SA_CONTINUATION('expand', ctx)
%       One continuation pass of the current phase: every existing lineage in
%       id order, then (continuation.starts = 'anchor_and_cold') new lineages
%       from accepted roots not yet in a lineage, ordered by (distinct by
%       membership from every lineage root at their node, index distance from
%       the anchor, lattice order, discovery order). Lineages are
%       expanded one after another (start order first).
%   [ctx, e] = SA_CONTINUATION('edge', ctx, spec)
%       One edge from a parent iterate at node p toward node q (also used by
%       the diagonal and deep rescue rungs): 'adaptive' or 'myxo_restart'.
%   L = SA_CONTINUATION('report_lines', ctx)     summary section
%   N = SA_CONTINUATION('node_classes', ctx)     per box node class
%   r = SA_CONTINUATION('row', Q, index)         lattice row of an index (0 = off)
%   n = SA_CONTINUATION('n_roots', ctx, q)       distinct roots stored at node q
%                                                (ctx.nroots(q) when present)
%   o = SA_CONTINUATION('owner', Q, root_id)     owning lineage (0 = none)
%   [root, att] = SA_CONTINUATION('node_parent', ctx, p)   preferred root at p
%   ctx = SA_CONTINUATION('attach', ctx, lin, q, root, att)
%   T = SA_CONTINUATION('targets', ctx)          the phase's target nodes
%   [root, q, att, fam] = SA_CONTINUATION('next_start', ctx)   next lineage start
%   R = SA_CONTINUATION('unowned_roots', ctx)    unowned non-pass roots at the targets
%   [q, p] = SA_CONTINUATION('next_target', ctx, lin)   lineage lin's next edge p -> q
%       (the last three are what 'expand' uses; exposed for the equivalence tests)
%   Cost: these helpers build their lookups (target mask, owned and no-slot
%   id sets, the lineage x node root table) once per call and pass no
%   anonymous function over the context (one costs tens of ms per call at
%   production size).
%
%   Lineage expansion: the next target q is the eligible node with the
%   smallest (L1 index distance from the lineage start, lattice order).
%   Eligible: q in the phase's target set, no success of this lineage at q,
%   fewer than continuation.max_edges_per_node edges of this lineage into q,
%   fewer than multiseed.max_roots_per_node roots stored at q, and an untried
%   cardinal neighbour p holding an expandable root of this lineage. Parent
%   preference: the axis-1 neighbour nearer the start, then axis 2, then the
%   farther ones (myxo's column-then-rows tree when every step succeeds). A
%   successful edge landing on a root owned by an earlier lineage records a
%   merge and is not expanded again; one landing on an unowned root (a cold
%   root) adopts it (a call-free merge for that root: it never consumes a
%   lineage slot). Non-anchor lineages: a slot is consumed at the first call;
%   at most continuation.max_lineages slots, continuation.max_calls_per_lineage
%   calls each.
%
%   Adaptive step rule: s = 1 (one lattice
%   step), r = 0. Refusal -> s = s/2, r = r+1; success -> advance s,
%   s = min(2s, remaining). The edge completes at remaining = 0 and is refused
%   when r > h; at most 2h+1 calls. myxo_restart: level l = 0..h restarts
%   from the parent with 2^l equal substeps (<= 2^(h+1)-1 calls).
%   Every substep passes the full boundary at its own coordinate before it
%   seeds the next; substeps are never stored. Each edge records the
%   scaled_rms_v1 distance parent -> child and the largest step-to-step
%   distance; either above identity.near_threshold(2) sets large_change
%   (report only, never connectivity).
%   Guess: every guess built from a stored iterate (the parent,
%   or the last valid substep) goes through sa_guess_mesh('from', sol, S);
%   continuation.guess_mesh = 'parent' passes the parent's own mesh.
%   continuation.region = 'cross' restricts the target set to
%   sa_modes('cross') (the lattice lines through the central node);
%   An adaptive spec may carry a start state (start_pos, start_step,
%   start_attempt: the deep rung resumes a refused edge) and max_refusals (2:
%   refused at the third refusal; the record keeps spec.h); every edge
%   records resumed_edge, start_fraction, start_step and, for edge_large_change, the
%   step-to-step distances d_sub with their fraction steps s_sub.
%   The per-root pass (sa_root_pass)
%   runs its edges through 'edge' with kind/route 'root_pass' (unit kind
%   'continuation', lineage 0, never counted as rescue calls); every edge
%   records the pass round ('round', 0 outside the pass); next_start skips
%   pass-created roots (a root whose creating discovery has route
%   'root_pass'), so no lineage opens from a root the pass steps from;
%   Q.root_pass holds the pass's rounds table.
%   Family fence (settings continuation.families; the run context carries
%   the list in force as ctx.fence, {} or absent = no fence): a root whose
%   observed family (its representative call's label) is outside the fence
%   is stored and reported, but a lineage landing on it adopts it without
%   making it expandable (counted: expansions withheld), next_start and
%   unowned_roots skip it (lineage starts withheld), the anchor lineage
%   starts from the first accepted cold or import root at the anchor node
%   that is inside the fence, and node_parent (the rescue parents) offers
%   only roots inside it. Q.fence (only when a fence is in force) holds the
%   list and the counts; a run without a fence keeps its queue state
%   unchanged.

switch action
    case 'init',         varargout{1} = init(varargin{:});
    case 'expand',       varargout{1} = expand(varargin{:});
    case 'edge',         [varargout{1}, varargout{2}] = run_edge(varargin{:});
    case 'report_lines', varargout{1} = report_lines(varargin{:});
    case 'node_classes', varargout{1} = node_classes(varargin{:});
    case 'row',          varargout{1} = row_of(varargin{:});
    case 'n_roots',      varargout{1} = n_roots(varargin{:});
    case 'owner',        varargout{1} = owner_of(varargin{:});
    case 'node_parent',  [varargout{1}, varargout{2}, varargout{3}] = node_parent(varargin{:});
    case 'attach',       varargout{1} = attach(varargin{:});
    case 'targets',      varargout{1} = phase_targets(varargin{:});
    case 'next_start',   [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = next_start(varargin{:});
    case 'unowned_roots', varargout{1} = unowned_roots(varargin{:});
    case 'next_target',  [varargout{1:max(nargout, 2)}] = next_target(varargin{:});
    case 'record_skips', varargout{1} = record_skips(varargin{:});
    otherwise, error('SA:CONTINUATION', 'Unknown sa_continuation action %s.', action);
end
end

% =========================================================================
% state

function Q = init(ctx)
L = ctx.L; S = ctx.S;
Q = struct();
Q.n = L.n_nodes;
Q.index = L.index;
Q.kmin = L.index_min; Q.kmax = L.index_max;
sizes = Q.kmax - Q.kmin + 1;
Q.strides = fliplr(cumprod([1, fliplr(sizes(2:end))]));
b = S.continuation.box;
if strcmp(S.continuation.region, 'cross')
    Q.box = sa_modes('cross', L);
elseif isempty(b)
    Q.box = true(Q.n, 1);
else
    Q.box = all(abs(L.index) <= b, 2);
end
Q.ring = ~ctx.C.in_range(:);
Q.anchor_node = find(all(L.index == 0, 2), 1);
Q.phase = ''; Q.pass = 0; Q.phases_run = {};
Q.lin = repmat(blank_lineage(Q.n), 1, 0);
Q.owner_root = {}; Q.owner_lin = zeros(1, 0);
Q.edges = repmat(blank_edge(), 1, 0);
Q.slots_used = 0;
Q.no_slot = {};
Q.rescue_done = false(Q.n, 1);
Q.rescue_calls = zeros(Q.n, 1);
Q.rescue_result = repmat({''}, Q.n, 1);
Q.root_pass = sa_root_pass('blank_rounds');      % the pass's rounds table
if knowledge_on(ctx)
    % source-edge skips (milestone 10; only when import.source_edges is not
    % 'none', so a run without them keeps its queue state unchanged)
    Q.source_skips = struct('category', {}, 'lineage', {}, 'round', {}, 'p', {}, 'q', {}, 'root', {});
    Q.source_skip_keys = {};
end
if fence_on(ctx)
    Q.fence = struct('labels', {reshape(ctx.fence, 1, [])}, 'expansions_withheld', 0, ...
        'starts_withheld', {{}}, 'pass_starts_withheld', 0, 'rescue_withheld', zeros(1, 0));
end
% the anchor lineage: the first accepted cold or import root at the anchor
% node, in attempt order (cold calls come first, so with cold on it is the
% first accepted cold root, as before); with a fence, the first one inside it
A = ctx.attempts;
k = find((strcmp({A.kind}, 'cold') | strcmp({A.kind}, 'import')) & strcmp({A.outcome}, 'accepted') & ...
    strcmp({A.target_key}, L.key{Q.anchor_node}), 1);
if ~isempty(k) && fence_on(ctx)
    cand = find((strcmp({A.kind}, 'cold') | strcmp({A.kind}, 'import')) & strcmp({A.outcome}, 'accepted') & ...
        strcmp({A.target_key}, L.key{Q.anchor_node}));
    k = [];
    for c = cand
        if root_in_fence(ctx, A(c).root_id), k = c; break; end
        if Q.box(Q.anchor_node), Q.fence.starts_withheld = union(Q.fence.starts_withheld, {A(c).root_id}); end
    end
end
if ~isempty(k) && Q.box(Q.anchor_node)
    Q = new_lineage(Q, 'anchor', Q.anchor_node, A(k).root_id, A(k).attempt_id, A(k).family);
end
end

function l = blank_lineage(n)
l = struct('id', 0, 'kind', '', 'start_node', 0, 'start_root', '', 'start_attempt', '', ...
    'family', '', 'calls', 0, 'slot', false, 'status', 'pending', 'merged_with', zeros(1, 0), ...
    'node_root', {repmat({''}, 1, n)}, 'node_attempt', {repmat({''}, 1, n)}, ...
    'expandable', false(1, n), 'edges_into', zeros(1, n), 'tried', zeros(0, 2));
end

function e = blank_edge()
e = struct('id', 0, 'kind', '', 'rule', '', 'h', 0, 'lineage', 0, 'p', 0, 'q', 0, ...
    'parent_root', '', 'parent_attempt', '', 'phase', '', 'pass', 0, 'status', '', ...
    'fraction_valid', 0, 'failing', '', 'calls', 0, 'refusals', 0, 'attempts', {{}}, ...
    'child_root', '', 'child_status', '', 'd_parent_child', NaN, 'd_max_substep', NaN, ...
    'large_change', false, 'resumed_edge', 0, 'start_fraction', 0, 'start_step', 1, ...
    'd_sub', zeros(1, 0), 's_sub', zeros(1, 0), 'round', 0);
end

function Q = new_lineage(Q, kind, node, root, att, family)
l = blank_lineage(Q.n);
l.id = numel(Q.lin) + 1;
l.kind = kind; l.start_node = node; l.start_root = root; l.start_attempt = att;
l.family = family; l.status = 'active';
l.node_root{node} = root; l.node_attempt{node} = att; l.expandable(node) = true;
Q.lin(l.id) = l;
Q = set_owner(Q, root, l.id);
end

function Q = set_owner(Q, root, lin)
Q.owner_root{end+1} = root;
Q.owner_lin(end+1) = lin;
end

function o = owner_of(Q, root)
k = find(strcmp(Q.owner_root, root), 1);
if isempty(k), o = 0; else, o = Q.owner_lin(k); end
end

function r = row_of(Q, idx)
if any(idx < Q.kmin) || any(idx > Q.kmax), r = 0; return; end
r = 1 + sum((idx - Q.kmin) .* Q.strides);
end

function n = n_roots(ctx, q)
% the node root-count vector when the context carries it (the executor's
% ctx; equal to the shards' counts, checked at every checkpoint), else the
% shards (a context built from a result or a run folder)
if isfield(ctx, 'nroots'), n = ctx.nroots(q); return; end
key = ctx.L.key{q};
if isKey(ctx.shards, key), n = numel(ctx.shards(key).roots); else, n = 0; end
end

function ids = roots_at(ctx, q)
key = ctx.L.key{q};
ids = {};
if isKey(ctx.shards, key)
    st = ctx.shards(key);
    ids = arrayfun(@(r) sprintf('%s/%s', key, r.root_id), st.roots, 'UniformOutput', false);
end
end

function m = target_mask(Q, S)
% The phase's target set as an n x 1 logical: interior = box nodes (outside
% the ring when search.ring_last), ring = box nodes in the ring, else the box
switch Q.phase
    case 'interior', m = Q.box(:) & (~S.search.ring_last | ~Q.ring(:));
    case 'ring',     m = Q.box(:) & Q.ring(:);
    otherwise,       m = Q.box(:);
end
m = logical(m);
end

function T = phase_targets(ctx)
% the phase's target nodes in lattice order (a column); no closure over ctx:
% an anonymous function capturing the whole context costs tens of ms per call
T = find(target_mask(ctx.Q, ctx.S));
end

function m = id_set(ids)
% a lookup set of root ids (isKey), built once per call
m = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for k = 1:numel(ids), m(ids{k}) = true; end
end

% =========================================================================
% expansion

function ctx = expand(ctx)
S = ctx.S;
for l = 1:numel(ctx.Q.lin)
    ctx = expand_lineage(ctx, l);
    if ctx.stop, return; end
end
if ~strcmp(S.continuation.starts, 'anchor_and_cold'), return; end
while true
    [root, q, att, fam] = next_start(ctx);
    if isempty(root)
        if isfield(ctx.Q, 'fence')                 % the unowned roots the fence keeps from starting
            ctx.Q.fence.starts_withheld = union(ctx.Q.fence.starts_withheld, fenced_unowned(ctx));
        end
        break
    end
    if ctx.Q.slots_used >= S.continuation.max_lineages
        % no slot: every remaining unowned root is recorded, none started
        ctx.Q.no_slot = union(ctx.Q.no_slot, unowned_roots(ctx));
        break
    end
    ctx.Q = new_lineage(ctx.Q, 'start', q, root, att, fam);
    ctx = expand_lineage(ctx, numel(ctx.Q.lin));
    if ctx.stop, return; end
end
end

function R = unowned_roots(ctx)
R = {};
T = phase_targets(ctx);
owned = id_set(ctx.Q.owner_root);              % owner_of(Q, id) ~= 0 (owner_lin >= 1)
fon = fence_on(ctx);
for q = T.'
    ids = roots_at(ctx, q);
    pr = pass_created(ctx, q);
    for k = 1:numel(ids)
        if fon && ~root_in_fence(ctx, ids{k}), continue; end     % outside the fence: never a start
        if ~isKey(owned, ids{k}) && ~pr(k), R{end+1} = ids{k}; end %#ok<AGROW>
    end
end
end

function R = fenced_unowned(ctx)
% unowned non-pass roots at the targets that the fence keeps from starting
R = {};
T = phase_targets(ctx);
owned = id_set(ctx.Q.owner_root);
for q = T.'
    ids = roots_at(ctx, q);
    pr = pass_created(ctx, q);
    for k = 1:numel(ids)
        if ~isKey(owned, ids{k}) && ~pr(k) && ~root_in_fence(ctx, ids{k}), R{end+1} = ids{k}; end %#ok<AGROW>
    end
end
end

function tf = fence_on(ctx)
% a family fence is in force (ctx.fence: the observed labels inside it)
tf = isfield(ctx, 'fence') && iscell(ctx.fence) && ~isempty(ctx.fence);
end

function tf = att_in_fence(ctx, att)
% the observed family of attempt att (ids are A###### = the row number)
A = ctx.attempts;
k = sscanf(att, 'A%d');
lab = '';
if isscalar(k) && k >= 1 && k <= numel(A) && strcmp(A(k).attempt_id, att)
    lab = A(k).observed_family;
else
    j = find(strcmp({A.attempt_id}, att), 1, 'last');
    if ~isempty(j), lab = A(j).observed_family; end
end
tf = any(strcmp(lab, ctx.fence));
end

function tf = root_in_fence(ctx, root)
% a root (<key>/R######) is inside the fence when its representative call's
% observed family is
key = regexprep(root, '/R\d+$', '');
tf = false;
if ~isKey(ctx.shards, key), return; end
st = ctx.shards(key);
for r = 1:numel(st.roots)
    if strcmp(sprintf('%s/%s', key, st.roots(r).root_id), root)
        tf = att_in_fence(ctx, st.roots(r).representative_attempt_id);
        return
    end
end
end

function [root, qbest, att, fam] = next_start(ctx)
% Unowned accepted roots in the phase's target set (any route, M-C).
root = ''; qbest = 0; att = ''; fam = '';
Q = ctx.Q;
best = [];
T = phase_targets(ctx);
a = Q.index(Q.anchor_node, :);
owned = id_set(Q.owner_root);                  % owner_of(Q, id) ~= 0 (owner_lin >= 1)
noslot = id_set(Q.no_slot);
NR = lineage_roots(Q);                         % lineage x node root ids
fon = fence_on(ctx);
for q = T.'
    key = ctx.L.key{q};
    if ~isKey(ctx.shards, key), continue; end
    st = ctx.shards(key);
    pr = pass_created(ctx, q);
    for r = 1:numel(st.roots)
        id = sprintf('%s/%s', key, st.roots(r).root_id);
        if isKey(owned, id) || isKey(noslot, id), continue; end
        if pr(r), continue; end                      % pass roots never seed a lineage
        if fon && ~att_in_fence(ctx, st.roots(r).representative_attempt_id), continue; end   % outside the fence
        % key 1: distinct by membership from every lineage root at q. A
        % redundant guard: the roots of one node's store are distinct by
        % membership by construction (each discovery joins exactly one
        % root; root ids are never reused or merged), and every lineage
        % root is owned, so an unowned root is never one of them and the
        % key is 0 on every reachable state
        same = any(strcmp(NR(:, q), id));
        d1 = sum(abs(Q.index(q, :) - a));
        dn = st.roots(r).created_discovery_number;
        k = [double(same), d1, q, dn];
        if isempty(best) || lexless(k, best)
            best = k; root = id; qbest = q;
            att = st.roots(r).representative_attempt_id;
            fam = st.roots(r).discoveries(1).seed_family;
        end
    end
end
end

function tf = lexless(a, b)
k = find(a ~= b, 1);
tf = ~isempty(k) && a(k) < b(k);
end

function NR = lineage_roots(Q)
% Every lineage's node_root as one row (numel(Q.lin) x Q.n cell), built once
% per call instead of an arrayfun over Q.lin for every root
if isempty(Q.lin), NR = cell(0, Q.n); else, NR = vertcat(Q.lin.node_root); end
end

function pr = pass_created(ctx, q)
% Per root of node q's shard, true when its creating discovery (the one
% numbered created_discovery_number) came from the per-root pass.
key = ctx.L.key{q};
pr = false(1, 0);
if ~isKey(ctx.shards, key), return; end
st = ctx.shards(key);
pr = false(1, numel(st.roots));
for r = 1:numel(st.roots)
    d = st.roots(r).discoveries;
    c = d([d.discovery_number] == st.roots(r).created_discovery_number);
    pr(r) = ~isempty(c) && strcmp(c(1).route, 'root_pass');
end
end

function ctx = expand_lineage(ctx, l)
S = ctx.S;
while true
    Q = ctx.Q;
    lin = Q.lin(l);
    if ~strcmp(lin.kind, 'anchor') && lin.calls >= S.continuation.max_calls_per_lineage
        ctx.Q.lin(l).status = 'lineage_cap';
        return
    end
    [q, p, skips] = next_target(ctx, l);
    if ~isempty(skips), ctx.Q = record_skips(ctx.Q, skips); end
    if q == 0
        ctx.Q.lin(l).status = 'done';
        return
    end
    if ~strcmp(lin.kind, 'anchor') && ~lin.slot
        ctx.Q.lin(l).slot = true;
        ctx.Q.slots_used = ctx.Q.slots_used + 1;
    end
    spec = struct('kind', 'continuation', 'route', 'continuation', ...
        'rule', S.continuation.step_rule, 'h', S.continuation.max_halvings, ...
        'lineage', l, 'p', p, 'q', q, 'parent_root', lin.node_root{p}, ...
        'parent_attempt', lin.node_attempt{p}, 'family', lin.family, 'cap', Inf);
    [ctx, e] = run_edge(ctx, spec);
    if e.calls > 0
        ctx.Q.lin(l).edges_into(q) = ctx.Q.lin(l).edges_into(q) + 1;
    end
    ctx.Q.lin(l).tried(end+1, :) = [p q];
    if strcmp(e.status, 'completed')
        ctx = attach(ctx, l, q, e.child_root, e.attempts{end});
    end
    if ctx.stop, return; end
end
end

function [qbest, pbest, skips] = next_target(ctx, l)
% skips (milestone 10): the candidate edges p -> q left out on source-edge
% knowledge (sa_import 'knowledge'): the lineage root at p maps to a source
% root whose edge into q the source refused with h >= continuation
% .max_halvings ('lineage_same_or_stronger'), or, under 'all', with a
% smaller h ('lineage_weaker')
Q = ctx.Q; S = ctx.S;
lin = Q.lin(l);
skips = struct('category', {}, 'lineage', {}, 'round', {}, 'p', {}, 'q', {}, 'root', {});
K = [];
if knowledge_on(ctx), K = ctx.ik; end
s = Q.index(lin.start_node, :);
qbest = 0; pbest = 0; best = [];
nc = size(Q.index, 2);
inT = target_mask(Q, S);                     % the target set, once per call
for q = 1:Q.n
    if ~inT(q) || ~isempty(lin.node_root{q}) || ...
            lin.edges_into(q) >= S.continuation.max_edges_per_node || ...
            n_roots(ctx, q) >= S.multiseed.max_roots_per_node
        continue
    end
    kq = Q.index(q, :);
    cand = zeros(0, 3);                       % [rank, p, 0]
    for c = 1:nc
        for sg = [-1 1]
            kp = kq; kp(c) = kp(c) + sg;
            p = row_of(Q, kp);
            if p == 0 || isempty(lin.node_root{p}) || ~lin.expandable(p), continue; end
            if any(lin.tried(:, 1) == p & lin.tried(:, 2) == q), continue; end
            if ~isempty(K)
                kk = sprintf('%s|%d', lin.node_root{p}, q);
                cat = '';
                if isKey(K.lin_ss, kk), cat = 'lineage_same_or_stronger';
                elseif strcmp(K.mode, 'all') && isKey(K.lin_weak, kk), cat = 'lineage_weaker';
                end
                if ~isempty(cat)
                    skips(end+1) = struct('category', cat, 'lineage', l, 'round', 0, 'p', p, 'q', q, ...
                        'root', lin.node_root{p}); %#ok<AGROW>
                    continue
                end
            end
            nearer = abs(kp(c) - s(c)) < abs(kq(c) - s(c));
            rank = c + nc * ~nearer;
            cand(end+1, :) = [rank, p, 0]; %#ok<AGROW>
        end
    end
    if isempty(cand), continue; end
    k = [sum(abs(kq - s)), q];
    if isempty(best) || lexless(k, best)
        best = k; qbest = q;
        cand = sortrows(cand, [1 2]);
        pbest = cand(1, 2);
    end
end
end

function tf = knowledge_on(ctx)
tf = isfield(ctx, 'ik') && isstruct(ctx.ik) && isfield(ctx.ik, 'active') && ctx.ik.active;
end

function Q = record_skips(Q, skips)
% each skipped (category, lineage or round, p, q, root) once
if ~isfield(Q, 'source_skips'), return; end
for k = 1:numel(skips)
    s = skips(k);
    key = sprintf('%s|%d|%d|%d|%d|%s', s.category, s.lineage, s.round, s.p, s.q, s.root);
    if any(strcmp(Q.source_skip_keys, key)), continue; end
    Q.source_skip_keys{end+1} = key;
    Q.source_skips(end+1) = s;
end
end

function ctx = attach(ctx, l, q, root, att)
% The lineage's success at q (one per lineage and node).
o = owner_of(ctx.Q, root);
ctx.Q.lin(l).node_root{q} = root;
ctx.Q.lin(l).node_attempt{q} = att;
if o == 0
    ctx.Q = set_owner(ctx.Q, root, l);
    if fence_on(ctx) && ~root_in_fence(ctx, root)
        % outside the fence: adopted and kept, never expanded
        ctx.Q.lin(l).expandable(q) = false;
        ctx.Q.fence.expansions_withheld = ctx.Q.fence.expansions_withheld + 1;
    else
        ctx.Q.lin(l).expandable(q) = true;
    end
elseif o ~= l
    ctx.Q.lin(l).merged_with = unique([ctx.Q.lin(l).merged_with, o]);
    ctx.Q.lin(l).expandable(q) = false;
end
end

function [root, att, lin] = node_parent(ctx, p)
% The preferred root at node p and the attempt whose iterate represents it:
% lineage-owned roots first (lineage id order), then unowned ones; with a
% fence, only roots inside it (none: root '').
root = ''; att = ''; lin = 0;
ids = roots_at(ctx, p);
if fence_on(ctx)
    keep = false(size(ids));                 % a loop: no closure over ctx
    for k = 1:numel(ids), keep(k) = root_in_fence(ctx, ids{k}); end
    ids = ids(keep);
end
if isempty(ids), return; end
own = zeros(size(ids));                      % a loop: no closure over ctx
for k = 1:numel(ids), own(k) = owner_of(ctx.Q, ids{k}); end
key = own; key(key == 0) = Inf;
[~, k] = min(key);
root = ids{k}; lin = own(k);
if lin > 0 && strcmp(ctx.Q.lin(lin).node_root{p}, root)
    att = ctx.Q.lin(lin).node_attempt{p};
else
    pk = ctx.L.key{p};
    st = ctx.shards(pk);
    rr = st.roots(strcmp(arrayfun(@(r) sprintf('%s/%s', pk, r.root_id), st.roots, ...
        'UniformOutput', false), root));
    att = rr.representative_attempt_id;
end
end

% =========================================================================
% one edge

function [ctx, e] = run_edge(ctx, spec)
S = ctx.S; Q = ctx.Q;
e = blank_edge();
e.id = numel(Q.edges) + 1;
e.kind = spec.kind; e.rule = spec.rule; e.h = spec.h; e.lineage = spec.lineage;
e.p = spec.p; e.q = spec.q; e.parent_root = spec.parent_root;
e.parent_attempt = spec.parent_attempt; e.phase = Q.phase; e.pass = Q.pass;
if isfield(spec, 'round') && ~isempty(spec.round), e.round = spec.round; end   % pass round
psol = sa_search('get_sol', ctx, spec.parent_attempt);
kp = Q.index(spec.p, :); kq = Q.index(spec.q, :);
thr = S.identity.near_threshold(2);
prev = psol; dmax = 0;
last_fail = '';
switch spec.rule
    case 'adaptive'
        pos = 0; s = 1; r = 0;
        guess = sa_guess_mesh('from', psol, S); gatt = spec.parent_attempt;
        hmax = spec.h;
        if isfield(spec, 'max_refusals') && ~isempty(spec.max_refusals), hmax = spec.max_refusals; end
        if isfield(spec, 'start_pos') && ~isempty(spec.start_pos)
            % resume a refused edge's adaptive state
            pos = spec.start_pos; s = spec.start_step;
            e.resumed_edge = spec.resumed_edge; e.start_fraction = pos; e.start_step = s;
            if ~isempty(spec.start_attempt) && ~strcmp(spec.start_attempt, spec.parent_attempt)
                prev = sa_search('get_sol', ctx, spec.start_attempt);
                guess = sa_guess_mesh('from', prev, S); gatt = spec.start_attempt;
            end
        end
        fprev = pos;
        while true
            f = pos + s;
            [ctx, A, R, why] = edge_call(ctx, spec, e, f, kp, kq, guess, gatt);
            if ~isempty(why), e.status = why; break; end
            e.calls = e.calls + 1; e.attempts{end+1} = A.attempt_id;
            if A.valid
                dd = sa_identity('profile_distance', ctx.M, prev, R.sol, ctx.K);
                dmax = max(dmax, dd);
                e.d_sub(end+1) = dd; e.s_sub(end+1) = f - fprev; fprev = f;
                prev = R.sol;
                pos = f;
                guess = sa_guess_mesh('from', R.sol, S); gatt = A.attempt_id;
                if pos == 1
                    e.status = 'completed'; e.child_root = A.root_id;
                    e.child_status = A.store_status;
                    break
                end
                s = min(2 * s, 1 - pos);
            else
                last_fail = fail_text(A);
                s = s / 2; r = r + 1;
                if r > hmax, e.status = 'refused'; break; end
            end
        end
        e.refusals = r; e.fraction_valid = pos;
    case 'myxo_restart'
        best = 0; r = 0;
        for level = 0:spec.h
            nsub = 2^level;
            guess = sa_guess_mesh('from', psol, S); gatt = spec.parent_attempt;
            prev = psol; ok = true;
            for j = 1:nsub
                f = j / nsub;
                [ctx, A, R, why] = edge_call(ctx, spec, e, f, kp, kq, guess, gatt);
                if ~isempty(why), e.status = why; ok = false; break; end
                e.calls = e.calls + 1; e.attempts{end+1} = A.attempt_id;
                if ~A.valid
                    last_fail = fail_text(A); ok = false; r = r + 1; break
                end
                dd = sa_identity('profile_distance', ctx.M, prev, R.sol, ctx.K);
                dmax = max(dmax, dd);
                if j == 1, e.d_sub = zeros(1, 0); e.s_sub = zeros(1, 0); end
                e.d_sub(end+1) = dd; e.s_sub(end+1) = 1 / nsub;
                prev = R.sol; best = max(best, f);
                guess = sa_guess_mesh('from', R.sol, S); gatt = A.attempt_id;
            end
            if ~isempty(e.status), break; end
            if ok
                e.status = 'completed'; e.child_root = A.root_id; e.child_status = A.store_status;
                best = 1;
                break
            end
        end
        if isempty(e.status), e.status = 'refused'; end
        e.refusals = r; e.fraction_valid = best;
    otherwise
        error('SA:CONTINUATION', 'Unknown step rule %s.', spec.rule);
end
if ~strcmp(e.status, 'completed'), e.failing = last_fail; end
e.d_max_substep = dmax;
if strcmp(e.status, 'completed')
    e.d_parent_child = sa_identity('profile_distance', ctx.M, psol, prev, ctx.K);
end
e.large_change = e.d_max_substep > thr || e.d_parent_child > thr;
ctx.Q.edges(e.id) = e;
end

function t = fail_text(A)
if any(strcmp(A.outcome, {'no_return'})) || (startsWith(A.outcome, 'substep') && ...
        isempty(A.failing) && ~A.valid)
    t = 'no_return';
else
    t = A.failing;
end
if isempty(t), t = A.outcome; end
end

function [ctx, A, R, why] = edge_call(ctx, spec, e, f, kp, kq, guess, gatt)
% Caps (never interrupted records: the scheduler counts its own units),
% then the executor.
A = []; R = []; why = '';
S = ctx.S;
if spec.lineage > 0 && strcmp(spec.kind, 'continuation') && ...
        ~strcmp(ctx.Q.lin(spec.lineage).kind, 'anchor') && ...
        ctx.Q.lin(spec.lineage).calls >= S.continuation.max_calls_per_lineage
    why = 'lineage_cap'; return
end
if e.calls >= spec.cap
    why = 'rescue_cap'; return
end
idx = kp + f * (kq - kp);
coord = S.lattice.anchor + S.lattice.step .* idx;      % never accumulated
q = spec.q;
if f == 1
    key = ctx.L.key{q};
else
    key = sprintf('%s<-%s@%s', ctx.L.key{q}, ctx.L.key{spec.p}, num2str(f, 17));
end
tg = struct('family', spec.family, 'left_type', '', 'right_type', '', 'node', q, ...
    'index', idx, 'coordinate', coord, 'key', key, 'cell_left', [NaN NaN], ...
    'cell_right', [NaN NaN]);
u = sa_search('unit');
u.route = spec.route;
if any(strcmp(spec.route, {'continuation', 'root_pass'})), u.kind = 'continuation'; else, u.kind = 'rescue'; end
u.phase = ctx.Q.phase; u.pass = ctx.Q.pass; u.family = spec.family;
u.key = key; u.index = idx; u.fraction = f; u.coordinate = coord;
u.lineage = spec.lineage; u.edge_id = e.id; u.parent_attempt = gatt;
u.parent_root = spec.parent_root; u.guess_source = ['attempt ' gatt sa_guess_mesh('tag', S)];
u.seed_left = sprintf('%s %s', spec.rule, ctx.L.key{spec.p});
u.seed_right = spec.parent_root;
[ctx, A, R] = sa_search('call', ctx, u, tg, guess, S.solver);
if isempty(A)
    why = 'unfinished'; return
end
if spec.lineage > 0 && strcmp(spec.kind, 'continuation')
    ctx.Q.lin(spec.lineage).calls = ctx.Q.lin(spec.lineage).calls + 1;
end
if startsWith(spec.route, 'rescue_')                 % rescue routes only (not the pass)
    ctx.Q.rescue_calls(q) = ctx.Q.rescue_calls(q) + 1;
end
end

% =========================================================================
% reporting

function N = node_classes(ctx)
% One class per box node: reached | unfinished | refused | unattempted |
% policy_skipped (a partition the tests assert).
Q = ctx.Q;
N = struct('node', {}, 'key', {}, 'ring', {}, 'class', {}, 'routes', {}, 'n_roots', {}, ...
    'reasons', {});
E = Q.edges;
Eq = [E.q]; Ec = [E.calls];          % once, not per node
for q = find(Q.box).'
    cls = 'unattempted'; routes = {}; reasons = {};
    nr = n_roots(ctx, q);
    if nr > 0
        cls = 'reached';
        st = ctx.shards(ctx.L.key{q});
        for r = 1:numel(st.roots)
            routes = union(routes, {st.roots(r).discoveries.route});
        end
        routes = reshape(routes, 1, []);
    else
        into = E(Eq == q & Ec > 0);
        if any(strcmp({into.status}, 'unfinished'))
            cls = 'unfinished';
        elseif ~isempty(into)
            cls = 'refused';
            reasons = unique({into.failing});
        elseif isfield(Q, 'source_skips') && any([Q.source_skips.q] == q)
            cls = 'refused';                         % milestone 10: only skipped edges
            reasons = {'refused in the source run (not posed here)'};
        elseif ~any(strcmp(Q.phases_run, phase_of(Q, ctx.S, q)))
            cls = 'policy_skipped';
            reasons = {'phase not run'};
        end
    end
    N(end+1) = struct('node', q, 'key', ctx.L.key{q}, 'ring', Q.ring(q), 'class', cls, ...
        'routes', {routes}, 'n_roots', nr, 'reasons', {reasons}); %#ok<AGROW>
end
end

function ph = phase_of(Q, S, q)
if ~S.search.ring_last, ph = 'all';
elseif Q.ring(q), ph = 'ring';
else, ph = 'interior';
end
end

function L = report_lines(ctx)
S = ctx.S; Q = ctx.Q;
L = {''};
L{end+1} = sprintf(['Search: routes cold=%d continuation=%d rescue=%d; step rule %s, ' ...
    'h = %d (max calls per edge 2h+1 = %d); phases run: %s'], S.routes.cold, ...
    S.routes.continuation, S.routes.rescue, S.continuation.step_rule, S.continuation.max_halvings, ...
    2 * S.continuation.max_halvings + 1, strjoin(Q.phases_run, ', '));
L{end+1} = sprintf('  route settings: rescue.rungs {%s}; rescue.ring %d; root_pass.enabled %d', ...
    strjoin(S.rescue.rungs, ','), S.rescue.ring, getf_pass(S));
if ~strcmp(S.continuation.guess_mesh, 'parent')
    L{end+1} = sprintf('Continuation guess mesh: %s, %d nodes', S.continuation.guess_mesh, ...
        S.continuation.guess_nodes);
end
N = node_classes(ctx);
for part = {'interior', 'ring'}
    if strcmp(part{1}, 'interior'), m = ~[N.ring]; else, m = [N.ring]; end
    if ~any(m), continue; end
    n = N(m);
    cls = {'reached', 'refused', 'unattempted', 'policy_skipped', 'unfinished'};
    cnt = cellfun(@(c) sum(strcmp({n.class}, c)), cls);
    if strcmp(part{1}, 'ring'), tag = ' (ours-only, no reference)'; else, tag = ''; end
    L{end+1} = sprintf('Node classes, %s box nodes%s: %d; %s', part{1}, tag, numel(n), ...
        strjoin(cellfun(@(c, k) sprintf('%s %d', c, k), cls, num2cell(cnt), ...
        'UniformOutput', false), ', ')); %#ok<AGROW>
    rr = n(strcmp({n.class}, 'reached'));
    rt = {};
    for k = 1:numel(rr), rt = [rt, reshape(rr(k).routes, 1, [])]; end %#ok<AGROW>
    [u, ~, j] = unique(rt);
    if ~isempty(u)
        L{end+1} = sprintf('  reached nodes holding a root of route: %s', strjoin(cellfun(@(a, b) ...
            sprintf('%s %d', a, b), u, num2cell(accumarray(j(:), 1).'), 'UniformOutput', false), ', ')); %#ok<AGROW>
    end
    L{end+1} = sprintf('  nodes with >= 2 distinct roots: %d', sum([n.n_roots] >= 2)); %#ok<AGROW>
end
E = Q.edges;
if ~isempty(E)
    kinds = unique({E.kind});
    for k = 1:numel(kinds)
        e = E(strcmp({E.kind}, kinds{k}));
        L{end+1} = sprintf(['Edges %s: %d (completed %d, refused %d, unfinished %d, caps %d); calls ' ...
            '%d, max per edge %d; large_change %d'], kinds{k}, numel(e), ...
            sum(strcmp({e.status}, 'completed')), sum(strcmp({e.status}, 'refused')), ...
            sum(strcmp({e.status}, 'unfinished')), sum(endsWith({e.status}, '_cap')), ...
            sum([e.calls]), max([e.calls]), sum([e.large_change])); %#ok<AGROW>
        rf = e(strcmp({e.status}, 'refused'));
        if ~isempty(rf)
            [u, ~, j] = unique({rf.failing});
            L{end+1} = sprintf('  refused: last valid fraction min/median/max %.4g/%.4g/%.4g; last failing: %s', ...
                min([rf.fraction_valid]), median([rf.fraction_valid]), max([rf.fraction_valid]), ...
                strjoin(cellfun(@(a, b) sprintf('%s %d', a, b), u, num2cell(accumarray(j(:), 1).'), ...
                'UniformOutput', false), '; ')); %#ok<AGROW>
        end
    end
end
if ~isempty(Q.lin)
    na = Q.lin(~strcmp({Q.lin.kind}, 'anchor'));
    L{end+1} = sprintf(['Lineages: %d (anchor %d); non-anchor slots used %d of %d; calls per ' ...
        'lineage %s; merges %d; unowned roots without a slot %d'], numel(Q.lin), ...
        sum(strcmp({Q.lin.kind}, 'anchor')), Q.slots_used, S.continuation.max_lineages, ...
        mat2str([Q.lin.calls]), sum(arrayfun(@(l) numel(l.merged_with), Q.lin)), numel(Q.no_slot));
    if ~isempty(na)
        L{end+1} = sprintf('  non-anchor lineage statuses: %s', strjoin({na.status}, ',')); %#ok<NASGU>
    end
end
if any(Q.rescue_done)
    L{end+1} = sprintf('Rescue: %d node(s) tried, %d calls, results: %s', sum(Q.rescue_done), ...
        sum(Q.rescue_calls), strjoin(Q.rescue_result(Q.rescue_done).', ','));
end
nb = sum(strcmp(Q.rescue_result, 'no_basis') & ~Q.rescue_done);
if nb > 0 || (S.routes.rescue && isequal(S.rescue.rungs, {'secant'}))
    L{end+1} = sprintf(['Rescue no_basis: %d node(s) whose ladder made no call (no lineage pair at ' ...
        'q-e, q-2e for secant); not marked done, eligible again in a later closure pass'], nb);
end
if isfield(Q, 'root_pass')
    L = [L, sa_root_pass('report_lines', ctx)];
end
if isfield(Q, 'fence')
    L = [L, fence_lines(ctx)];
end
end

function L = fence_lines(ctx)
% The fence section: the roots stored outside it (by label), the
% expansions, lineage starts and per-root pass starts it withheld
F = ctx.Q.fence;
labs = {};
keys = ctx.shards.keys;
n = 0;
for k = 1:numel(keys)
    st = ctx.shards(keys{k});
    for r = 1:numel(st.roots)
        n = n + 1;
        a = st.roots(r).representative_attempt_id;
        if ~att_in_fence(ctx, a)
            j = sscanf(a, 'A%d');
            if isscalar(j) && j <= numel(ctx.attempts), labs{end+1} = ctx.attempts(j).observed_family; end %#ok<AGROW>
        end
    end
end
L = {sprintf('Fence (continuation.families): {%s}', strjoin(F.labels, ', '))};
by = '';
if ~isempty(labs)
    [u, ~, j] = unique(labs);
    c = accumarray(j(:), 1).';
    [~, o] = sort(-c);
    o = o(1:min(end, 12));
    by = [' (' strjoin(cellfun(@(a, b) sprintf('%s %d', a, b), u(o), num2cell(c(o)), 'UniformOutput', false), ...
        ', ') tern(numel(u) > 12, ', ...', '') ')'];
end
L{end+1} = sprintf('  roots stored outside the fence: %d of %d%s', numel(labs), n, by);
L{end+1} = sprintf(['  withheld by the fence: lineage expansions %d (roots adopted, not expanded); lineage ' ...
    'starts %d (roots); per-root pass starts %d (root-rounds)'], F.expansions_withheld, ...
    numel(F.starts_withheld), F.pass_starts_withheld);
if isfield(F, 'rescue_withheld')
    Q = ctx.Q;
    nb = find(strcmp(Q.rescue_result(:).', 'no_basis') & ~Q.rescue_done(:).');
    fw = intersect(nb, F.rescue_withheld);
    L{end+1} = sprintf(['  rescue: %d node(s) whose ladder had no parent or basis inside the fence (recorded no_basis; ' ...
        'genuine no_basis %d)'], numel(fw), numel(nb) - numel(fw));
end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function v = getf_pass(S)
v = isfield(S, 'root_pass') && isfield(S.root_pass, 'enabled') && S.root_pass.enabled;
end
