function rows = sa_t83_start_lookups(E)
%SA_T83_START_LOOKUPS  The lineage-start lookups (sa_continuation 'targets',
%   'next_start', 'unowned_roots', 'next_target', 'node_parent') against the
%   per-element versions they replace (the reference functions below are
%   the earlier per-element code, verbatim apart from their names). Zero bvp4c calls.
%   A synthetic 9 x 9 lattice (index -4..4 per axis; ring = the outer
%   band), 0-3 roots per node (some created by the per-root pass), and a
%   queue grown as 'expand' grows it: a lineage per start, successes
%   attached at neighbouring nodes (adopted or merged), random tried pairs
%   and edge counts, roots recorded without a slot. For the phases
%   interior, ring and all, search.ring_last true and false, a box limit and
%   none, and contexts with and without the node root-count vector: every
%   output is isequal (class and size included) to the reference's at every
%   step, through the no-slot branch to an empty next_start.
T = 'T83 lineage-start lookups equal the per-element versions';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.pkg);
rng(83, 'twister');
nchk = 0; bad = {}; nstart = 0; nsame = 0; npass = 0; nnoslot = 0;
for variant = 1:4
    ring_last = variant ~= 2;
    with_counts = variant >= 3;
    boxed = variant == 4;
    ctx = build(ring_last, with_counts, boxed);
    for ph = {'interior', 'ring', 'all'}
        c = ctx; c.Q.phase = ph{1};
        nloop = 0;
        while true
            nloop = nloop + 1;
            if nloop > 1000                                 % a fault fails here instead of hanging
                bad{end+1} = sprintf('v%d %s: no empty next_start after 1000 steps', variant, ph{1}); %#ok<AGROW>
                break
            end
            nchk = nchk + 1;
            tag = sprintf('v%d %s lin %d', variant, ph{1}, numel(c.Q.lin));
            if ~same(ref_targets(c), sa_continuation('targets', c)), bad{end+1} = [tag ' targets']; end %#ok<AGROW>
            [r0, q0, a0, f0] = ref_next_start(c);
            [r1, q1, a1, f1] = sa_continuation('next_start', c);
            if ~same({r0, q0, a0, f0}, {r1, q1, a1, f1}), bad{end+1} = [tag ' next_start']; end %#ok<AGROW>
            U = ref_unowned(c);
            if ~same(U, sa_continuation('unowned_roots', c)), bad{end+1} = [tag ' unowned_roots']; end %#ok<AGROW>
            for l = 1:numel(c.Q.lin)
                [qa, pa] = ref_next_target(c, l); [qb, pb] = sa_continuation('next_target', c, l);
                if ~same([qa pa], [qb pb]), bad{end+1} = sprintf('%s next_target %d', tag, l); end %#ok<AGROW>
            end
            for p = 1:c.Q.n
                [x0, y0, z0] = ref_node_parent(c, p); [x1, y1, z1] = sa_continuation('node_parent', c, p);
                if ~same({x0, y0, z0}, {x1, y1, z1}), bad{end+1} = sprintf('%s node_parent %d', tag, p); end %#ok<AGROW>
            end
            if isempty(r1), break; end
            if c.Q.slots_used >= c.S.continuation.max_lineages
                c.Q.no_slot = union(c.Q.no_slot, U);       % expand's no-slot branch
                nnoslot = nnoslot + 1;
                continue                                    % the next next_start returns nothing
            end
            nstart = nstart + 1;
            % a start at a node that already holds another lineage's root (key 1 compares the start
            % with those roots; it is 0 here as on every reachable state: lineage roots are owned and
            % owned roots never start, and one node's roots are distinct by membership)
            nsame = nsame + any(cellfun(@(l) ~isempty(l{q1}), {c.Q.lin.node_root}));
            c.Q = new_lineage(c.Q, q1, r1, a1, f1);
            l = numel(c.Q.lin);
            c.Q.lin(l).slot = true; c.Q.slots_used = c.Q.slots_used + 1;
            % successes at a few nodes near the start: unowned roots adopted, owned ones merged
            kq = c.Q.index(q1, :);
            for dd = [1 0; -1 0; 0 1; 0 -1; 1 1; -1 -1; 2 0].'
                p = sa_continuation('row', c.Q, kq + dd.');
                if p == 0 || rand < 0.3, continue; end
                ids = roots_of(c, p);
                if isempty(ids), continue; end
                k = randi(numel(ids));
                st = c.shards(c.L.key{p});
                c = sa_continuation('attach', c, l, p, ids{k}, st.roots(k).representative_attempt_id);
            end
            % random tried pairs and edge counts (next_target's filters)
            for t = 1:3
                q = randi(c.Q.n);
                kp = c.Q.index(q, :); kp(randi(2)) = kp(randi(2)) + 1;
                p = sa_continuation('row', c.Q, kp);
                if p > 0, c.Q.lin(randi(l)).tried(end+1, :) = [p q]; end
                c.Q.lin(randi(l)).edges_into(randi(c.Q.n)) = randi([1 4]);
            end
            if rand < 0.25 && ~isempty(U), c.Q.no_slot = union(c.Q.no_slot, U(randi(numel(U)))); end
        end
    end
end
npass = sum(cellfun(@(k) any(arrayfun(@(r) is_pass(r), ctx.shards(k).roots)), ctx.shards.keys));
rows = sa_check(rows, T, sprintf(['targets, next_start, unowned_roots, next_target and node_parent equal the ' ...
    'per-element versions at every step (%d states, 4 variants x 3 phases; %d lineage starts, %d of them ' ...
    'at a node already holding another lineage''s root; the no-slot branch %d times)'], nchk, nstart, nsame, nnoslot), ...
    isempty(bad) && nstart >= 40 && nsame >= 1 && nnoslot >= 6 && npass >= 1, strjoin(bad(1:min(end, 8)), '; '));
src = fileread(fullfile(E.pkg, 'sa_continuation.m'));
body = regexprep(src, '(?m)^\s*%[^\n]*$', '');               % comment lines out (no dot: line-bounded)
rows = sa_check(rows, T, 'sa_continuation.m passes no anonymous function over ctx (each call costs tens of ms at production size)', ...
    isempty(regexp(body, '@\([^)]*\)[^\n]*\<ctx\>', 'once')), '');
% the whole package: no anonymous function names ctx after its argument list (continuation lines joined,
% comment lines out); in sa_report.m the run context is called st, so st is checked there too
g = dir(fullfile(E.pkg, '**', '*.m'));
hits = {};
for k = 1:numel(g)
    f = fullfile(g(k).folder, g(k).name);
    b = regexprep(fileread(f), '\.\.\.[^\n]*\n', ' ');      % join continuation lines
    b = regexprep(b, '(?m)^\s*%[^\n]*$', '');                    % comment lines out (line-bounded)
    big = '\<ctx\>';
    if strcmp(g(k).name, 'sa_report.m'), big = '\<(ctx|st)\>'; end
    m = regexp(b, ['@\([^)]*\)[^\n]*' big], 'match');
    for q = 1:numel(m), hits{end+1} = sprintf('%s: %s', g(k).name, m{q}(1:min(end, 120))); end %#ok<AGROW>
end
rows = sa_check(rows, T, sprintf(['no anonymous function in the package (%d files) captures the run context ' ...
    '(ctx; st in sa_report.m)'], numel(g)), numel(g) > 50 && isempty(hits), strjoin(hits, ' || '));
end

% =========================================================================
% the synthetic state

function ctx = build(ring_last, with_counts, boxed)
[i, j] = ndgrid(-4:4, -4:4);
idx = sortrows([i(:) j(:)]);
n = size(idx, 1);
keys = arrayfun(@(k) sprintf('k:%d,%d', idx(k, 1), idx(k, 2)), (1:n).', 'UniformOutput', false);
L = struct('n_nodes', n, 'index', idx, 'key', {keys});
S = struct('search', struct('ring_last', ring_last), ...
    'continuation', struct('max_edges_per_node', 2, 'max_lineages', 5), ...
    'multiseed', struct('max_roots_per_node', 3));
Q = struct();
Q.n = n; Q.index = idx; Q.kmin = [-4 -4]; Q.kmax = [4 4]; Q.strides = [9 1];
if boxed, Q.box = all(abs(idx) <= 3, 2); else, Q.box = true(n, 1); end
Q.ring = any(abs(idx) >= 3, 2);
Q.anchor_node = find(all(idx == 0, 2));
Q.phase = ''; Q.pass = 1;
Q.lin = repmat(blank_lineage(n), 1, 0);
Q.owner_root = {}; Q.owner_lin = zeros(1, 0); Q.slots_used = 0; Q.no_slot = {};
shards = containers.Map('KeyType', 'char', 'ValueType', 'any');
nroots = zeros(n, 1);
fams = {'HSS_HSS', 'HSS_HPSS', 'HPSS_HSS'};
for q = 1:n
    nr = randi([0 3]);
    if q == Q.anchor_node, nr = 2; end
    if nr == 0, continue; end
    dn = randperm(6, nr);                                  % creation numbers, not in root order
    for r = 1:nr
        nd = randi([1 2]);
        route = 'cold'; if rand < 0.15, route = 'root_pass'; end
        d = struct('discovery_number', num2cell(1:nd), 'route', 'continuation', ...
            'seed_family', fams{randi(3)});
        d(1).discovery_number = dn(r); d(1).route = route;
        roots(r) = struct('root_id', sprintf('R%06d', r), 'created_discovery_number', dn(r), ...
            'representative_attempt_id', sprintf('A%06d', 1000 * q + r), 'discoveries', d); %#ok<AGROW>
    end
    shards(keys{q}) = struct('roots', roots(1:nr));
    nroots(q) = nr;
    clear roots
end
ctx = struct('S', S, 'L', L, 'Q', Q, 'shards', shards);
if with_counts, ctx.nroots = nroots; end
% the anchor lineage, as init: the anchor node's first root
ids = roots_of(ctx, Q.anchor_node);
st = shards(keys{Q.anchor_node});
ctx.Q = new_lineage(ctx.Q, Q.anchor_node, ids{1}, st.roots(1).representative_attempt_id, 'HSS_HSS');
ctx.Q.lin(1).kind = 'anchor';
end

function tf = is_pass(r)
d = r.discoveries([r.discoveries.discovery_number] == r.created_discovery_number);
tf = strcmp(d(1).route, 'root_pass');
end

function l = blank_lineage(n)
l = struct('id', 0, 'kind', '', 'start_node', 0, 'start_root', '', 'start_attempt', '', ...
    'family', '', 'calls', 0, 'slot', false, 'status', 'pending', 'merged_with', zeros(1, 0), ...
    'node_root', {repmat({''}, 1, n)}, 'node_attempt', {repmat({''}, 1, n)}, ...
    'expandable', false(1, n), 'edges_into', zeros(1, n), 'tried', zeros(0, 2));
end

function Q = new_lineage(Q, node, root, att, family)
l = blank_lineage(Q.n);
l.id = numel(Q.lin) + 1;
l.kind = 'start'; l.start_node = node; l.start_root = root; l.start_attempt = att;
l.family = family; l.status = 'active';
l.node_root{node} = root; l.node_attempt{node} = att; l.expandable(node) = true;
Q.lin(l.id) = l;
Q.owner_root{end+1} = root;
Q.owner_lin(end+1) = l.id;
end

function ids = roots_of(ctx, q)
key = ctx.L.key{q};
ids = {};
if isKey(ctx.shards, key)
    st = ctx.shards(key);
    ids = arrayfun(@(r) sprintf('%s/%s', key, r.root_id), st.roots, 'UniformOutput', false);
end
end

function tf = same(a, b)
tf = isequal(a, b) && strcmp(class(a), class(b)) && isequal(size(a), size(b));
if tf && iscell(a)
    tf = all(cellfun(@(x, y) strcmp(class(x), class(y)) && isequal(size(x), size(y)), a, b));
end
end

% =========================================================================
% reference: the per-element versions (the earlier sa_continuation.m)

function tf = in_targets(Q, S, q)
switch Q.phase
    case 'interior', tf = Q.box(q) && (~S.search.ring_last || ~Q.ring(q));
    case 'ring',     tf = Q.box(q) && Q.ring(q);
    otherwise,       tf = Q.box(q);
end
end

function T = ref_targets(ctx)
T = find(arrayfun(@(q) in_targets(ctx.Q, ctx.S, q), (1:ctx.Q.n).'));
end

function R = ref_unowned(ctx)
R = {};
T = ref_targets(ctx);
for q = T.'
    ids = roots_of(ctx, q);
    pr = pass_created(ctx, q);
    for k = 1:numel(ids)
        if sa_continuation('owner', ctx.Q, ids{k}) == 0 && ~pr(k), R{end+1} = ids{k}; end %#ok<AGROW>
    end
end
end

function [root, qbest, att, fam] = ref_next_start(ctx)
root = ''; qbest = 0; att = ''; fam = '';
Q = ctx.Q;
best = [];
T = ref_targets(ctx);
a = Q.index(Q.anchor_node, :);
for q = T.'
    key = ctx.L.key{q};
    if ~isKey(ctx.shards, key), continue; end
    st = ctx.shards(key);
    pr = pass_created(ctx, q);
    for r = 1:numel(st.roots)
        id = sprintf('%s/%s', key, st.roots(r).root_id);
        if sa_continuation('owner', Q, id) ~= 0 || any(strcmp(Q.no_slot, id)), continue; end
        if pr(r), continue; end
        same = any(arrayfun(@(l) strcmp(l.node_root{q}, id), Q.lin));
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

function pr = pass_created(ctx, q)
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

function [qbest, pbest] = ref_next_target(ctx, l)
Q = ctx.Q; S = ctx.S;
lin = Q.lin(l);
s = Q.index(lin.start_node, :);
qbest = 0; pbest = 0; best = [];
nc = size(Q.index, 2);
for q = 1:Q.n
    if ~in_targets(Q, S, q) || ~isempty(lin.node_root{q}) || ...
            lin.edges_into(q) >= S.continuation.max_edges_per_node || ...
            sa_continuation('n_roots', ctx, q) >= S.multiseed.max_roots_per_node
        continue
    end
    kq = Q.index(q, :);
    cand = zeros(0, 3);
    for c = 1:nc
        for sg = [-1 1]
            kp = kq; kp(c) = kp(c) + sg;
            p = sa_continuation('row', Q, kp);
            if p == 0 || isempty(lin.node_root{p}) || ~lin.expandable(p), continue; end
            if any(lin.tried(:, 1) == p & lin.tried(:, 2) == q), continue; end
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

function [root, att, lin] = ref_node_parent(ctx, p)
root = ''; att = ''; lin = 0;
ids = roots_of(ctx, p);
if isempty(ids), return; end
own = cellfun(@(r) sa_continuation('owner', ctx.Q, r), ids);
key = own; key(key == 0) = Inf;
[~, k] = min(key);
root = ids{k}; lin = own(k);
if lin > 0 && strcmp(ctx.Q.lin(lin).node_root{p}, root)
    att = ctx.Q.lin(lin).node_attempt{p};
else
    st = ctx.shards(ctx.L.key{p});
    rr = st.roots(strcmp(arrayfun(@(r) sprintf('%s/%s', ctx.L.key{p}, r.root_id), st.roots, ...
        'UniformOutput', false), root));
    att = rr.representative_attempt_id;
end
end
