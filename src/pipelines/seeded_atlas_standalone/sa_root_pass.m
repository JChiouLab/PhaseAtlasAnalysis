function varargout = sa_root_pass(action, varargin)
%SA_ROOT_PASS  Per-root continuation pass.
%
%   ctx = SA_ROOT_PASS('run', ctx)
%       Rounds of the pass for the current phase and search pass, placed
%       after the continuation lineages and before rescue (sa_search). Every
%       distinct root r at a node p of the continuation region (Q.box) takes
%       one adaptive step (h = root_pass.max_halvings, at most 2h+1 calls)
%       to each cardinal neighbour q of p in the phase's target set where no
%       matching root is known; whatever valid root the step lands on enters
%       the store under the usual identity rule (scaled_rms_v1 at
%       identity.similarity_tolerance) and is new only if it differs from
%       every stored root at q. Edges go through sa_continuation('edge')
%       with kind = route = 'root_pass', rule 'adaptive', lineage 0, the
%       root's representative attempt as the parent iterate (never
%       re-solved), the library family of the parent's observed mode name
%       as the unit family (family_lookup below), and the edge field
%       round = k.
%       Round 1 starts from every distinct root in the region; round k+1
%       from the NEW roots of round k only (a child matching a stored root
%       is not stepped from again). The pass stops when a round's target set
%       is empty, or after root_pass.max_rounds rounds. Reaching
%       root_pass.max_calls_per_round (checked before every edge) ends the
%       pass for the current search pass (status root_pass_cap; the untried
%       targets are reported as unattempted, never refused; they are posed
%       again only if a later closure pass runs). A round's target list is
%       enumerated once, at the round's start (as the per-pass list), and
%       every target of it is posed; across rounds and search passes rule
%       (a) keeps a posed connection from being posed again.
%       multiseed.max_roots_per_node is not
%       applied (the pass is bounded per root). Pass-created roots carry
%       lineage 0, never seed a lineage (sa_continuation next_start) and are
%       never a secant basis (sa_rescue needs one lineage at q-e and q-2e).
%   [T, D] = SA_ROOT_PASS('targets', ctx [, start])
%       The target list (pure function of the store and the edge table):
%       (r, q) is a target iff r is a distinct root at node p in Q.box (or
%       in the cell list `start` of root ids, when given), q is a cardinal
%       neighbour of p in the phase's target set (sa_continuation
%       'targets'), and NEITHER (a) an edge of any route or status with
%       calls > 0 already went from r into q, NOR (b) r is the completed
%       child of an edge from q. Order: p in lattice order, then the root's
%       creation (discovery) number at p, then the directions +axis1,
%       -axis1, +axis2, -axis2. T: struct array (p, q, root, attempt,
%       family, dn, dir). D: counts (start_roots, start_nodes, candidates,
%       known_a, known_b_only, posed).
%   R = SA_ROOT_PASS('blank_rounds')     1x0 rounds table (Q.root_pass)
%   L = SA_ROOT_PASS('report_lines', ctx)  summary.txt lines

switch action
    case 'run',          varargout{1} = run(varargin{:});
    case 'targets',      [varargout{1}, varargout{2}] = targets(varargin{:});
    case 'blank_rounds', varargout{1} = repmat(blank_row(), 1, 0);
    case 'report_lines', varargout{1} = report_lines(varargin{:});
    otherwise, error('SA:ROOT_PASS', 'Unknown sa_root_pass action %s.', action);
end
end

% =========================================================================

function r = blank_row()
r = struct('search_pass', 0, 'phase', '', 'round', 0, 'start_roots', 0, 'candidates', 0, ...
    'known_a', 0, 'known_b_only', 0, 'targets', 0, 'edges', 0, 'calls', 0, ...
    'completed', 0, 'new_roots', 0, 'matched', 0, 'refused', 0, 'unattempted', 0, 'max_calls_edge', 0, ...
    'status', '');
end

function ctx = run(ctx)
S = ctx.S;
P = S.root_pass;
start = [];                                   % round 1: every root in the region
for k = 1:P.max_rounds
    [T, D] = targets(ctx, start);
    row = blank_row();
    row.search_pass = ctx.Q.pass; row.phase = ctx.Q.phase; row.round = k;
    row.start_roots = D.start_roots; row.candidates = D.candidates; row.known_a = D.known_a;
    row.known_b_only = D.known_b_only; row.targets = numel(T);
    if isempty(T)
        row.status = 'empty';
        ctx.Q.root_pass(end+1) = row;
        return
    end
    new = {};
    capped = false;
    for t = 1:numel(T)
        if row.calls >= P.max_calls_per_round
            row.unattempted = numel(T) - t + 1;
            capped = true;
            break
        end
        tg = T(t);
        spec = struct('kind', 'root_pass', 'route', 'root_pass', 'rule', 'adaptive', ...
            'h', P.max_halvings, 'lineage', 0, 'p', tg.p, 'q', tg.q, 'parent_root', tg.root, ...
            'parent_attempt', tg.attempt, 'family', tg.family, 'cap', Inf, 'round', k);
        [ctx, e] = sa_continuation('edge', ctx, spec);
        row.calls = row.calls + e.calls;
        if e.calls > 0, row.edges = row.edges + 1; end
        row.max_calls_edge = max(row.max_calls_edge, e.calls);
        if ctx.stop
            row.status = 'unfinished';
            row.unattempted = numel(T) - t;
            ctx.Q.root_pass(end+1) = row;
            return
        end
        switch e.status
            case 'completed'
                row.completed = row.completed + 1;
                if strcmp(e.child_status, 'new_root')
                    new{end+1} = e.child_root; %#ok<AGROW>
                else
                    row.matched = row.matched + 1;
                end
            case 'refused'
                row.refused = row.refused + 1;
        end
    end
    row.new_roots = numel(new);
    if capped
        row.status = 'root_pass_cap';
        ctx.Q.root_pass(end+1) = row;
        return
    end
    if k == P.max_rounds && ~isempty(new)
        row.status = 'max_rounds';
    else
        row.status = 'done';
    end
    ctx.Q.root_pass(end+1) = row;
    if isempty(new)
        % the next round's target set is empty by construction (round k+1
        % starts from the new roots of round k only); record it
        r2 = blank_row(); r2.search_pass = ctx.Q.pass; r2.phase = ctx.Q.phase; r2.round = k + 1;
        r2.status = 'empty';
        if k < P.max_rounds, ctx.Q.root_pass(end+1) = r2; end
        return
    end
    start = new;
end
end

% =========================================================================
% targets

function [T, D] = targets(ctx, start)
if nargin < 2, start = []; end
Q = ctx.Q;
T = struct('p', {}, 'q', {}, 'root', {}, 'attempt', {}, 'family', {}, 'dn', {}, 'dir', {});
D = struct('start_roots', 0, 'start_nodes', 0, 'candidates', 0, 'known_a', 0, 'known_b_only', 0, 'posed', 0);
inT = false(Q.n, 1);
inT(sa_continuation('targets', ctx)) = true;
nc = size(Q.index, 2);
dirs = zeros(2 * nc, nc);
for c = 1:nc
    dirs(2 * c - 1, c) = 1; dirs(2 * c, c) = -1;      % +axis1, -axis1, +axis2, -axis2
end
E = Q.edges;
[KA, KB] = known_sets(E);
fam_of = family_lookup(ctx);
nodes = 0;
for p = find(Q.box(:)).'
    key = ctx.L.key{p};
    if ~isKey(ctx.shards, key), continue; end
    st = ctx.shards(key);
    if isempty(st.roots), continue; end
    [~, o] = sort([st.roots.created_discovery_number]);
    hit = false;
    for r = o
        id = sprintf('%s/%s', key, st.roots(r).root_id);
        if ~isempty(start) && ~any(strcmp(start, id)), continue; end
        hit = true;
        D.start_roots = D.start_roots + 1;
        att = st.roots(r).representative_attempt_id;
        for d = 1:size(dirs, 1)
            q = sa_continuation('row', Q, Q.index(p, :) + dirs(d, :));
            if q == 0 || ~inT(q), continue; end
            D.candidates = D.candidates + 1;
            a = isKey(KA, sprintf('%s|%d', id, q));
            b = isKey(KB, sprintf('%s|%d', id, q));
            if a, D.known_a = D.known_a + 1; elseif b, D.known_b_only = D.known_b_only + 1; end
            if a || b, continue; end
            T(end+1) = struct('p', p, 'q', q, 'root', id, 'attempt', att, 'family', fam_of(st.roots(r)), ...
                'dn', st.roots(r).created_discovery_number, 'dir', d); %#ok<AGROW>
        end
    end
    nodes = nodes + hit;
end
D.start_nodes = nodes;
D.posed = numel(T);
end

function [KA, KB] = known_sets(E)
% (a) parent_root|q of every edge with calls > 0 (any route, any status);
% (b) child_root|p of every completed edge (the child at E.q, its parent
% node E.p).
KA = containers.Map('KeyType', 'char', 'ValueType', 'logical');
KB = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for k = 1:numel(E)
    if E(k).calls > 0 && ~isempty(E(k).parent_root)
        KA(sprintf('%s|%d', E(k).parent_root, E(k).q)) = true;
    end
    if strcmp(E(k).status, 'completed') && ~isempty(E(k).child_root)
        KB(sprintf('%s|%d', E(k).child_root, E(k).p)) = true;
    end
end
end

function fn = family_lookup(ctx)
% The library family of the parent's observed mode name (its
% representative attempt's shape label, mapped by sa_mode_name
% 'library_family' with the module's mode_families: in both models the flat
% half 0, a pole peak of mode 1 and the middle peak 2- have a library
% type); the creating discovery's seed family when the name maps to
% no library family (a seam peak, 2+, 3+/-, ..., 'unlabelled') or the
% attempt carries none. The family is the work unit's family (unit key,
% attempt row, the child discovery's seed_family); the edge's guess is the
% parent's own iterate whatever the family.
ids = containers.Map('KeyType', 'char', 'ValueType', 'double');
A = ctx.attempts;
for k = 1:numel(A), ids(A(k).attempt_id) = k; end
if isfield(ctx, 'M'), M = ctx.M; else, M = sa_model(ctx.S.model); end   % a context rebuilt from a folder
types = M.library_types;
map = [];
if isfield(M, 'mode_families'), map = M.mode_families; end
fn = @(root) fam(root, A, ids, types, map);
end

function f = fam(root, A, ids, types, map)
f = '';
if isKey(ids, root.representative_attempt_id)
    f = sa_mode_name('library_family', A(ids(root.representative_attempt_id)).observed_family, types, map);
end
if isempty(f)
    d = root.discoveries([root.discoveries.discovery_number] == root.created_discovery_number);
    f = d(1).seed_family;
end
end

% =========================================================================
% reporting

function L = report_lines(ctx)
L = {};
Q = ctx.Q;
if ~isfield(Q, 'root_pass'), return; end
S = ctx.S;
en = isfield(S, 'root_pass') && S.root_pass.enabled;
R = Q.root_pass;
E = Q.edges;
if isempty(E), pe = E; else, pe = E(strcmp({E.kind}, 'root_pass')); end
if ~en && isempty(R), return; end
nn = sum([R.new_roots]);
L{end+1} = sprintf(['Per-root pass (root_pass.enabled %d, h = %d, max rounds %d): %d round row(s); targets %d; ' ...
    'edges %d (completed %d, refused %d, unfinished %d); calls %d, max per edge %d; new roots %d, matched %d; ' ...
    'calls per new root %.3g'], en, S.root_pass.max_halvings, S.root_pass.max_rounds, numel(R), ...
    sum([R.targets]), numel(pe), sum(strcmp({pe.status}, 'completed')), sum(strcmp({pe.status}, 'refused')), ...
    sum(strcmp({pe.status}, 'unfinished')), sum([pe.calls]), max([0 pe.calls]), nn, sum([R.matched]), ...
    sum([pe.calls]) / max(nn, 1));
for k = 1:numel(R)
    r = R(k);
    L{end+1} = sprintf(['  pass %d %s round %d: start roots %d, candidates %d, known (a) %d, (b) only %d, ' ...
        'targets %d, calls %d, new %d, matched %d, refused %d, unattempted %d; %s'], ...
        r.search_pass, r.phase, r.round, r.start_roots, r.candidates, r.known_a, r.known_b_only, r.targets, ...
        r.calls, r.new_roots, r.matched, r.refused, r.unattempted, r.status); %#ok<AGROW>
end
end
