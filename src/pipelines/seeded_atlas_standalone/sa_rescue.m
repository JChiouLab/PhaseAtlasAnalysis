function varargout = sa_rescue(action, varargin)
%SA_RESCUE  Rescue ladder for unreached nodes.
%
%   [ctx, added] = SA_RESCUE('run', ctx)
%       After the continuation queue of the pass empties: every node of the
%       phase's target set (lattice order) that holds no root, was not rescued
%       before, and has a cardinal neighbour holding a root, gets the ladder
%       in settings order, stopping at the first success:
%         diagonal     a diagonal-neighbour root, adaptive two-axis edge,
%                      h = continuation.max_halvings          (<= 2h+1 calls)
%         secant       collinear roots of ONE lineage at q-e (P1) and q-2e
%                      (P2): guess 2*deval(P1, x1) - deval(P2, x1) on P1's
%                      mesh x1 (the guess mesh of P1 under
%                      continuation.guess_mesh; default = P1's mesh), all
%                      states; one call per direction, in the
%                      order +axis1, -axis1, +axis2, -axis2 (P1 = q - e)
%         deep         RESUMES the refused cardinal
%                      adaptive edge into q with the largest last valid
%                      fraction (ties: the edge's lineage, 0 last; then the
%                      parent's lattice order; then edge id) from its last
%                      valid substep's stored iterate (the parent's when that
%                      fraction is 0), at pos = fraction_valid with step
%                      s_fail/2 (half the last refused step), refused at its
%                      third refusal (2 halvings; <= 10 calls at h = 4); the
%                      record keeps h = max_halvings + 2 and gains
%                      resumed_edge, start_fraction, start_step. Without a
%                      refused cardinal adaptive edge into q: the plain rule
%                      (the best cardinal parent, adaptive from the parent,
%                      h = max_halvings + 2)
%       Cap: rescue.max_calls_per_node over the whole ladder (each node is
%       rescued at most once per run). Each rung is a fresh attempt under the
%       full boundary; a rescued root joins its parent's lineage (an unowned
%       parent leaves it unowned, so it may start a lineage) and re-opens
%       continuation. added = number of nodes rescued.
%       rescue.ring = false skips the ring band (nodes
%       outside the library's native range): continuation only there.
%       A node whose ladder made NO call (e.g. secant
%       only, and no lineage holds roots at q-e and q-2e: cold and pass roots
%       are never a secant basis) is not marked rescue_done; it records
%       rescue_result 'no_basis' and stays eligible for the next closure
%       pass. At most one ladder WITH calls per node, as before.
%   g = SA_RESCUE('secant_guess', solP1, solP2 [, S])    the secant guess (test hook)

switch action
    case 'run',          [varargout{1}, varargout{2}] = run(varargin{:});
    case 'secant_guess', varargout{1} = secant_guess(varargin{:});
    otherwise, error('SA:RESCUE', 'Unknown sa_rescue action %s.', action);
end
end

function g = secant_guess(P1, P2, S)
% 2*deval(P1, x1) - deval(P2, x1) on P1's mesh x1; x1 =
% the guess mesh of P1 under continuation.guess_mesh (default 'parent').
if nargin < 3, S = []; end
g = sa_guess_mesh('secant', P1, P2, S);
end

function [ctx, added] = run(ctx)
S = ctx.S;
added = 0;
T = sa_continuation('targets', ctx);
nc = size(ctx.Q.index, 2);
% the rescue set is fixed when the continuation queue has emptied: nodes that
% gain a neighbour root during this pass wait for the next continuation pass
todo = [];
for q = T.'
    if ctx.Q.rescue_done(q) || sa_continuation('n_roots', ctx, q) > 0, continue; end
    if ctx.Q.ring(q) && ~S.rescue.ring, continue; end            % ring rescue off
    if ~isempty(card_roots(ctx, q, nc)), todo(end+1) = q; end %#ok<AGROW>
end
for q = todo
    Q = ctx.Q;
    if sa_continuation('n_roots', ctx, q) > 0, continue; end
    card = card_roots(ctx, q, nc);
    ctx.Q.rescue_done(q) = true;
    cap = S.rescue.max_calls_per_node;
    used = 0; ok = false; result = 'refused';
    for rg = S.rescue.rungs
        if used >= cap, break; end
        switch rg{1}
            case 'diagonal'
                [ctx, ok, n, lin] = rung_edge(ctx, q, diagonal_parent(ctx, q), ...
                    'rescue_diagonal', S.continuation.max_halvings, cap - used);
            case 'secant'
                [ctx, ok, n, lin] = rung_secant(ctx, q, cap - used);
            case 'deep'
                Er = resumed_edge(ctx, q);
                if isempty(Er)
                    [ctx, ok, n, lin] = rung_edge(ctx, q, deep_parent(ctx, q, card), ...
                        'rescue_deep', S.continuation.max_halvings + 2, cap - used);
                else
                    [ctx, ok, n, lin] = rung_resume(ctx, q, Er, cap - used);
                end
        end
        used = used + n;
        if ctx.stop, ctx.Q.rescue_result{q} = 'unfinished'; return; end
        if ok
            result = ['rescued_' rg{1}];
            added = added + 1;
            A = ctx.attempts(end);
            if lin > 0
                ctx = sa_continuation('attach', ctx, lin, q, A.root_id, A.attempt_id);
            end
            break
        end
    end
    if ~ok && used >= cap, result = 'refused_cap'; end
    if used == 0 && ~ok
        % A zero-call ladder leaves the node eligible (done is set
        % before the rungs, so the state every call of a ladder
        % sees is unchanged; it is undone here only when no call was made)
        ctx.Q.rescue_done(q) = false;
        result = 'no_basis';
    end
    ctx.Q.rescue_result{q} = result;
end
end

function card = card_roots(ctx, q, nc)
Q = ctx.Q; kq = Q.index(q, :);
card = [];
for c = 1:nc
    for sg = [-1 1]
        k = kq; k(c) = k(c) + sg;
        p = sa_continuation('row', Q, k);
        if p > 0 && sa_continuation('n_roots', ctx, p) > 0, card(end+1) = p; end %#ok<AGROW>
    end
end
end

function P = diagonal_parent(ctx, q)
% First diagonal neighbour holding a root, by (owning lineage, lattice order).
Q = ctx.Q; P = [];
kq = Q.index(q, :);
best = [];
if size(Q.index, 2) < 2, return; end
for s1 = [-1 1]
    for s2 = [-1 1]
        k = kq; k(1) = k(1) + s1; k(2) = k(2) + s2;
        p = sa_continuation('row', Q, k);
        if p == 0 || sa_continuation('n_roots', ctx, p) == 0, continue; end
        [root, att, lin] = sa_continuation('node_parent', ctx, p);
        lk = lin; if lk == 0, lk = Inf; end
        key = [lk, p];
        if isempty(best) || key(1) < best(1) || (key(1) == best(1) && key(2) < best(2))
            best = key; P = struct('p', p, 'root', root, 'att', att, 'lin', lin);
        end
    end
end
end

function P = deep_parent(ctx, q, card)
Q = ctx.Q; P = []; best = [];
E = Q.edges;
for p = card
    [root, att, lin] = sa_continuation('node_parent', ctx, p);
    m = [E.q] == q & [E.p] == p & strcmp({E.status}, 'refused');
    if any(m), fr = max([E(m).fraction_valid]); else, fr = -1; end
    lk = lin; if lk == 0, lk = Inf; end
    key = [-fr, lk, p];
    if isempty(best) || lexless(key, best)
        best = key; P = struct('p', p, 'root', root, 'att', att, 'lin', lin);
    end
end
end

function tf = lexless(a, b)
k = find(a ~= b, 1);
tf = ~isempty(k) && a(k) < b(k);
end

function [ctx, ok, n, lin] = rung_edge(ctx, q, P, kind, h, cap)
ok = false; n = 0; lin = 0;
if isempty(P), return; end
fam = family_of(ctx, P.root);
spec = struct('kind', kind, 'route', kind, 'rule', 'adaptive', 'h', h, 'lineage', P.lin, ...
    'p', P.p, 'q', q, 'parent_root', P.root, 'parent_attempt', P.att, 'family', fam, 'cap', cap);
[ctx, e] = sa_continuation('edge', ctx, spec);
n = e.calls; ok = strcmp(e.status, 'completed'); lin = P.lin;
end

function Er = resumed_edge(ctx, q)
% The refused cardinal adaptive edge into q with the largest
% fraction_valid; ties: the edge's lineage (0 last), the parent's lattice
% order, the edge id.
Er = [];
E = ctx.Q.edges;
if isempty(E), return; end
kq = ctx.Q.index(q, :);
best = [];
for k = find([E.q] == q & strcmp({E.status}, 'refused') & strcmp({E.rule}, 'adaptive') & [E.calls] > 0)
    e = E(k);
    if sum(abs(ctx.Q.index(e.p, :) - kq)) ~= 1, continue; end      % cardinal only
    lk = e.lineage; if lk == 0, lk = Inf; end
    key = [-e.fraction_valid, lk, e.p, e.id];
    if isempty(best) || lexless(key, best)
        best = key; Er = e;
    end
end
end

function [ctx, ok, n, lin] = rung_resume(ctx, q, Er, cap)
% The deep rung resumes the refused edge Er.
S = ctx.S;
ids = {ctx.attempts.attempt_id};
last = ctx.attempts(strcmp(ids, Er.attempts{end}));
pos = Er.fraction_valid;
s_fail = last.fraction - pos;
satt = Er.parent_attempt;
if pos > 0
    for j = numel(Er.attempts):-1:1
        a = ctx.attempts(strcmp(ids, Er.attempts{j}));
        if a.valid && a.fraction == pos, satt = a.attempt_id; break; end
    end
end
fam = family_of(ctx, Er.parent_root);
spec = struct('kind', 'rescue_deep', 'route', 'rescue_deep', 'rule', 'adaptive', ...
    'h', S.continuation.max_halvings + 2, 'lineage', Er.lineage, 'p', Er.p, 'q', q, ...
    'parent_root', Er.parent_root, 'parent_attempt', Er.parent_attempt, 'family', fam, ...
    'cap', cap, 'start_pos', pos, 'start_step', s_fail / 2, 'start_attempt', satt, ...
    'max_refusals', 2, 'resumed_edge', Er.id);
[ctx, e] = sa_continuation('edge', ctx, spec);
n = e.calls; ok = strcmp(e.status, 'completed'); lin = Er.lineage;
if lin > 0 && ~strcmp(ctx.Q.lin(lin).node_root{Er.p}, Er.parent_root), lin = 0; end
end

function [ctx, ok, n, lin] = rung_secant(ctx, q, cap)
ok = false; n = 0; lin = 0;
Q = ctx.Q; S = ctx.S;
kq = Q.index(q, :);
nc = size(Q.index, 2);
for c = 1:nc
    for sg = [1 -1]
        if n >= cap || ctx.stop, return; end
        e1 = zeros(1, nc); e1(c) = sg;
        p1 = sa_continuation('row', Q, kq - e1);
        p2 = sa_continuation('row', Q, kq - 2 * e1);
        if p1 == 0 || p2 == 0, continue; end
        l = find(arrayfun(@(L) ~isempty(L.node_root{p1}) && ~isempty(L.node_root{p2}), Q.lin), 1);
        if isempty(l), continue; end
        L = Q.lin(l);
        s1 = sa_search('get_sol', ctx, L.node_attempt{p1});
        s2 = sa_search('get_sol', ctx, L.node_attempt{p2});
        guess = secant_guess(s1, s2, S);
        e = struct('id', numel(ctx.Q.edges) + 1);
        tg = struct('family', L.family, 'left_type', '', 'right_type', '', 'node', q, ...
            'index', kq, 'coordinate', ctx.L.coordinate(q, :), 'key', ctx.L.key{q}, ...
            'cell_left', [NaN NaN], 'cell_right', [NaN NaN]);
        u = sa_search('unit');
        u.route = 'rescue_secant'; u.kind = 'rescue'; u.phase = Q.phase; u.pass = Q.pass;
        u.family = L.family; u.key = ctx.L.key{q}; u.index = kq; u.fraction = 1;
        u.coordinate = ctx.L.coordinate(q, :); u.lineage = l; u.edge_id = e.id;
        u.parent_attempt = L.node_attempt{p1}; u.parent_root = L.node_root{p1};
        u.guess_source = sprintf('secant 2*%s - %s%s', L.node_attempt{p1}, L.node_attempt{p2}, ...
            sa_guess_mesh('tag', S));
        u.seed_left = sprintf('secant %s %s', ctx.L.key{p1}, ctx.L.key{p2});
        u.seed_right = L.node_root{p2};
        [ctx, A] = sa_search('call', ctx, u, tg, guess, S.solver);
        E = struct('id', e.id, 'kind', 'rescue_secant', 'rule', 'secant', 'h', 0, 'lineage', l, ...
            'p', p1, 'q', q, 'parent_root', L.node_root{p1}, 'parent_attempt', L.node_attempt{p1}, ...
            'phase', Q.phase, 'pass', Q.pass, 'status', '', 'fraction_valid', 0, 'failing', '', ...
            'calls', 0, 'refusals', 0, 'attempts', {{}}, 'child_root', '', 'child_status', '', ...
            'd_parent_child', NaN, 'd_max_substep', NaN, 'large_change', false, 'resumed_edge', 0, ...
            'start_fraction', 0, 'start_step', 1, 'd_sub', zeros(1, 0), 's_sub', zeros(1, 0), ...
            'round', 0);   % edge fields and the pass round
        if isempty(A)
            E.status = 'unfinished'; ctx.Q.edges(E.id) = E; return
        end
        n = n + 1;
        ctx.Q.rescue_calls(q) = ctx.Q.rescue_calls(q) + 1;
        E.calls = 1; E.attempts = {A.attempt_id};
        if A.valid
            E.status = 'completed'; E.fraction_valid = 1; E.child_root = A.root_id;
            E.child_status = A.store_status;
            E.d_parent_child = sa_identity('profile_distance', ctx.M, s1, ...
                sa_search('get_sol', ctx, A.attempt_id), ctx.K);
            E.d_max_substep = E.d_parent_child;
            E.d_sub = E.d_parent_child; E.s_sub = 1;
            E.large_change = E.d_parent_child > S.identity.near_threshold(2);
            ctx.Q.edges(E.id) = E;
            ok = true; lin = l; return
        end
        E.status = 'refused'; E.refusals = 1;
        if isempty(A.failing), E.failing = A.outcome; else, E.failing = A.failing; end
        ctx.Q.edges(E.id) = E;
    end
end
end

function fam = family_of(ctx, root)
key = regexprep(root, '/R\d+$', '');
st = ctx.shards(key);
r = st.roots(strcmp(arrayfun(@(x) sprintf('%s/%s', key, x.root_id), st.roots, ...
    'UniformOutput', false), root));
fam = r.discoveries(1).seed_family;
end
