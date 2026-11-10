function varargout = sa_record_call(varargin)
%SA_RECORD_CALL  Validity, root store and per-call record for one corrector call.
%
%   [ctx, A] = SA_RECORD_CALL(ctx, R, unit, tg, target, V, write)
%       Evaluates the validity boundary on R (or takes V, the re-derived
%       record of a replayed attempt), classifies the outcome, adds a valid
%       NODE return of any route (cold, continuation, rescue_*) to its node's
%       root-store shard (an import return too: route 'import', milestone
%       10; refinement returns and substeps are never
%       discoveries), and appends the attempt row. With write = true the
%       attempt file attempts/A######.mat is written atomically (returned
%       iterate x, y, yp, stats; options; warning; validity record; target;
%       work unit; provenance). Replayed attempts are never rewritten.
%   [ok, why] = SA_RECORD_CALL('budget', st)
%       Call and wall budget check (kept for older callers); the
%       executor uses sa_search('budget').
%   tf = SA_RECORD_CALL('cut_by_cap', kind, opts, solver, failing, reason, returned)
%       The "cut by cap" label of one
%       call. True iff kind is 'cold', the call ran tighter than the run's
%       solver (opts.nmax < solver.nmax or opts.mesh_reject <
%       solver.mesh_reject), and it was refused on a mesh bound: a RETURNED
%       call whose failing checks include convergence_proposed_mesh or mesh,
%       or a call whose corrector reason is the module's mesh blow-up (the
%       only test for a call that returned nothing, whose every check is
%       NaN and therefore listed as failing). Independent of the returned
%       mesh; a positivity-only refusal is not a cut; a SUPERSET of the
%       calls the cap changed ("refused on a mesh bound while the cap was
%       tighter", not "would have converged"). Outcome and verdict are
%       never changed by it.
%   Attempt-row fields: nmax (the call's own opts.nmax) and
%   cut_by_cap.
%   Row field max_proposed_mesh
%   = P, the largest mesh bvp4c proposed in the call (R.max_proposed_mesh: 0
%   when no proposal was made, NaN when nothing was recorded -- a stub
%   corrector, a legacy record, an interrupted call); the attempt file gains
%   mesh_trace (R.mesh_trace without the volatile per-process call token).
%   On every call (live and replay) a recorded
%   trace must be consistent, else SA:INTEGRITY "mesh trace inconsistent":
%   n_proposals > 0 -> P >= stats.nmeshpoints and the last proposal =
%   stats.nmeshpoints; n_proposals = 0 and returned -> stats.nmeshpoints =
%   N0 = m; returned -> m <= max(N0, P); at most one initial mesh; P equals
%   its proposals. Live label check: for a cold call run
%   tighter than the run's solver, cut_by_cap must equal the exact rule at
%   the call's own options, cut iff P > opts.nmax or (returned and m >=
%   opts.mesh_reject).
%
%   Light shards: between calls ctx.shards holds each node's store
%   light (sa_store_io 'light': no n_dense returned mesh / profile); before
%   the kernel add the node's store is made full again from its
%   discoveries' attempt files (sa_store_io 'full'), and after the add it
%   stays full until the executor's shard flush (write, or the replay
%   compare), which makes it light again. ctx.nroots(node) (the
%   node root-count vector, ctx.node_of: key -> node) is set here, the
%   only place a store gains a root.
%
%   Outcomes: accepted (valid node return, stored), invalid, no_return,
%   refined_valid (refinement call that passed), substep_valid /
%   substep_invalid (non-lattice coordinate, never stored), interrupted
%   (written by the resume reconciliation, sa_resume).

if ischar(varargin{1}) && strcmp(varargin{1}, 'budget')
    [varargout{1}, varargout{2}] = budget_ok(varargin{2});
    return
end
if ischar(varargin{1}) && strcmp(varargin{1}, 'cut_by_cap')
    varargout{1} = cut_rule(varargin{2:end});
    return
end
[varargout{1}, varargout{2}] = record_call(varargin{:});
end

function [ok, why] = budget_ok(st)
ok = true; why = '';
if st.calls >= st.S.budget.max_calls
    ok = false; why = 'budget_calls';
elseif toc(st.t_run) >= st.S.budget.max_wall_s
    ok = false; why = 'budget_wall';
end
end

function [ctx, A] = record_call(ctx, R, unit, tg, target, V, write)
M = ctx.M; S = ctx.S;
n = numel(ctx.attempts) + 1;
id = sprintf('A%06d', n);
if isempty(V)
    V = sa_validity('evaluate', M, R, ctx.theta, target, S);
end
A = sa_report('blank');
A.attempt_id = id;
A.kind = unit.kind;
A.parent_attempt = unit.parent_attempt;
A.target_number = unit.target_number;
A.family = tg.family;
A.target_key = tg.key;
A.index = tg.index;
A.coordinate = tg.coordinate;
A.seed_left = unit.seed_left;
A.seed_right = unit.seed_right;
A.corrector_class = R.corrector_class;
A.corrector_reason = R.corrector_reason;
A.corrector_nodal_min = R.corrector_nodal_min;
A.convergence_class = V.convergence_class;
A.wall_s = R.wall_s;
A.reltol = R.opts.reltol;
A.failing = strjoin(V.failing, ';');
A.flags = strjoin(V.flags, ';');
A.valid = V.valid;
A.file = fullfile('attempts', [id '.mat']);
v = @(name) value_of(V, name);
A.returned_mesh = v('mesh');
A.maxres = v('convergence_maxres');
A.positivity_min = v('positivity');
A.residual_scaled_sup = v('residual');
A.seam_concentration = v('seam_concentration');
A.seam_flux = v('seam_flux');
A.s_max = v('s_max');
if isfield(V.positivity, 'nodal_min')
    A.nodal_min = V.positivity.nodal_min;
    A.positivity_row = V.positivity.row;
    A.positivity_xi = V.positivity.xi;
    A.uniform_grid_min = V.positivity.uniform_grid_min;
end
if isfield(V.residual, 'scaled_rms')
    A.residual_scaled_rms = V.residual.scaled_rms;
    A.residual_unscaled_sup = V.residual.unscaled_sup;
end
if V.returned && v('finite') == 0
    try
        A.observed_family = shape_label(M, R.sol, ctx);
    catch ME
        A.observed_family = 'unlabelled';
        fprintf('%s: family label failed, recorded as unlabelled (%s: %s)\n', id, ME.identifier, ME.message);
    end
end
% ---- route fields --------------------------------------------------------
A.route = unit.route;
A.phase = unit.phase;
A.lineage = unit.lineage;
A.edge_id = unit.edge_id;
A.fraction = unit.fraction;
A.parent_root = unit.parent_root;
A.guess_source = unit.guess_source;
A.retry_of = unit.retry_of;
A.unit_key = sa_search('unit_key', unit);
A.positivity_sampled_min = v('positivity_sampled');
A.residual_bvp4c = v('residual_bvp4c');
% ---- cold-cap fields ------------------------------------------------------
A.nmax = R.opts.nmax;
A.cut_by_cap = ~V.valid && cut_rule(unit.kind, R.opts, S.solver, V.failing, R.corrector_reason, ...
    ~isempty(R.sol));
% ---- mesh-trace fields ----------------------------------------------------
tr = trace_of(R);
A.max_proposed_mesh = NaN;
if isfield(R, 'max_proposed_mesh'), A.max_proposed_mesh = double(R.max_proposed_mesh); end
why = trace_problem(R, tr, A.max_proposed_mesh, unit.kind, S.solver, A.cut_by_cap);
if ~isempty(why)
    try
        sa_journal('note', ctx.state_dir, sprintf('integrity mesh trace inconsistent at %s: %s', id, why));
    catch
    end
    error('SA:INTEGRITY', 'Integrity stop at attempt count %d: mesh trace inconsistent at %s: %s', ...
        n - 1, id, why);
end

node_return = unit.fraction == 1 && ~strcmp(unit.kind, 'refinement');
if isempty(R.sol)
    if node_return, A.outcome = 'no_return'; else, A.outcome = sub_or(unit, 'no_return'); end
elseif ~V.valid
    if node_return, A.outcome = 'invalid'; else, A.outcome = sub_or(unit, 'invalid'); end
elseif node_return
    store = shard(ctx, tg.key);              % full again (refilled from the attempt files)
    cand = sa_identity('candidate', M, R.sol, tg, id, ...
        nonempty(A.observed_family, 'unlabelled'), V, S, unit.route, ctx.K);
    [store, out] = sa_root_store('add', store, cand);
    ctx.shards(tg.key) = store;              % full until the executor's flush, then light
    ctx.dirty{end+1} = tg.key;
    if isfield(ctx, 'nroots') && isfield(ctx, 'node_of') && isKey(ctx.node_of, tg.key)
        ctx.nroots(ctx.node_of(tg.key)) = numel(store.roots);   % the node root-count vector
    end
    A.outcome = 'accepted';
    A.root_id = sprintf('%s/%s', tg.key, out.root_id);
    A.store_status = out.status;
    A.representative_distance = out.representative_distance;
elseif strcmp(unit.kind, 'refinement')
    A.outcome = 'refined_valid';
else
    A.outcome = 'substep_valid';
end
if write
    rec = struct('attempt', A, 'sol', R.sol, 'opts', R.opts, 'warning_id', R.warning_id, ...
        'warning_msg', R.warning_msg, 'validity', V, 'target', target, 'target_row', tg, ...
        'unit', unit, 'provenance', ctx.prov, 'mesh_trace', stored_trace(tr));
    sa_write('mat', fullfile(ctx.state_dir, A.file), rec);
end
ctx.attempts(n) = A;
ctx.calls = n;
end

function tf = cut_rule(kind, opts, solver, failing, reason, returned)
if ischar(failing), failing = strsplit(failing, ';'); end
tighter = opts.nmax < solver.nmax || opts.mesh_reject < solver.mesh_reject;
blowup = ischar(reason) && startsWith(strtrim(reason), 'mesh blow-up');
mesh_bound = blowup || (returned && any(ismember({'convergence_proposed_mesh', 'mesh'}, failing)));
tf = strcmp(kind, 'cold') && tighter && mesh_bound;
end

function tr = trace_of(R)
% the call's mesh trace; the "not recorded" trace when R carries none
if isfield(R, 'mesh_trace') && isstruct(R.mesh_trace) && isscalar(R.mesh_trace)
    tr = R.mesh_trace;
else
    tr = sa_mesh_trace('blank');
end
end

function tr = stored_trace(tr)
% the per-process call token is volatile (a counter of the recorder's resets
% in this MATLAB process): it is checked at 'take' and never stored, so
% records stay a pure function of the work
if isfield(tr, 'token'), tr = rmfield(tr, 'token'); end
end

function why = trace_problem(R, tr, P, kind, solver, cut)
% Mesh-trace consistency checks; '' when consistent or when
% nothing was recorded (P NaN)
why = '';
if ~isfield(tr, 'proposals') || ~isfield(tr, 'initial') || ~isfield(tr, 'n_initial')
    why = 'malformed mesh trace'; return
end
recorded = isfinite(tr.initial) || ~isempty(tr.proposals);
if ~recorded
    if ~isnan(P), why = sprintf('P %g with nothing recorded', P); end
    return
end
if ~isequaln(P, sa_mesh_trace('max_of', tr))
    why = sprintf('P %g differs from its proposals (%g)', P, sa_mesh_trace('max_of', tr)); return
end
if tr.n_initial > 1
    why = sprintf('%d initial meshes recorded in one call', tr.n_initial); return
end
returned = ~isempty(R.sol);
if returned
    m = numel(R.sol.x);
    ns = NaN;
    if isfield(R.sol, 'stats') && isstruct(R.sol.stats) && isfield(R.sol.stats, 'nmeshpoints')
        ns = double(R.sol.stats.nmeshpoints);
    end
    if ~isempty(tr.proposals)
        if ~(P >= ns && tr.proposals(end) == ns)
            why = sprintf('P %g, last proposal %g, stats.nmeshpoints %g', P, tr.proposals(end), ns); return
        end
    elseif ~(ns == tr.initial && tr.initial == m)
        why = sprintf('no proposal but stats.nmeshpoints %g, initial %g, returned mesh %d', ns, tr.initial, m); return
    end
    if ~(m <= max(tr.initial, P))
        why = sprintf('returned mesh %d above max(initial %g, P %g)', m, tr.initial, P); return
    end
end
o = R.opts;
if strcmp(kind, 'cold') && (o.nmax < solver.nmax || o.mesh_reject < solver.mesh_reject)
    exact = P > o.nmax || (returned && numel(R.sol.x) >= o.mesh_reject);
    if ~isequal(logical(cut), exact)
        why = sprintf(['cut_by_cap %d disagrees with the exact rule %d (P %g, nmax %g, returned %d, ' ...
            'mesh_reject %g)'], cut, exact, P, o.nmax, returned, o.mesh_reject); return
    end
end
end

function o = sub_or(unit, node_outcome)
if strcmp(unit.kind, 'refinement'), o = node_outcome; else, o = 'substep_invalid'; end
end

function store = shard(ctx, key)
% The node's store in full form: the executor keeps shards light between
% calls (sa_store_io 'light'); the kernel needs every discovery's returned
% mesh and profile, recomputed here from each discovery's attempt file.
% The function handle captures only the state folder (all that get_sol
% reads): a handle over ctx walks the whole context at every call, tens of
% ms each at production size; the same file is loaded the same way.
if isKey(ctx.shards, key)
    sd = struct('state_dir', ctx.state_dir);
    store = sa_store_io('full', ctx.shards(key), ctx.M, ctx.S, ctx.K, @(id) sa_search('get_sol', sd, id));
else
    store = sa_root_store('empty', ctx.K);
end
end

function x = value_of(V, name)
k = strcmp({V.checks.name}, name);
if any(k), x = V.checks(find(k, 1)).value; else, x = NaN; end
end

function s = nonempty(s, d)
if isempty(s), s = d; end
end

function lab = shape_label(M, sol, ctx)
% The module's family label. A module whose shape_label takes a fourth
% argument gets the run's field scale (ctx.field_scale: the anchor HSS
% values, the mode names' uniform-state reference); a three-argument one
% (a test stub) does not.
if abs(nargin(M.shape_label)) >= 4
    lab = M.shape_label(sol, ctx.theta, ctx.S.model_options, ctx.field_scale);
else
    lab = M.shape_label(sol, ctx.theta, ctx.S.model_options);
end
end
