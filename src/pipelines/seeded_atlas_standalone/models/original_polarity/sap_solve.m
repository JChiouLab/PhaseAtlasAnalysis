function R = sap_solve(theta, guess, target, solver, mo, corrector)
%SAP_SOLVE  One coupled-BVP solver call, returned in the pipeline's common form.
%
%   R = SAP_SOLVE(theta, guess, target, solver, model_options, corrector)
%   builds the solver options from the shared solver settings and the
%   model options, makes exactly ONE call
%       [sol, info] = corrector(theta, guess, target, opts)
%   and returns
%       R.sol          the returned iterate: sol, or the solver's retained
%                      rejected iterate (mesh blow-up / non-finite), or []
%       R.opts         the options actually passed (reltol, abstol, nmax, ...)
%       R.warning_id   bvp4c warning identifier ('' if none; NaN-free text)
%       R.warning_msg
%       R.corrector_class, R.corrector_reason, R.wall_s
%       R.corrector_nodal_min   the solver's own nodal minimum (mesh nodes)
%       R.mesh_trace   the mesh trace of the call
%                      (sa_mesh_trace 'take': initial mesh, every proposal,
%                      max_proposed, over_nmax); NaN fields when nothing was
%                      recorded (a stub corrector, or a throw before N0)
%       R.max_proposed_mesh   P = R.mesh_trace.max_proposed
%   `corrector` defaults to the copied sap_p28_solve. Tests inject a stub with
%   the same signature; this is the only place a solver is called. The
%   copy's bvp4c resolves to models/original_polarity/private/bvp4c.m, the
%   forwarder to the traced copy; the
%   recorder is reset before the solver call and taken after it on every path
%   (a stale token is SA:TRACE:E_STALE, raised here, outside the solver's
%   own try/catch).

if nargin < 6 || isempty(corrector), corrector = @sap_p28_solve; end
opts = struct('reltol', solver.reltol, 'abstol', solver.abstol, ...
    'nmax', solver.nmax, 'mesh_reject', solver.mesh_reject, ...
    'zero_forcing_root', false, ...
    'seam_tol', mo.seam_tol, 'g6a_tol', mo.g6a_tol, 'g6b_tol', mo.g6b_tol, ...
    'g6b_refine_tol', mo.g6b_refine_tol, 'antisym_abs', mo.antisym_abs, ...
    'antisym_rel', mo.antisym_rel, 'vel_abs', mo.vel_abs, 'vel_rel', mo.vel_rel);

tok = sa_mesh_trace('reset');
t0 = tic;
try
    [sol, info] = corrector(theta, guess, target, opts);
catch ME
    sa_mesh_trace('reset');            % discard: nothing leaks into the next call
    rethrow(ME);
end
R = struct();
R.wall_s = toc(t0);
tr = sa_mesh_trace('take', tok, opts.nmax);
R.opts = opts;
if isempty(sol) && isstruct(info) && isfield(info, 'rejected_sol') && ~isempty(info.rejected_sol)
    sol = info.rejected_sol;
end
R.sol = sol;
R.warning_id = ''; R.warning_msg = '';
if isfield(info, 'convergence') && isstruct(info.convergence)
    R.warning_id  = char(info.convergence.warning_id);
    R.warning_msg = char(info.convergence.warning_msg);
elseif ~isempty(sol)
    % no classifier ran (rejected before classification): the warning state is
    % unknown, so the convergence warning check must fail closed
    R.warning_id = '<unknown>';
end
R.corrector_class  = getf(info, 'class', 'UNKNOWN');
R.corrector_reason = getf(info, 'reason', '');
R.corrector_nodal_min = getf(info, 'minimum', NaN);
R.mesh_trace = tr;
R.max_proposed_mesh = tr.max_proposed;
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
