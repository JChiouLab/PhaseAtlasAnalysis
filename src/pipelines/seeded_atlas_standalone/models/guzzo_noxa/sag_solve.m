function R = sag_solve(theta, guess, target, solver, mo, corrector)
%SAG_SOLVE  One bvp4c call on the coupled M. xanthus BVP, in the common form.
%
%   R = SAG_SOLVE(theta, guess, target, solver, model_options, corrector)
%   makes exactly ONE solver call (the coupled-BVP solve) and returns the shared R struct
%       R.sol          the returned iterate (whatever bvp4c returned; [] when
%                      it threw). A return is NOT convergence: the shared
%                      validity checks decide.
%       R.opts         reltol, abstol, nmax, mesh_reject (+ the validity tolerances)
%       R.warning_id, R.warning_msg   lastwarn after the call (reset before;
%                      warnings are never suppressed)
%       R.corrector_class   'RETURNED' | 'NO_RETURN'
%       R.corrector_reason  '' | 'bvp4c_error: <id> <msg>' |
%                      'singularity_guard: <n> hit(s); <id>'
%       R.wall_s, R.corrector_nodal_min (min of the 8 concentration rows at
%                      the returned mesh nodes; NaN without a return)
%   `corrector` (tests) is a stub with the signature
%       [sol, info] = corrector(theta, guess, target, opts)
%   Default: bvp4c with bvpset(RelTol, AbsTol, NMax), no other option.
%       R.mesh_trace, R.max_proposed_mesh   the
%                      call's mesh trace (sa_mesh_trace 'take') and its P; NaN
%                      when nothing was recorded (a stub, or a throw before N0)
%   n = SAG_SOLVE('count') / SAG_SOLVE('reset'): real bvp4c invocations in
%   this process.
%   bvp4c here resolves to models/guzzo_noxa/private/bvp4c.m, the forwarder
%   to the traced copy; the recorder is reset before the call and taken after
%   it on every path (bvp4c errors are caught below; a stale token is
%   SA:TRACE:E_STALE, raised outside that catch).

persistent n_calls
if isempty(n_calls), n_calls = 0; end
if ischar(theta)
    switch theta
        case 'count', R = n_calls;
        case 'reset', n_calls = 0; R = 0;
        otherwise, error('SAG:SOLVE', 'Unknown action %s.', theta);
    end
    return
end
if nargin < 6, corrector = []; end
opts = struct('reltol', solver.reltol, 'abstol', solver.abstol, ...
    'nmax', solver.nmax, 'mesh_reject', solver.mesh_reject, ...
    'zero_forcing_root', false, 'seam_tol', mo.seam_tol, 'acc_tol', mo.acc_tol, ...
    'mass_dense_tol', mo.mass_dense_tol, 'mass_refine_tol', mo.mass_refine_tol, ...
    'antisym_abs', mo.antisym_abs, 'antisym_rel', mo.antisym_rel, ...
    'vel_abs', mo.vel_abs, 'vel_rel', mo.vel_rel);
R = struct('sol', [], 'opts', opts, 'warning_id', '', 'warning_msg', '', ...
    'corrector_class', 'NO_RETURN', 'corrector_reason', '', 'wall_s', NaN, ...
    'corrector_nodal_min', NaN, 'mesh_trace', [], 'max_proposed_mesh', NaN);

tok = sa_mesh_trace('reset');
if ~isempty(corrector)
    t0 = tic;
    try
        [sol, info] = corrector(theta, guess, target, opts);
    catch ME
        sa_mesh_trace('reset');        % discard: nothing leaks into the next call
        rethrow(ME);
    end
    R.wall_s = toc(t0);
    R.sol = sol;
    R.warning_id = getf(info, 'warning_id', '');
    R.warning_msg = getf(info, 'warning_msg', '');
    R.corrector_reason = getf(info, 'reason', '');
else
    prevw = warning;
    warning('on', 'MATLAB:bvp4c:RelTolNotMet'); warning('on', 'MATLAB:bvp4c:IllCondJac');
    cl = onCleanup(@() warning(prevw)); %#ok<NASGU>
    sag_coupled24('guard_reset');
    lastwarn('', '');
    n_calls = n_calls + 1;
    t0 = tic;
    sol = [];
    err = '';
    try
        sol = bvp4c(@(x, y) sag_coupled24('ode', x, y, theta), ...
                    @(ya, yb) sag_coupled24('bc', ya, yb, theta, target), guess, ...
                    bvpset('RelTol', opts.reltol, 'AbsTol', opts.abstol, 'NMax', opts.nmax));
    catch ME
        err = sprintf('%s %s', ME.identifier, ME.message);
    end
    R.wall_s = toc(t0);
    [wmsg, wid] = lastwarn;                % lastwarn returns [message, identifier]
    R.warning_id = wid; R.warning_msg = wmsg;
    hits = sag_coupled24('guard_count');
    R.sol = sol;
    if hits > 0
        R.corrector_reason = sprintf('singularity_guard: %d hit(s); %s', hits, err);
    elseif ~isempty(err)
        R.corrector_reason = ['bvp4c_error: ' err];
    end
end
R.mesh_trace = sa_mesh_trace('take', tok, opts.nmax);
R.max_proposed_mesh = R.mesh_trace.max_proposed;
if ~isempty(R.sol)
    R.corrector_class = 'RETURNED';
    try
        R.corrector_nodal_min = min(R.sol.y([1 3 5 7 9 11 13 15], :), [], 'all');
    catch
    end
end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
