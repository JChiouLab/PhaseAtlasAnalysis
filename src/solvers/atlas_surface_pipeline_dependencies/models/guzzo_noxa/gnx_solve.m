function [sol, srec] = gnx_solve(solinit, par, nR, nB, opts)
%GNX_SOLVE  One bvp4c call on the Guzzo_noXA BVP, counted (packet §4g, §8).
%
%   [sol, srec] = GNX_SOLVE(solinit, par, nR, nB)
%   [sol, srec] = GNX_SOLVE(solinit, par, nR, nB, opts)
%       solinit: struct x (1 x n, from 0 to ell), y (12 x n), see gnx_guess.
%       opts (optional): RelTol (1e-8), AbsTol (1e-6), NMax (20000) and
%       bvp_fun (default @bvp4c; tests inject a stub). Every invocation of
%       bvp_fun increments the process counter by exactly one.
%       sol is whatever bvp4c returned ([] on a caught exception); a return
%       is NOT convergence: srec.converged applies the gate
%         stats present, stats.maxres < RelTol, stats.nmeshpoints <= NMax,
%         numel(x) <= NMax, no RelTolNotMet / IllCondJac / SingJac warning
%       (lastwarn reset before the call; warnings are never suppressed).
%       srec.reason: '' | 'bvp4c_error' | 'singularity_guard' |
%       'not_converged' | 'no_stats'.
%   n = GNX_SOLVE('count')          bvp_fun invocations in this process
%   n = GNX_SOLVE('nmax_count')     invocations whose mesh reached NMax
%   GNX_SOLVE('reset')              zero both counters

persistent n_calls n_nmax
if isempty(n_calls), n_calls = 0; n_nmax = 0; end
if ischar(solinit)
    switch solinit
        case 'count', sol = n_calls;
        case 'nmax_count', sol = n_nmax;
        case 'reset', n_calls = 0; n_nmax = 0; sol = 0;
        otherwise, error('GNX:SOLVE', 'Unknown action %s.', solinit);
    end
    return;
end
if nargin < 5 || isempty(opts), opts = struct(); end
RelTol = getf(opts, 'RelTol', 1e-8);
AbsTol = getf(opts, 'AbsTol', 1e-6);
NMax = getf(opts, 'NMax', 20000);
bvp_fun = getf(opts, 'bvp_fun', @bvp4c);

bopts = bvpset('RelTol', RelTol, 'AbsTol', AbsTol, 'NMax', NMax);
srec = struct('RelTol', RelTol, 'AbsTol', AbsTol, 'NMax', NMax, ...
    'threshold', AbsTol / RelTol, 'returned', false, 'error_id', '', ...
    'error_msg', '', 'warn_id', '', 'warn_msg', '', 'guard_hits', 0, ...
    'maxres', NaN, 'nmeshpoints', NaN, 'n_x', NaN, 'nmax_reached', false, ...
    'converged', false, 'reason', '', 'wall_s', NaN, 'call_index', NaN);

gnx_bvp_ode('guard_reset');
lastwarn('', '');
n_calls = n_calls + 1;
srec.call_index = n_calls;
t0 = tic;
sol = [];
try
    sol = bvp_fun(@(x, y) gnx_bvp_ode(x, y, par), ...
                  @(ya, yb) gnx_bvp_bc(ya, yb, nR, nB), solinit, bopts);
    srec.returned = true;
catch ME
    srec.error_id = ME.identifier;
    srec.error_msg = ME.message;
end
srec.wall_s = toc(t0);
[wmsg, wid] = lastwarn;
srec.warn_id = wid; srec.warn_msg = wmsg;
srec.guard_hits = gnx_bvp_ode('guard_count');

if srec.returned && isstruct(sol) && isfield(sol, 'stats') && isstruct(sol.stats) ...
        && isfield(sol.stats, 'maxres') && isfield(sol.stats, 'nmeshpoints')
    srec.maxres = sol.stats.maxres;
    srec.nmeshpoints = sol.stats.nmeshpoints;
    srec.n_x = numel(sol.x);
    srec.nmax_reached = srec.nmeshpoints > NMax || srec.n_x > NMax || ...
        strcmp(wid, 'MATLAB:bvp4c:RelTolNotMet');
    bad_warn = any(strcmp(wid, {'MATLAB:bvp4c:RelTolNotMet', ...
        'MATLAB:bvp4c:IllCondJac', 'MATLAB:bvp4c:SingJac'}));
    srec.converged = srec.maxres < RelTol && srec.nmeshpoints <= NMax && ...
        srec.n_x <= NMax && ~bad_warn;
    if ~srec.converged
        srec.reason = 'not_converged';
        if srec.guard_hits > 0, srec.reason = 'singularity_guard'; end
    end
elseif srec.returned
    srec.reason = 'no_stats';
else
    srec.reason = 'bvp4c_error';
    if srec.guard_hits > 0, srec.reason = 'singularity_guard'; end
end
if srec.nmax_reached, n_nmax = n_nmax + 1; end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
