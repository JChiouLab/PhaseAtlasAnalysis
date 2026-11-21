function C = sap_bvp_convergence(sol, opts, wid, wmsg)
%BVP_CONVERGENCE  One shared bvp4c convergence classifier for the coupled-BVP
%   solves.
%
%   CONVERGED  <=>  finite return AND stats.maxres < RelTol AND
%                   NOT(stats.nmeshpoints > NMax).
%   Missing stats are UNRESOLVED_NO_STATS and are NEVER inferred to pass.
%
%   Returns a struct; it changes no equation, boundary condition or gate.
if nargin < 3, wid = ''; end
if nargin < 4, wmsg = ''; end
C = struct('classifier','sap_bvp_convergence/1', ...
    'gate','finite AND stats.maxres < RelTol AND NOT(stats.nmeshpoints > NMax)', ...
    'returned_n_mesh',NaN,'stats_nmeshpoints',NaN,'stats_maxres',NaN,'has_stats',false, ...
    'finite',false,'maxres_lt_reltol',false,'proposed_over_nmax',false, ...
    'warn_reltol_not_met',strcmp(wid,'MATLAB:bvp4c:RelTolNotMet'), ...
    'warn_illcond_jac',strcmp(wid,'MATLAB:bvp4c:IllCondJac'), ...
    'warning_id',wid,'warning_msg',wmsg,'converged',false,'class','NO_RETURN');
if isempty(sol) || ~isstruct(sol) || ~isfield(sol,'x')
    C.class = 'NO_RETURN'; return;
end
C.returned_n_mesh = numel(sol.x);
C.finite = isfield(sol,'y') && all(isfinite(sol.y(:)));
C.has_stats = isfield(sol,'stats') && isstruct(sol.stats) && ...
              isfield(sol.stats,'maxres') && isfield(sol.stats,'nmeshpoints');
rt = getf(opts,'reltol',NaN); nm = getf(opts,'nmax',NaN);
C.reltol = rt; C.nmax = nm;
if C.has_stats
    C.stats_nmeshpoints = sol.stats.nmeshpoints;
    C.stats_maxres = sol.stats.maxres;
    C.maxres_lt_reltol = isfinite(rt) && C.stats_maxres < rt;
    C.proposed_over_nmax = isfinite(nm) && C.stats_nmeshpoints > nm;
end
if ~C.has_stats || ~isfinite(rt)
    C.class = 'UNRESOLVED_NO_STATS';
elseif ~C.finite
    C.class = 'NON_FINITE';
elseif C.warn_reltol_not_met || C.proposed_over_nmax
    C.class = 'UNRESOLVED_RELTOL_NMAX';
elseif ~C.maxres_lt_reltol
    C.class = 'UNRESOLVED_RELTOL';
elseif C.warn_illcond_jac
    C.class = 'UNRESOLVED_ILLCOND';
else
    C.class = 'CONVERGED'; C.converged = true;
end
end
function v = getf(s,f,d), if isstruct(s)&&isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
