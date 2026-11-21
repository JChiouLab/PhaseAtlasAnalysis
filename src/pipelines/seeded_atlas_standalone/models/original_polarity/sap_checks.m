function C = sap_checks(sol, theta, target, mo)
%SAP_CHECKS  Validity checks of the S. pombe coupled solution.
%
%   C = SAP_CHECKS(sol, theta, target, model_options) recomputes the
%   diagnostics of the coupled solve from the returned solution (never trusting
%   the solver's own flags) and returns the common check records
%       name, value, tolerance, comparison, gating, owner (the defining module)
%   Gating: middle concentration, middle flux, accumulator mass,
%   dense mass and its 20001/40001 refinement, source antisymmetry
%   and flux-vs-source velocity consistency. Recorded, not gating: s_max
%   and the four sources (forced constrained states).
%   The sap_forced28 diagnostics evaluate the dense masses with deval, so sol
%   must carry its bvp4c interpolant data (x, y, yp).

D = sap_forced28('diagnostics', sol, theta, target);
antisym_gate = max(mo.antisym_abs, mo.antisym_rel * D.s_max);
vel_gate     = max(mo.vel_abs, mo.vel_rel * D.vel_max);

C = [ ...
  rec('seam_concentration',    max(abs(D.seam_conc_mismatch)), mo.seam_tol, '<', true)
  rec('seam_flux',             D.max_seam_flux,                 mo.seam_tol, '<', true)
  rec('mass_accumulator',      D.G6a_residual,                  mo.g6a_tol,  '<', true)
  rec('mass_dense',            D.G6b_residual,                  mo.g6b_tol,  '<', true)
  rec('mass_dense_refinement', D.G6b_refine,                    mo.g6b_refine_tol, '<', true)
  rec('shift_antisymmetry',    D.antisym_max,                   antisym_gate, '<', true)
  rec('velocity_consistency',  D.vel_discrepancy_max,           vel_gate,     '<', true)
  rec('s_max',                 D.s_max,                         NaN, 'record', false)
  rec('s_GL',                  D.s_GL,                          NaN, 'record', false)
  rec('s_PL',                  D.s_PL,                          NaN, 'record', false)
  rec('s_GR',                  D.s_GR,                          NaN, 'record', false)
  rec('s_PR',                  D.s_PR,                          NaN, 'record', false)
  ].';
end

function r = rec(name, value, tol, cmp, gating)
r = struct('name', name, 'value', double(value), 'tolerance', double(tol), ...
    'comparison', cmp, 'gating', gating, 'owner', 'module');
end
