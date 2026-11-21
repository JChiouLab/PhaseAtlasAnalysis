function C = sag_checks(sol, theta, target, mo)
%SAG_CHECKS  Validity checks of the M. xanthus coupled solution.
%
%   C = SAG_CHECKS(sol, theta, target, model_options) recomputes every value
%   from the returned iterate (x, y, yp, the bvp4c interpolant) and returns
%   the common check records  name, value, tolerance, comparison, gating, owner (the module that defines the check).
%
%   Gating (all checks below are defined by this model module)
%     real_valued            count of entries of y, yp, parameters with a
%                            nonzero imaginary part == 0;
%                            evaluated first, on the raw arrays
%     positivity_exact       exact per-interval Hermite minimum of the 8
%                            concentration rows > 0 (sag_hermite_min)
%     seam_concentration     max |Y_L - Y_R| at the seam over 4 fields < seam_tol
%     seam_flux              max |D (Y_L' + Y_R')| at the seam < seam_tol
%     mass_accumulator       max over the 4 half masses |I(Lh) - n| / n <= acc_tol
%     mass_dense             exact Hermite integral of each half-mass density,
%                            max |int - n| / n <= mass_dense_tol (a consistency
%                            check of the interpolant)
%     mass_dense_refinement  the same integral on the bisected mesh (deval
%                            values and slopes at midpoints), max relative
%                            difference <= mass_refine_tol
%     shift_antisymmetry     max(|s_RL + s_RR|, |s_BL + s_BR|)
%                            < max(antisym_abs, antisym_rel * s_max)
%     velocity_consistency   max over the four half masses h |J_h + Lh s_h|
%                            < max(vel_abs, vel_rel * max |J_h|), J_h =
%                            DM Rm_h'(Lh) + DC Rc_h'(Lh) (same for B); the
%                            right-half terms make a solution with s_RR forced
%                            equal to s_RL refusable
%   Recorded, not gating: s_max, the four sources, accumulator start, the
%   positivity_exact row index, the two velocity readings.

L = sag_coupled24('layout');
nimag = count_imag(sol, 'y') + count_imag(sol, 'yp') + count_imag(sol, 'parameters');
C = rec('real_valued', nimag, 0, '==', true);
try
    C = [C, value_checks(sol, theta, target, mo, L)];
catch
    % a check that cannot be evaluated is absent, so the verdict refuses the
    % return (missing required check); real_valued is still recorded
end
end

function C = value_checks(sol, theta, target, mo, L)
x = real(sol.x); y = sol.y; yp = sol.yp;
rows = [L.conc_left L.conc_right];
mins = nan(1, 8); xat = nan(1, 8);
for k = 1:8
    [mins(k), xat(k)] = sag_hermite_min(x, real(y(rows(k), :)), real(yp(rows(k), :)));
end
if any(isnan(mins)), pmin = NaN; kmin = NaN; else, [pmin, kmin] = min(mins); end
D = sag_coupled24('diagnostics', sol, theta, target);
antisym_gate = max(mo.antisym_abs, mo.antisym_rel * D.s_max);
vel_gate = max(mo.vel_abs, mo.vel_rel * D.J_max);
C = [ ...
  rec('positivity_exact',     pmin,                         0,                  '>',  true), ...
  rec('seam_concentration',    max(abs(D.seam_conc_mismatch)), mo.seam_tol,      '<',  true), ...
  rec('seam_flux',             max(abs(D.seam_flux)),        mo.seam_tol,        '<',  true), ...
  rec('mass_accumulator',      max(D.acc_rel),               mo.acc_tol,         '<=', true), ...
  rec('mass_dense',            max(D.mass_dense_rel),        mo.mass_dense_tol,  '<=', true), ...
  rec('mass_dense_refinement', max(D.mass_refine_rel),       mo.mass_refine_tol, '<=', true), ...
  rec('shift_antisymmetry',    D.antisym_max,                antisym_gate,       '<',  true), ...
  rec('velocity_consistency',  D.vel_discrepancy_max,        vel_gate,           '<',  true), ...
  rec('positivity_exact_row',  kmin,                         NaN, 'record', false), ...
  rec('positivity_exact_xi',   xat_of(xat, kmin),            NaN, 'record', false), ...
  rec('accumulator_start',     D.acc_start_max,              NaN, 'record', false), ...
  rec('s_max',                 D.s_max,                      NaN, 'record', false), ...
  rec('sigma_RL',              D.sigma(1),                   NaN, 'record', false), ...
  rec('sigma_BL',              D.sigma(2),                   NaN, 'record', false), ...
  rec('sigma_RR',              D.sigma(3),                   NaN, 'record', false), ...
  rec('sigma_BR',              D.sigma(4),                   NaN, 'record', false), ...
  rec('velocity_R_slope',      D.vel_from_slope(1),          NaN, 'record', false), ...
  rec('velocity_B_slope',      D.vel_from_slope(2),          NaN, 'record', false)];
end

function n = count_imag(sol, f)
if isfield(sol, f) && ~isempty(sol.(f)), v = sol.(f); n = nnz(imag(v(:)) ~= 0); else, n = 0; end
end

function v = xat_of(xat, k)
if isnan(k), v = NaN; else, v = xat(k); end
end

function r = rec(name, value, tol, cmp, gating)
r = struct('name', name, 'value', double(value), 'tolerance', double(tol), ...
    'comparison', cmp, 'gating', gating, 'owner', 'module');
end
