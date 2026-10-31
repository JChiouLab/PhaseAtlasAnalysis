function rec = gnx_validity(kind, varargin)
%GNX_VALIDITY  Guzzo_noXA validity record (packet §4; missing check = invalid).
%
%   rec = GNX_VALIDITY('hpss',    sol, srec, ctx)   half-domain BVP struct
%   rec = GNX_VALIDITY('quarter', sol, srec, ctx)   quarter-domain BVP struct
%       ctx: par, ell, nR, nB (segment masses). Gated checks:
%         solver_converged  srec.converged (stats.maxres < RelTol, mesh <= NMax,
%                           no RelTolNotMet/IllCondJac/SingJac warning)
%         finite            x, y, yp finite and real
%         positivity        exact per-interval Hermite minimum of Rm Rc Bm Bc > 0
%         mass              dense (n_dense-point deval) trapezoid mass,
%                           relative error <= tol.mass_rel per pool
%         accumulator_end   |mR(ell) - nR| <= acc_end_rel*nR, same for mB
%         orientation       Bm(0) >= Bm(ell) (peak at x = 0)
%         morphology        gnx_morphology on deval at morph_nodes uniform
%                           nodes: class HPSS (amplitude = relative contrast);
%                           stores far_trough_gap_rel, autocorr_lag1, mesh
%                           (packet v1.2 condition 1)
%       Recorded, not gated: sigma (|sigmaR|, |sigmaB|), bvp4c threshold
%       AbsTol/RelTol, per-row max|f| on the mesh, residual diagnostic
%       |y' - F| / max(|F|, threshold) sup and RMS per row on the dense grid,
%       the dense-grid positivity cross-check.
%   rec = GNX_VALIDITY('mpss', qrec, solq, ctx)
%       Stitched MPSS from a quarter struct whose record is qrec. ctx: par,
%       ell_q, nR_half, nB_half. Gated: quarter_valid, seam_slope
%       (|Bm'(0)| <= mpss_slope), symmetry (<= mpss_sym), mass_doubling
%       (dense, relative <= mass_rel), peak_fraction in (lo, hi), contrast
%       ((max-min)/max of Bm >= contrast). Also returns rec.x, rec.y (stitched).
%   rec = GNX_VALIDITY('hss', root, info, ctx)
%       Uniform root [Rm Rc Bm Bc]; info from gnx_hss_roots; ctx: par.
%       Gated: finite, positivity (> 0), residual (|sR| <= hss_resid_rel*sum|terms|,
%       same for sB), route_agreement, single_root.
%   ok = GNX_VALIDITY('verdict', rec)
%       Recomputes the verdict from the record: every required check of
%       rec.kind present with pass == true (fail closed).
%   names = GNX_VALIDITY('required', kind)

switch kind
    case 'verdict',  rec = verdict(varargin{1}); return;
    case 'required', rec = required(varargin{1}); return;
    case {'hpss', 'quarter'}, rec = bvp_record(kind, varargin{:});
    case 'mpss', rec = mpss_record(varargin{:});
    case 'hss',  rec = hss_record(varargin{:});
    otherwise, error('GNX:VALIDITY', 'Unknown kind %s.', kind);
end
rec.valid = verdict(rec);
end

% -------------------------------------------------------------------------
function names = required(kind)
switch kind
    case {'hpss', 'quarter'}
        names = {'solver_converged', 'finite', 'positivity', 'mass', ...
                 'accumulator_end', 'orientation', 'morphology'};
    case 'mpss'
        names = {'quarter_valid', 'seam_slope', 'symmetry', 'mass_doubling', ...
                 'peak_fraction', 'contrast'};
    case 'hss'
        names = {'finite', 'positivity', 'residual', 'route_agreement', 'single_root'};
    otherwise
        error('GNX:VALIDITY', 'Unknown kind %s.', kind);
end
end

function ok = verdict(rec)
ok = false;
if ~isstruct(rec) || ~isfield(rec, 'kind') || ~isfield(rec, 'checks'), return; end
names = required(rec.kind);
for k = 1:numel(names)
    if ~isfield(rec.checks, names{k}), return; end
    c = rec.checks.(names{k});
    if ~isstruct(c) || ~isfield(c, 'pass') || ~isequal(c.pass, true), return; end
end
ok = true;
end

function c = chk(pass, value, tol, varargin)
c = struct('pass', logical(pass), 'value', value, 'tol', tol);
for k = 1:2:numel(varargin), c.(varargin{k}) = varargin{k+1}; end
end

% -------------------------------------------------------------------------
function rec = bvp_record(kind, sol, srec, ctx)
D = gnx_model('descriptor');
nm = D.numerics; tol = D.tol;
rec = struct('kind', kind, 'checks', struct(), 'diag', struct(), 'valid', false);
rec.ell = ctx.ell; rec.nR = ctx.nR; rec.nB = ctx.nB;
conv = isstruct(srec) && isfield(srec, 'converged') && srec.converged;
rec.checks.solver_converged = chk(conv, getf(srec, 'maxres', NaN), getf(srec, 'RelTol', NaN), ...
    'nmeshpoints', getf(srec, 'nmeshpoints', NaN), 'NMax', getf(srec, 'NMax', NaN), ...
    'warn_id', getf(srec, 'warn_id', ''), 'reason', getf(srec, 'reason', ''), ...
    'threshold', getf(srec, 'threshold', NaN));
if ~(isstruct(sol) && isfield(sol, 'x') && isfield(sol, 'y') && isfield(sol, 'yp') ...
        && size(sol.y, 1) == 12 && size(sol.yp, 1) == 12 && ...
        numel(sol.x) == size(sol.y, 2) && numel(sol.x) == size(sol.yp, 2) && numel(sol.x) >= 2)
    rec.checks.finite = chk(false, NaN, NaN, 'reason', 'malformed_solution');
    return;
end
fin = all(isfinite(sol.x(:))) && all(isfinite(sol.y(:))) && all(isfinite(sol.yp(:))) && ...
      isreal(sol.x) && isreal(sol.y) && isreal(sol.yp);
rec.checks.finite = chk(fin, double(fin), NaN);
if ~fin, return; end
x = sol.x; ell = x(end);
fields = D.fields; rows = D.conc_rows;
% positivity: exact Hermite minimum per field
mins = nan(1, 4); xat = nan(1, 4);
for f = 1:4
    [mins(f), xat(f)] = gnx_hermite_min(x, sol.y(rows(f), :), sol.yp(rows(f), :));
end
[mn, fi] = min(mins);
% dense grid (positivity cross-check, mass, residual diagnostic)
xd = linspace(0, ell, nm.n_dense);
[Yd, YPd] = deval(sol, xd);
dense_min = min(Yd(rows, :), [], 2).';
rec.checks.positivity = chk(mn > 0, mn, 0, 'field', fields{fi}, 'x', xat(fi), ...
    'per_field_min', mins, 'dense_min', dense_min);
mR = trapz(xd, Yd(1, :) + Yd(3, :));
mB = trapz(xd, Yd(5, :) + Yd(7, :));
eR = abs(mR - ctx.nR) / abs(ctx.nR); eB = abs(mB - ctx.nB) / abs(ctx.nB);
rec.checks.mass = chk(max(eR, eB) <= tol.mass_rel, [eR eB], tol.mass_rel, ...
    'dense_mass', [mR mB], 'target', [ctx.nR ctx.nB]);
aR = abs(sol.y(11, end) - ctx.nR) / abs(ctx.nR);
aB = abs(sol.y(12, end) - ctx.nB) / abs(ctx.nB);
a0 = [sol.y(11, 1) sol.y(12, 1)];
rec.checks.accumulator_end = chk(max(aR, aB) <= tol.acc_end_rel && all(a0 == 0), ...
    [aR aB], tol.acc_end_rel, 'start', a0);
bm0 = sol.y(5, 1); bml = sol.y(5, end);
rec.checks.orientation = chk(bm0 >= bml, [bm0 bml], NaN);
xg = linspace(0, ell, nm.morph_nodes);
yg = deval(sol, xg);
[acc, info] = gnx_morphology(xg, yg(rows, :));
% v1.2 condition 1 (Z7 ruling): the record names the evaluation mesh and
% carries far_trough_gap/amp and autocorr_lag1 explicitly (NaN when the
% classifier stopped before computing them, e.g. at the amplitude gate).
rec.checks.morphology = chk(acc, string(info.class), "HPSS", 'reason', string(info.reason), ...
    'info', info, 'far_trough_gap_rel', info.far_trough_gap / info.amp, ...
    'autocorr_lag1', info.autocorr_lag1, ...
    'mesh', sprintf('deval_uniform_%d', nm.morph_nodes), 'mesh_nodes', nm.morph_nodes);
% diagnostics (recorded, not gated)
rec.diag.sigma = abs(sol.y([9 10], 1)).';
thr = getf(srec, 'threshold', NaN);
rec.diag.threshold = thr;
rec.diag.row_max_abs_f = max(abs(sol.yp), [], 2).';
Fd = gnx_bvp_ode(xd, Yd, ctx.par);
den = max(abs(Fd), thr);
r = abs(YPd - Fd) ./ den;
rec.diag.resid_sup = max(r, [], 2).';
rec.diag.resid_rms = sqrt(mean(r.^2, 2)).';
rec.diag.maxres = getf(srec, 'maxres', NaN);
rec.diag.n_mesh = numel(x);
rec.diag.cyt_mean_121 = [mean(yg(3, :)) mean(yg(7, :))];   % myxo :253 convention
rec.diag.seam = sol.y(rows, end).';
end

% -------------------------------------------------------------------------
function rec = mpss_record(qrec, solq, ctx)
D = gnx_model('descriptor'); tol = D.tol; nm = D.numerics;
rec = struct('kind', 'mpss', 'checks', struct(), 'diag', struct(), 'valid', false, ...
    'x', [], 'y', []);
qok = isstruct(qrec) && isfield(qrec, 'kind') && strcmp(qrec.kind, 'quarter') && verdict(qrec);
rec.checks.quarter_valid = chk(qok, double(qok), NaN);
if ~qok, return; end
xq = solq.x(:).'; ellq = xq(end);
yq = solq.y([1 3 5 7], :);
% peak-to-peak reflection (build_mpss_branch.m:625-650 convention): the
% quarter's peak (x = 0) goes to the seam x = ell_q of the stitched half.
x_left = ellq - fliplr(xq);  x_right = ellq + xq;
xs = [x_left(1:end-1), x_right];
ys = [fliplr(yq(:, 2:end)), yq];
rec.x = xs; rec.y = ys;
seam = abs(solq.y(6, 1));
rec.checks.seam_slope = chk(seam <= tol.mpss_slope, seam, tol.mpss_slope);
sym = max(abs(ys - fliplr(ys)), [], 'all');
rec.checks.symmetry = chk(sym <= tol.mpss_sym, sym, tol.mpss_sym);
% dense mass of the stitched profile from the quarter's Hermite interpolant
xd = linspace(0, 2 * ellq, 2 * nm.n_dense - 1);
xm = abs(xd - ellq);
xm(xm > ellq) = ellq;
Yd = deval(solq, xm);
mR = trapz(xd, Yd(1, :) + Yd(3, :));
mB = trapz(xd, Yd(5, :) + Yd(7, :));
eR = abs(mR - ctx.nR_half) / ctx.nR_half; eB = abs(mB - ctx.nB_half) / ctx.nB_half;
rec.checks.mass_doubling = chk(max(eR, eB) <= tol.mass_rel, [eR eB], tol.mass_rel, ...
    'dense_mass', [mR mB]);
[bmax, ipk] = max(ys(3, :));
pf = xs(ipk) / xs(end);
rec.checks.peak_fraction = chk(pf > tol.mpss_peak_lo && pf < tol.mpss_peak_hi, pf, ...
    [tol.mpss_peak_lo tol.mpss_peak_hi]);
ctr = (bmax - min(ys(3, :))) / max(bmax, 1e-12);
rec.checks.contrast = chk(ctr >= tol.contrast, ctr, tol.contrast);
rec.diag.outer_trough = ys(:, end).';
end

% -------------------------------------------------------------------------
function rec = hss_record(root, info, ctx)
D = gnx_model('descriptor'); tol = D.tol;
rec = struct('kind', 'hss', 'checks', struct(), 'diag', struct(), 'valid', false);
n_roots = getf(info, 'n_roots', NaN);
agree = isfield(info, 'flag') && strcmp(info.flag, 'ok');
rec.checks.route_agreement = chk(agree, getf(info, 'route_max_rel', NaN), tol.hss_route_rel, ...
    'flag', getf(info, 'flag', ''));
rec.checks.single_root = chk(n_roots == 1, n_roots, 1);
if isempty(root) || numel(root) ~= 4
    rec.checks.finite = chk(false, NaN, NaN, 'reason', 'no_representative');
    return;
end
u = root(:);
fin = all(isfinite(u)) && isreal(u);
rec.checks.finite = chk(fin, double(fin), NaN);
if ~fin, return; end
rec.checks.positivity = chk(all(u > 0), min(u), 0);
[sR, sB, aR, aB] = gnx_reaction(u, ctx.par);
rR = abs(sR) / aR; rB = abs(sB) / aB;
rec.checks.residual = chk(abs(sR) <= tol.hss_resid_rel * aR && abs(sB) <= tol.hss_resid_rel * aB, ...
    [rR rB], tol.hss_resid_rel, 'sR', sR, 'sB', sB);
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f), v = s.(f); else, v = d; end
end
