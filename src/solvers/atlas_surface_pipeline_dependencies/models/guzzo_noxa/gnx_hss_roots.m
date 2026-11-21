function [roots4, info] = gnx_hss_roots(densR, densB, par)
%GNX_HSS_ROOTS  Homogeneous steady states of the M. xanthus (Guzzo_noXA) model by two routes.
%
%   [roots4, info] = GNX_HSS_ROOTS(densR, densB, par)
%   densR, densB: mean densities of the segment. Uniform states satisfy
%   sR = sB = 0 with Rc = densR - Rm, Bc = densB - Bm; RomR's balance gives
%   Rm = g*densR/(g + dR), g = kR + kRB*Bm, and one scalar equation
%   F(Bm) = 0 on [0, densB] remains.
%   Route (i): F*(1 + Bm/K)*(1 + (Bm/KA)^2)*(g + dR) is a degree-6
%     polynomial; its roots with |Im| <= hss_imag*max(1,|Re|) and real part
%     in [0, densB] are polished on F (fzero on a sign-change bracket around
%     the root), never on the polynomial.
%   Route (ii): the reference M. xanthus code's 2001-point scan plus fzero,
%     re-implemented here.
%   Accepted (info.flag = 'ok') when both routes give the same count, the
%   roots agree to hss_route_rel relative, and each root satisfies
%   |sR| <= hss_resid_rel*sum|terms|, same for sB. Otherwise info.flag =
%   'hss_route_disagreement' (existence false at the caller).
%   roots4 rows are [Rm Rc Bm Bc] (field order; the reference M. xanthus code stores [Bm Rm Bc Rc]),
%   sorted by Bm; they are route (i)'s polished roots.

D = gnx_model('descriptor'); tol = D.tol;
F = @(bm) balB(bm, densR, densB, par);
% ---- route (i): polynomial ----
pc = poly_coeffs(densR, densB, par);
r = roots(pc);
isre = abs(imag(r)) <= tol.hss_imag * max(1, abs(real(r)));
rr = sort(real(r(isre)));
rr = rr(rr >= 0 & rr <= densB);
poly_bm = zeros(0, 1);
polish_ok = true;
for k = 1:numel(rr)
    [bm, okp] = polish(F, rr(k), densB);
    polish_ok = polish_ok && okp;
    poly_bm(end+1, 1) = bm; %#ok<AGROW>
end
poly_bm = sort(poly_bm);
% ---- route (ii): scan + fzero ----
n_scan = 2001;
bms = linspace(0, densB, n_scan);
vals = arrayfun(F, bms);
scan_bm = zeros(0, 1);
for k = 1:n_scan-1
    if vals(k) == 0, vals(k) = eps; end
    if sign(vals(k)) ~= sign(vals(k+1))
        scan_bm(end+1, 1) = fzero(F, [bms(k), bms(k+1)]); %#ok<AGROW>
    end
end
scan_bm = sort(scan_bm);
% ---- assemble ----
roots4 = zeros(numel(poly_bm), 4);
resid = zeros(numel(poly_bm), 2);
resid_ok = true;
for k = 1:numel(poly_bm)
    bm = poly_bm(k);
    rm = rmOf(bm, densR, par);
    roots4(k, :) = [rm, densR - rm, bm, densB - bm];
    [sR, sB, aR, aB] = gnx_reaction(roots4(k, :).', par);
    resid(k, :) = [abs(sR) / aR, abs(sB) / aB];
    resid_ok = resid_ok && abs(sR) <= tol.hss_resid_rel * aR && abs(sB) <= tol.hss_resid_rel * aB;
end
info = struct('n_roots', numel(poly_bm), 'n_scan', numel(scan_bm), ...
    'poly_bm', poly_bm.', 'scan_bm', scan_bm.', 'n_poly_real_in_range', numel(rr), ...
    'route_max_rel', NaN, 'resid_rel', resid, 'polish_ok', polish_ok, 'flag', '');
same_count = numel(poly_bm) == numel(scan_bm);
if same_count && ~isempty(poly_bm)
    info.route_max_rel = max(abs(poly_bm - scan_bm) ./ max(abs(scan_bm), realmin));
elseif same_count
    info.route_max_rel = 0;
end
if same_count && info.route_max_rel <= tol.hss_route_rel && resid_ok && polish_ok
    info.flag = 'ok';
else
    info.flag = 'hss_route_disagreement';
end
end

% -------------------------------------------------------------------------
function pc = poly_coeffs(densR, densB, par)
% F * (1 + b/K) * (1 + b^2/KA^2) * (g + dR), descending powers of b:
%   (kB + kBB b)(densB - b)(1 + b/K)(1 + b^2/KA^2)(g + dR)
%   - dB b (1 + b^2/KA^2)(g + dR)
%   - dBA densR g b^2 (1 + b/K)
A  = [par.kBB, par.kB];
Bp = [-1, densB];
Kp = [1 / par.K, 1];
Qp = [1 / par.KA^2, 0, 1];
G  = [par.kRB, par.kR];
GD = [par.kRB, par.kR + par.dR];
t1 = conv(conv(conv(conv(A, Bp), Kp), Qp), GD);
t2 = par.dB * conv(conv([1 0], Qp), GD);
t3 = par.dBA * densR * conv(conv(G, [1 0 0]), Kp);
n = max([numel(t1), numel(t2), numel(t3)]);
pad = @(v) [zeros(1, n - numel(v)), v];
pc = pad(t1) - pad(t2) - pad(t3);
end

function [bm, ok] = polish(F, r0, densB)
% fzero on a sign-change bracket around r0, inside [0, densB]
ok = false; bm = r0;
for delta = 10.^(-12:-2)
    a = max(0, r0 - delta * max(1, abs(r0)));
    b = min(densB, r0 + delta * max(1, abs(r0)));
    fa = F(a); fb = F(b);
    if fa == 0, bm = a; ok = true; return; end
    if fb == 0, bm = b; ok = true; return; end
    if sign(fa) ~= sign(fb)
        bm = fzero(F, [a, b]); ok = true; return;
    end
end
end

function rm = rmOf(bm, densR, par)
g = par.kR + par.kRB * bm;
rm = g * densR / (g + par.dR);
end

function v = balB(bm, densR, densB, par)
rm = rmOf(bm, densR, par);
aeff = rm / (1 + (bm / par.KA)^2);
v = (par.kB + par.kBB * bm) * (densB - bm) ...
    - par.dB * bm / (1 + bm / par.K) ...
    - par.dBA * aeff * bm^2;
end
