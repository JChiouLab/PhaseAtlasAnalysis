function H = gnx_sweep_hss(P)
%GNX_SWEEP_HSS  Uniform no-peak (HSS) roots on the (dn_R, dn_B) grid (zero bvp4c calls; M. xanthus model).
%   dn_R, dn_B are density differences (n_left - n_right)/(L/2); the paper's fraction is dn/(2*rho).
%
%   H = GNX_SWEEP_HSS(P)     P = gnx_model('parameters', ...)
%   Grid dn = linspace(-h, h, res); rows i index dn_B, columns j index dn_R;
%   half-domain densities rho + dn/2. Per cell:
%   gnx_hss_roots, all roots kept; the representative (1 x 4 [Rm Rc Bm Bc])
%   exists only when the validity record passes (one root, both routes agree,
%   residual, positivity); otherwise none (no root is reported when a check fails).

res = P.res;
dn = linspace(-P.dn_half_range, P.dn_half_range, res);
H = struct();
H.dn = dn;
H.roots = cell(res, res);
H.n_roots = zeros(res, res);
H.SOL = cell(res, res);
H.existence_mask = false(res, res);
H.validity = cell(res, res);
H.flag = strings(res, res);
for i = 1:res
    for j = 1:res
        densR = P.rho_R + dn(j) / 2;
        densB = P.rho_B + dn(i) / 2;
        [r4, info] = gnx_hss_roots(densR, densB, P);
        H.roots{i, j} = r4;
        H.n_roots(i, j) = size(r4, 1);
        H.flag(i, j) = string(info.flag);
        root = [];
        if size(r4, 1) == 1, root = r4(1, :); end
        rec = gnx_validity('hss', root, info, struct('par', P));
        rec.diag.info = info;
        rec.diag.dens = [densR densB];
        H.validity{i, j} = rec;
        if rec.valid
            H.existence_mask(i, j) = true;
            H.SOL{i, j} = root;
        end
    end
end
end
