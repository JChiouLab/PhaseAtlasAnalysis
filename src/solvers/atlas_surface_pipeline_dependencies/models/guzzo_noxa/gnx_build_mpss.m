function M = gnx_build_mpss(P, opts)
%GNX_BUILD_MPSS  MPSS = reflected quarter-domain half-peak (packet §4i).
%
%   M = GNX_BUILD_MPSS(P)
%   M = GNX_BUILD_MPSS(P, opts)     opts passed to gnx_sweep_hpss
%   Quarter sweep on ell_q = L/4 at quarter masses (half masses / 2, same
%   densities, so the grid cell is unchanged), then for every accepted quarter
%   a peak-to-peak reflection (build_mpss_branch.m:625-650 convention) onto
%   [0, L/2] and the MPSS record of gnx_validity('mpss'): quarter valid
%   (positivity, orientation on the quarter, ...), seam slope, symmetry,
%   dense mass doubling, peak fraction in (0.4, 0.6), relative contrast.
%   An all-false MPSS result is a valid output.

if nargin < 2 || isempty(opts), opts = struct(); end
ell_q = P.L / 4;
Lh = P.L / 2;
Q = gnx_sweep_hpss(P, ell_q, 'quarter', opts);
res = P.res;
M = struct();
M.quarter = Q;
M.ell_q = ell_q;
M.existence_mask = false(res, res);
M.SOL = cell(res, res);           % stitched {x, y} (4 rows Rm Rc Bm Bc)
M.SOL_quarter = cell(res, res);   % quarter bvp4c struct
M.validity = cell(res, res);
M.fail_reason = strings(res, res);
for i = 1:res
    for j = 1:res
        if ~Q.ok_grid(i, j)
            M.fail_reason(i, j) = "quarter_not_found";
            continue;
        end
        ctx = struct('par', P, 'ell_q', ell_q, ...
            'nR_half', (P.rho_R + Q.dn(j) / 2) * Lh, ...
            'nB_half', (P.rho_B + Q.dn(i) / 2) * Lh);
        rec = gnx_validity('mpss', Q.validity{i, j}, Q.SOLS{i, j}, ctx);
        M.validity{i, j} = rec;
        if rec.valid
            M.existence_mask(i, j) = true;
            M.SOL{i, j} = struct('x', rec.x, 'y', rec.y);
            M.SOL_quarter{i, j} = Q.SOLS{i, j};
        else
            names = gnx_validity('required', 'mpss');
            for k = 1:numel(names)
                if ~isfield(rec.checks, names{k}) || ~rec.checks.(names{k}).pass
                    M.fail_reason(i, j) = string(names{k}); break;
                end
            end
        end
    end
end
M.n_calls = Q.n_calls;
M.n_nmax = Q.n_nmax;
end
