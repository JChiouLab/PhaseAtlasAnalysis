function out_path = gnx_write_library(out_dir, type, data, P, meta)
%GNX_WRITE_LIBRARY  Write one labelled Guzzo_noXA seed library (packet §2c, §4j).
%
%   out_path = GNX_WRITE_LIBRARY(out_dir, type, data, P, meta)
%   type: 'HSS' (data from gnx_sweep_hss), 'HPSS' (gnx_sweep_hpss on Lh) or
%   'MPSS' (gnx_build_mpss). meta: run_id, model_version, and optional
%   provenance fields (copied into `provenance`). Writes
%   branch_<TYPE>_<run_id>.mat and copies it to branch_<TYPE>_latest.mat.
%   Labels: model_id, model_version, schema, fields, state_layout,
%   orientation, domain, storage_kind; axes dn_R, dn_B (1 x res) beside
%   R_DIFF, B_DIFF (res x res, = dn*Lh) and Lh. No PAK_DIFF, GEF_DIFF, M_BG,
%   M_P or a..f are written; `par` is the normalized noxa struct, so the
%   original reader stops in check_par (error SAP:LIBRARY).
%   Refuses (error GNX:LIBRARY) before writing if: a stored BVP struct has
%   accumulators that do not start at 0 or end off nR/nB beyond acc_end_rel
%   (RC-7), existence and SOL / flux disagree, or a stored cell's validity
%   record does not pass on recomputation.

D = gnx_model('descriptor');
res = P.res;
Lh = P.L / 2;
dn = linspace(-P.dn_half_range, P.dn_half_range, res);
[R_DIFF, B_DIFF] = meshgrid(dn * Lh, dn * Lh);
L = struct();
L.model_id = D.id;
L.model_version = meta.model_version;
L.schema = D.schema;
L.fields = D.fields;
L.state_layout = D.state_layout;
L.branch_type = type;
L.dn_R = dn;
L.dn_B = dn;
L.R_DIFF = R_DIFF;
L.B_DIFF = B_DIFF;
L.Lh = Lh;
L.res = res;
L.par = P;
L.existence_mask = logical(data.existence_mask);
L.SOL = data.SOL;
L.validity = data.validity;
L.flux_values = nan(res, res, 4);
L.flux_fields = D.fields;
L.flux_source_location = strings(res, res);
L.mean_cytosol = nan(res, res, 2);
L.mean_cytosol_fields = {'Rc', 'Bc'};
L.is_stable = false(res, res);
L.is_stable_computed = false(res, res);
L.stability_source = 'not_computed';
L.run_id = meta.run_id;
L.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
L.provenance = meta;
switch type
    case 'HSS'
        L.orientation = 'uniform';
        L.domain = [0 Lh];
        L.storage_kind = 'uniform_state';
        L.construction = 'algebraic_two_route';
        L.roots = data.roots;
        L.n_roots = data.n_roots;
        L.hss_flag = data.flag;
        L.n_calls = 0; L.n_nmax = 0;
        for i = 1:res
            for j = 1:res
                if ~L.existence_mask(i, j), continue; end
                u = L.SOL{i, j};
                L.flux_values(i, j, :) = u;
                L.flux_source_location(i, j) = "uniform";
                L.mean_cytosol(i, j, :) = [u(2) u(4)];
            end
        end
    case 'HPSS'
        L.orientation = 'peak_at_0';
        L.domain = [0 Lh];
        L.storage_kind = 'bvp4c_struct';
        L.construction = 'direct_BVP_myxo_search';
        L.last_reason = data.last_reason;
        L.n_attempts = data.n_attempts;
        L.phase_log = data.phase_log;
        L.n_calls = data.n_calls; L.n_nmax = data.n_nmax;
        L.cap_reached = data.cap_reached;
        for i = 1:res
            for j = 1:res
                if ~L.existence_mask(i, j), continue; end
                s = L.SOL{i, j};
                check_struct(s, (P.rho_R + dn(j)/2) * Lh, (P.rho_B + dn(i)/2) * Lh, D, i, j);
                L.flux_values(i, j, :) = s.y([1 3 5 7], end);
                L.flux_source_location(i, j) = "seam_x_Lh";
                L.mean_cytosol(i, j, :) = L.validity{i, j}.diag.cyt_mean_121;
            end
        end
    case 'MPSS'
        L.orientation = 'interior_peak';
        L.domain = [0 Lh; 0 P.L/4];
        L.storage_kind = 'quarter_bvp4c+stitched_xy';
        L.construction = 'reflected_quarter_half_peak';
        L.SOL_quarter = data.SOL_quarter;
        L.fail_reason = data.fail_reason;
        L.quarter_ok_grid = data.quarter.ok_grid;
        L.quarter_last_reason = data.quarter.last_reason;
        L.quarter_validity = data.quarter.validity;
        L.phase_log = data.quarter.phase_log;
        L.n_calls = data.n_calls; L.n_nmax = data.n_nmax;
        L.cap_reached = data.quarter.cap_reached;
        ell_q = P.L / 4;
        for i = 1:res
            for j = 1:res
                if ~L.existence_mask(i, j), continue; end
                sq = L.SOL_quarter{i, j};
                check_struct(sq, (P.rho_R + dn(j)/2) * ell_q, (P.rho_B + dn(i)/2) * ell_q, D, i, j);
                L.flux_values(i, j, :) = L.SOL{i, j}.y(:, end);
                L.flux_source_location(i, j) = "outer_trough_x_Lh";
                L.mean_cytosol(i, j, :) = data.quarter.validity{i, j}.diag.cyt_mean_121;
            end
        end
    otherwise
        error('GNX:LIBRARY', 'Unknown library type %s.', type);
end
% consistency: existence <=> non-empty SOL <=> finite flux; valid records
ex = L.existence_mask;
nonempty = ~cellfun(@isempty, L.SOL);
fv = reshape(L.flux_values, res * res, 4);
finite_flux = reshape(all(isfinite(fv), 2), res, res);
if ~isequal(ex, nonempty) || ~isequal(ex, finite_flux)
    error('GNX:LIBRARY', '%s: existence, SOL and flux disagree.', type);
end
for k = find(ex)'
    if ~gnx_validity('verdict', L.validity{k})
        error('GNX:LIBRARY', '%s: stored cell %d has a failing validity record.', type, k);
    end
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
out_path = fullfile(out_dir, sprintf('branch_%s_%s.mat', type, meta.run_id));
save(out_path, '-struct', 'L');
copyfile(out_path, fullfile(out_dir, sprintf('branch_%s_latest.mat', type)));
end

function check_struct(s, nR, nB, D, i, j)
if ~(isstruct(s) && isfield(s, 'y') && size(s.y, 1) == 12)
    error('GNX:LIBRARY', 'cell (%d,%d): stored struct is not a 12-row BVP struct.', i, j);
end
a0 = s.y([11 12], 1);
aend = s.y([11 12], end);
rel = abs(aend - [nR; nB]) ./ [nR; nB];
if any(a0 ~= 0) || any(aend == 0) || any(rel > D.tol.acc_end_rel)
    error('GNX:LIBRARY', ['cell (%d,%d): accumulator rows start at [%g %g] and end at ' ...
        '[%g %g] (targets %g %g); refusing a zeroed or off-mass struct (RC-7).'], ...
        i, j, a0(1), a0(2), aend(1), aend(2), nR, nB);
end
end
