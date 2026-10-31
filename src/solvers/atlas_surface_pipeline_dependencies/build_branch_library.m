function report = build_branch_library(hpss_path, hss_path, out_dir)
    % BUILD_BRANCH_LIBRARY  Packet 2 adapter: HPSS/HSS -> branch-library schema.
    %
    % Additive, READ-ONLY adapter (docs/ai_handoff/atlas_plan.md, "Half-Domain
    % Branch Library Schema"). Reads the existing canonical HPSS and HSS .mat
    % files and writes NEW branch_HPSS_<id>.mat / branch_HSS_<id>.mat files
    % (plus branch_<TYPE>_latest.mat canonical copies, matching the project's
    % timestamped-then-copyfile convention). NEVER writes to, or modifies, the
    % source files. MPSS is explicitly out of scope (Packet 1).
    %
    % ---- PACKET 2 REPAIR (2026-07-09, user-directed, on Opus) ---------------
    % This version supersedes the first Packet 2 pass, which set the HPSS
    % existence_mask to ~failure_mask alone. That was unsafe for Packet 3:
    % on the real source file 1000/1595 (62.7%) of "solved" HPSS cells are in
    % fact flat / HSS-collapsed profiles (amp < 0.5), which are NOT boundary
    % peaks and would silently contaminate an "HPSS+HPSS" sheet. The repair
    % makes HPSS existence MORPHOLOGY-FIRST: a cell counts as HPSS only if the
    % solver converged AND the returned profile is actually a boundary half-peak.
    %
    % Two schema notes this repair deliberately changes vs. atlas_plan.md:
    %  * atlas_plan.md's "Minimum vs extended schema" note said to DEFER
    %    flux_values / flux_source_location to Packets 3-4. The user explicitly
    %    directed adding branch-level flux provenance here (task item 2), so
    %    those fields ARE emitted now. (Follow-up: sync that note in
    %    atlas_plan.md; not edited here to keep this pass adapter-only.)
    %  * Stability fields now include stability_source and are documented as
    %    ANNOTATION ONLY -- they never enter any existence mask (Invariant 3).
    %
    % Returns a report struct with per-type provenance, counts, and the
    % assertion results (Level 1/2 of validation_protocol.md), so it can be
    % driven from a smoke test, not just read off the console.
    %
    % ========================================================================
    % FIELD-MAPPING TABLE
    % ========================================================================
    %
    % HPSS  (source: bvp_sweep_full_results.mat style)
    %   branch_type            <- fixed 'HPSS'
    %   construction           <- fixed 'direct_BVP'
    %   PAK_DIFF, GEF_DIFF     <- copied unchanged
    %   solver_existence_mask  <- ~failure_mask (the solver converged here)
    %   branch_class           <- res x res string: per solver-solved cell one of
    %                             'HSS' (amp<0.5, collapsed), 'HPSS' (boundary
    %                             half-peak), 'MPSS_LIKE' (single interior peak),
    %                             'MULTIPEAK' (multiple significant crowns,
    %                             including two boundary peaks separated by an
    %                             interior trough), 'UNKNOWN'
    %                             (polarized but non-boundary / ill-formed);
    %                             'FAILED' where the solver did not converge.
    %                             Classified from the profile u(x), never from
    %                             which solver produced it (Invariant 1).
    %   existence_mask         <- solver_existence_mask & branch_class=='HPSS'
    %                             (MORPHOLOGY-FIRST; excludes amp<0.5 by design)
    %   mask_hss_collapsed / mask_mpss_like / mask_multipeak / mask_unknown /
    %   mask_failed            <- res x res logical diagnostics; together with
    %                             existence_mask they PARTITION the grid.
    %   SOL                    <- SOL_grid bvp4c struct at HPSS cells; [] elsewhere
    %                             (Invariant 2: empty where not this branch type)
    %   peak_x                 <- sol.x(argmax(u)) at HPSS cells; NaN elsewhere
    %   amp                    <- max(u)-min(u) at HPSS cells; NaN elsewhere
    %                             (recomputed from the profile, not copied)
    %   mass_residual          <- max(|trapz(u+q+v+w)-nG|,|trapz(p+q)-nP|) at HPSS
    %                             cells (recomputed per Invariant 4); NaN elsewhere
    %   flux_values            <- res x res x 5, [u,q,v,w,p] at the TROUGH
    %                             (x=L_half, sol.y(:,end)) at HPSS cells; NaN elsewhere
    %   flux_source_location   <- res x res string: 'trough' at HPSS cells; '' elsewhere
    %   is_stable              <- is_stable_grid copied (ANNOTATION ONLY)
    %   is_stable_computed     <- solver_existence_mask if Phase 5 ran; else false
    %   stability_source       <- 'bvp_jacobian_eig_500mesh' or 'not_computed'
    %   par, p, L_total, M_BG, M_P, computation_id
    %                          <- copied unchanged (Invariant 5)
    %   res, source_path, source_computation_id, adapter_run_id  <- added provenance
    %
    % HSS  (source: hss_sweep_multibranch_results.mat style)
    %   branch_type            <- fixed 'HSS'
    %   construction           <- fixed 'algebraic'
    %   PAK_DIFF, GEF_DIFF     <- copied unchanged
    %   existence_mask         <- stored tracked branch OR an unambiguous
    %                             recovered raw unique root (HSS is uniform by
    %                             definition -- no morphology class needed)
    %   SOL                    <- {[u,q,v,w,p]} 1x5 at the representative branch;
    %                             [] where absent (uniform 5-vector per atlas_plan)
    %   peak_x                 <- NaN where existent (uniform: no peak location)
    %   amp                    <- 0 where existent (uniform: zero spatial amplitude)
    %   mass_residual          <- max(|L_half*(u+q+v+w)-nG|,|L_half*(p+q)-nP|)
    %                             at representative branch (closed-form)
    %   flux_values            <- res x res x 5, [u,q,v,w,p] at representative
    %                             branch where existent; NaN elsewhere
    %   flux_source_location   <- 'uniform' where existent; '' elsewhere
    %   is_stable              <- existence_mask (solved-defaults-stable, ANNOTATION)
    %   is_stable_computed     <- false(res,res) (no HSS stability check exists)
    %   stability_source       <- 'not_computed'
    %   hss_multibranch_preserved       <- true (full branch arrays carried over)
    %   representative_branch_index     <- res x res, first non-NaN k;
    %                             -1 for a recovered raw unique root, 0 absent
    %   representative_source  <- stored_tracked_branch or
    %                             recovered_raw_unique_root
    %   raw_unique_recovery_mask <- legacy cells recovered only because
    %                             roots_cell contains exactly one accepted root
    %   n_hss_branches_available        <- source tracked n_branches_grid
    %   U/Q/V/W/P_hss_branches, n_branches_grid, branch_id_grid
    %                          <- FULL multibranch arrays copied unchanged so no
    %                             coexisting HSS root is lost (Packet 3/bistable use)
    %   par, L_total, M_BG, M_P, computation_id  <- copied unchanged
    %   res, source_path, source_computation_id, adapter_run_id  <- added provenance

    if nargin < 1 || isempty(hpss_path)
        hpss_path = 'data/outputs/bvp_sweep_full_results.mat';
    end
    if nargin < 2 || isempty(hss_path)
        hss_path = 'data/outputs/hss_sweep_multibranch_results.mat';
    end
    if nargin < 3 || isempty(out_dir)
        out_dir = 'data/outputs';
    end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    adapter_run_id = datestr(now, 'yyyymmdd_HHMMSS');
    report = struct();
    report.adapter_run_id = adapter_run_id;

    fprintf('=== build_branch_library (Packet 2 adapter, morphology-first) ===\n');
    fprintf('adapter_run_id = %s\n', adapter_run_id);

    % ------------------------------------------------------------------ HPSS
    if exist(hpss_path, 'file')
        fprintf('\n--- HPSS source: %s ---\n', hpss_path);
        src = load(hpss_path);
        fprintf('  computation_id=%s  res=%d  L_total=%g  M_BG=%g  M_P=%g\n', ...
                src.computation_id, size(src.PAK_DIFF,1), src.L_total, src.M_BG, src.M_P);
        [hpss, checks] = adapt_hpss(src, hpss_path, adapter_run_id);
        out_path = fullfile(out_dir, sprintf('branch_HPSS_%s.mat', hpss.computation_id));
        save_branch_file(out_path, hpss);
        copyfile(out_path, fullfile(out_dir, 'branch_HPSS_latest.mat'));
        fprintf('  wrote %s\n  wrote %s\n', out_path, fullfile(out_dir, 'branch_HPSS_latest.mat'));
        report.hpss = struct('source_path', hpss_path, 'out_path', out_path, ...
                             'computation_id', hpss.computation_id, ...
                             'res', size(hpss.PAK_DIFF,1), 'checks', checks);
        print_checks('HPSS', checks);
    else
        warning('build_branch_library: HPSS source not found: %s (skipped)', hpss_path);
        report.hpss = struct('source_path', hpss_path, 'skipped', true);
    end

    % ------------------------------------------------------------------ HSS
    if exist(hss_path, 'file')
        fprintf('\n--- HSS source: %s ---\n', hss_path);
        src = load(hss_path);
        fprintf('  computation_id=%s  res=%d  L_total=%g  M_BG=%g  M_P=%g\n', ...
                src.computation_id, size(src.PAK_DIFF,1), src.L_total, src.M_BG, src.M_P);
        [hss, checks] = adapt_hss(src, hss_path, adapter_run_id);
        out_path = fullfile(out_dir, sprintf('branch_HSS_%s.mat', hss.computation_id));
        save_branch_file(out_path, hss);
        copyfile(out_path, fullfile(out_dir, 'branch_HSS_latest.mat'));
        fprintf('  wrote %s\n  wrote %s\n', out_path, fullfile(out_dir, 'branch_HSS_latest.mat'));
        report.hss = struct('source_path', hss_path, 'out_path', out_path, ...
                            'computation_id', hss.computation_id, ...
                            'res', size(hss.PAK_DIFF,1), 'checks', checks);
        print_checks('HSS', checks);
    else
        warning('build_branch_library: HSS source not found: %s (skipped)', hss_path);
        report.hss = struct('source_path', hss_path, 'skipped', true);
    end

    % --------------------------------------------------- cross-type grid check
    if isfield(report.hpss, 'out_path') && isfield(report.hss, 'out_path')
        a = load(report.hpss.out_path, 'PAK_DIFF', 'GEF_DIFF');
        b = load(report.hss.out_path,  'PAK_DIFF', 'GEF_DIFF');
        grids_match = isequal(a.PAK_DIFF, b.PAK_DIFF) && isequal(a.GEF_DIFF, b.GEF_DIFF);
        report.grids_match_across_types = grids_match;
        fprintf('\n--- Cross-type check ---\n');
        fprintf('  PAK_DIFF/GEF_DIFF identical between branch_HPSS and branch_HSS: %d\n', grids_match);
        if ~grids_match
            warning(['build_branch_library: HPSS and HSS grids do NOT match. Sheet ', ...
                     'stitching (Packet 3) is invalid until sources share (res,M_BG,M_P,L).']);
        end
    end

    fprintf('\n=== build_branch_library complete ===\n');
end


%% ========================================================================
%  HPSS ADAPTER  (morphology-first)
%  ========================================================================

function [out, checks] = adapt_hpss(src, source_path, adapter_run_id)
    res = size(src.PAK_DIFF, 1);
    n_GEF_bar = src.M_BG / 2;
    n_PAK_bar = src.M_P / 2;

    out.branch_type = 'HPSS';
    out.construction = 'direct_BVP';
    out.PAK_DIFF = src.PAK_DIFF;
    out.GEF_DIFF = src.GEF_DIFF;

    solver_existence_mask = ~src.failure_mask;
    out.solver_existence_mask = solver_existence_mask;

    % --- classify every solver-solved profile ---
    branch_class = strings(res, res);
    branch_class(src.failure_mask) = "FAILED";
    morphology_reason = strings(res, res);
    morphology_reason(src.failure_mask) = "solver_failed";
    peak_frac = nan(res, res);          % peak_x / domain length, for boundary assert
    n_negative_found = 0;
    for i = 1:res
        for j = 1:res
            if ~solver_existence_mask(i,j), continue; end
            sol = src.SOL_grid{i,j};
            morphology = classify_hpss_morphology(sol);
            branch_class(i,j) = morphology.class;
            morphology_reason(i,j) = morphology.reason;
            peak_frac(i,j) = morphology.peak_frac;
            if any(sol.y([1,3,5,7,9],:) < -1e-9, 'all')
                n_negative_found = n_negative_found + 1;
            end
        end
    end
    out.branch_class = branch_class;
    out.morphology_reason = morphology_reason;

    % morphology-first existence + diagnostic masks (these PARTITION the grid)
    out.existence_mask    = solver_existence_mask & (branch_class == "HPSS");
    out.mask_hss_collapsed = branch_class == "HSS";
    out.mask_mpss_like     = branch_class == "MPSS_LIKE";
    out.mask_multipeak     = branch_class == "MULTIPEAK";
    out.mask_unknown       = branch_class == "UNKNOWN";
    out.mask_failed        = branch_class == "FAILED";

    % --- per-cell branch quantities, populated only at valid HPSS cells ---
    out.SOL = cell(res, res);
    out.peak_x = nan(res, res);
    out.amp = nan(res, res);
    out.mass_residual = nan(res, res);
    out.flux_values = nan(res, res, 5);
    out.flux_source_location = strings(res, res);   % "" everywhere by default

    n_nonboundary_hpss = 0;
    for i = 1:res
        for j = 1:res
            if ~out.existence_mask(i,j), continue; end
            sol = src.SOL_grid{i,j};
            u = sol.y(1,:);
            [umax, ipk] = max(u);
            out.SOL{i,j} = sol;
            out.peak_x(i,j) = sol.x(ipk);
            out.amp(i,j) = umax - min(u);
            if peak_frac(i,j) >= 0.10, n_nonboundary_hpss = n_nonboundary_hpss + 1; end

            gef_int = trapz(sol.x, sol.y(1,:) + sol.y(3,:) + sol.y(5,:) + sol.y(7,:));
            pak_int = trapz(sol.x, sol.y(9,:) + sol.y(3,:));
            nG = n_GEF_bar + src.GEF_DIFF(i,j)/2;
            nP = n_PAK_bar + src.PAK_DIFF(i,j)/2;
            out.mass_residual(i,j) = max(abs(gef_int - nG), abs(pak_int - nP));

            % flux provenance: fast-diffusing values at the trough (x = L_half = end)
            out.flux_values(i,j,:) = sol.y([1,3,5,7,9], end);
            out.flux_source_location(i,j) = "trough";
        end
    end

    % --- stability (ANNOTATION ONLY; never gates any existence mask) ---
    skip_stability = isfield(src, 'p') && isfield(src.p, 'skip_stability') && src.p.skip_stability;
    out.is_stable = src.is_stable_grid;
    if skip_stability
        out.is_stable_computed = false(res, res);
        out.stability_source = 'not_computed';
    else
        out.is_stable_computed = solver_existence_mask;   % Phase 5 runs on all solved cells
        out.stability_source = 'bvp_jacobian_eig_500mesh';
    end

    % provenance
    out.par = src.par;
    if isfield(src, 'p'), out.p = src.p; end
    out.L_total = src.L_total;  out.M_BG = src.M_BG;  out.M_P = src.M_P;
    out.computation_id = src.computation_id;
    out.res = res;  out.source_path = source_path;
    out.source_computation_id = src.computation_id;  out.adapter_run_id = adapter_run_id;

    % ---- Level 1 / Level 2 checks ----
    ex = out.existence_mask;
    checks = struct();
    checks.grid_partitioned_by_classes = ...
        isequal(ex | out.mask_hss_collapsed | out.mask_mpss_like | ...
                out.mask_multipeak | out.mask_unknown | out.mask_failed, true(res,res)) && ...
        (sum(ex(:)) + sum(out.mask_hss_collapsed(:)) + sum(out.mask_mpss_like(:)) + ...
         sum(out.mask_multipeak(:)) + sum(out.mask_unknown(:)) + sum(out.mask_failed(:)) == res*res);
    checks.existence_excludes_amp_below_0p5 = all(out.amp(ex) >= 0.5);
    checks.valid_hpss_are_boundary_peaks = (n_nonboundary_hpss == 0);
    checks.valid_cells_have_SOL = all(~cellfun(@isempty, out.SOL(ex)));
    checks.invalid_cells_empty_SOL = all(cellfun(@isempty, out.SOL(~ex)));
    checks.invalid_cells_nan_scalars = all(isnan(out.amp(~ex))) && ...
                                       all(isnan(out.mass_residual(~ex))) && ...
                                       all(isnan(out.peak_x(~ex)));
    fv = reshape(out.flux_values, res*res, 5);
    checks.flux_finite_only_at_valid = all(isfinite(fv(ex(:),:)), 'all') && ...
                                       all(isnan(fv(~ex(:),:)), 'all');
    checks.flux_labeled_only_at_valid = all(out.flux_source_location(ex) == "trough") && ...
                                        all(out.flux_source_location(~ex) == "");
    checks.no_negative_concentrations = (n_negative_found == 0);
    checks.n_negative_concentration_cells = n_negative_found;
    checks.is_stable_is_annotation_not_existence = true; % is_stable never used in existence_mask above
    checks.grid_shape_matches = isequal(size(out.PAK_DIFF), [res res]);
    checks.n_solver_success = sum(solver_existence_mask(:));
    checks.n_valid_HPSS = sum(ex(:));
    checks.n_HSS_collapsed = sum(out.mask_hss_collapsed(:));
    checks.n_MPSS_like = sum(out.mask_mpss_like(:));
    checks.n_multipeak = sum(out.mask_multipeak(:));
    checks.n_unknown = sum(out.mask_unknown(:));
    checks.n_failed = sum(out.mask_failed(:));
    checks.n_total = res*res;
    checks.info_only_stability_was_computed = ~skip_stability;
end


%% ========================================================================
%  HSS ADAPTER
%  ========================================================================

function [out, checks] = adapt_hss(src, source_path, adapter_run_id)
    res = size(src.PAK_DIFF, 1);
    L_half = src.L_total / 2;
    n_GEF_bar = src.M_BG / 2;
    n_PAK_bar = src.M_P / 2;
    K = size(src.U_hss_branches, 3);

    out.branch_type = 'HSS';
    out.construction = 'algebraic';
    out.PAK_DIFF = src.PAK_DIFF;
    out.GEF_DIFF = src.GEF_DIFF;

    % Legacy tracked-HSS files can contain a raw, accepted unique root that was
    % omitted from U_hss_branches only because its artificial tracked ID fell
    % below max_branches during coverage pruning. Recover such a cell only when
    % roots_cell proves there is exactly ONE raw root, so no branch choice is
    % inferred. Multiple-root cells still require an explicitly stored branch.
    stored_mask = src.n_branches_grid > 0;
    raw_unique_mask = false(res, res);
    if isfield(src, 'roots_cell') && isfield(src, 'n_raw_roots_grid')
        raw_unique_mask = (src.n_raw_roots_grid == 1) & ~stored_mask;
    end
    out.existence_mask = stored_mask | raw_unique_mask;
    out.raw_unique_recovery_mask = raw_unique_mask;
    out.SOL = cell(res, res);
    out.peak_x = nan(res, res);
    out.amp = nan(res, res);
    out.mass_residual = nan(res, res);
    out.flux_values = nan(res, res, 5);
    out.flux_source_location = strings(res, res);
    out.representative_branch_index = zeros(res, res);
    out.representative_source = strings(res, res);

    n_negative_found = 0;
    for i = 1:res
        for j = 1:res
            k = find(~isnan(squeeze(src.U_hss_branches(i,j,:))), 1);
            if ~isempty(k)
                out.representative_branch_index(i,j) = k;
                out.representative_source(i,j) = "stored_tracked_branch";
                u = src.U_hss_branches(i,j,k); q = src.Q_hss_branches(i,j,k);
                v = src.V_hss_branches(i,j,k); w = src.W_hss_branches(i,j,k);
                pv = src.P_hss_branches(i,j,k);
            elseif raw_unique_mask(i,j)
                roots = src.roots_cell{i,j};
                if size(roots,1) ~= 1 || size(roots,2) ~= 5
                    error('build_branch_library:badRawUniqueRoot', ...
                        'Cell (%d,%d) marked raw-unique but roots_cell is not 1x5.', i, j);
                end
                out.representative_branch_index(i,j) = -1; % explicit raw-unique sentinel
                out.representative_source(i,j) = "recovered_raw_unique_root";
                u = roots(1,1); q = roots(1,2); v = roots(1,3);
                w = roots(1,4); pv = roots(1,5);
            else
                continue;
            end
            out.SOL{i,j} = [u, q, v, w, pv];
            out.peak_x(i,j) = NaN;      % uniform state: no peak location
            out.amp(i,j) = 0;           % uniform state: zero spatial amplitude by definition

            nG = n_GEF_bar + src.GEF_DIFF(i,j)/2;
            nP = n_PAK_bar + src.PAK_DIFF(i,j)/2;
            out.mass_residual(i,j) = max(abs(L_half*(u+q+v+w) - nG), abs(L_half*(pv+q) - nP));

            out.flux_values(i,j,:) = [u, q, v, w, pv];
            out.flux_source_location(i,j) = "uniform";

            if any([u q v w pv] < -1e-9), n_negative_found = n_negative_found + 1; end
        end
    end

    % --- stability (ANNOTATION ONLY) ---
    out.is_stable = out.existence_mask;        % solved-defaults-stable convention
    out.is_stable_computed = false(res, res);  % no HSS stability check exists anywhere
    out.stability_source = 'not_computed';

    % --- preserve ALL HSS branches (practical: small file) ---
    out.hss_multibranch_preserved = true;
    out.n_hss_branches_available = src.n_branches_grid;
    out.U_hss_branches = src.U_hss_branches;
    out.Q_hss_branches = src.Q_hss_branches;
    out.V_hss_branches = src.V_hss_branches;
    out.W_hss_branches = src.W_hss_branches;
    out.P_hss_branches = src.P_hss_branches;
    out.n_branches_grid = src.n_branches_grid;
    out.n_raw_roots_grid = get_or(src, 'n_raw_roots_grid', src.n_branches_grid);
    out.roots_cell = get_or(src, 'roots_cell', cell(res,res));
    if isfield(src, 'branch_id_grid'), out.branch_id_grid = src.branch_id_grid; end

    % provenance
    out.par = src.par;
    out.L_total = src.L_total;  out.M_BG = src.M_BG;  out.M_P = src.M_P;
    out.computation_id = src.computation_id;
    out.res = res;  out.source_path = source_path;
    out.source_computation_id = src.computation_id;  out.adapter_run_id = adapter_run_id;

    % ---- Level 1 / Level 2 checks ----
    ex = out.existence_mask;
    checks = struct();
    checks.existence_equals_stored_or_raw_unique = isequal(ex, stored_mask | raw_unique_mask);
    checks.no_ambiguous_raw_root_recovered = ...
        all(out.n_raw_roots_grid(raw_unique_mask) == 1);
    checks.valid_cells_have_SOL = all(~cellfun(@isempty, out.SOL(ex)));
    checks.invalid_cells_empty_SOL = all(cellfun(@isempty, out.SOL(~ex)));
    checks.invalid_cells_nan_scalars = all(isnan(out.mass_residual(~ex))) && all(isnan(out.peak_x(~ex)));
    checks.amp_zero_where_existent = all(out.amp(ex) == 0);
    fv = reshape(out.flux_values, res*res, 5);
    checks.flux_finite_only_at_valid = all(isfinite(fv(ex(:),:)), 'all') && all(isnan(fv(~ex(:),:)), 'all');
    checks.flux_labeled_only_at_valid = all(out.flux_source_location(ex) == "uniform") && ...
                                        all(out.flux_source_location(~ex) == "");
    checks.representative_k_set_iff_existent = ...
        isequal(out.representative_branch_index ~= 0, ex);
    checks.representative_source_set_iff_existent = ...
        isequal(out.representative_source ~= "", ex);
    checks.no_negative_concentrations = (n_negative_found == 0);
    checks.n_negative_concentration_cells = n_negative_found;
    checks.is_stable_is_annotation_not_existence = true;
    checks.grid_shape_matches = isequal(size(out.PAK_DIFF), [res res]);
    checks.multibranch_preserved = out.hss_multibranch_preserved && (K == size(out.U_hss_branches,3));
    checks.n_valid_HSS = sum(ex(:));
    checks.n_recovered_raw_unique = sum(raw_unique_mask(:));
    checks.n_cells_with_multiple_branches = sum(src.n_branches_grid(:) > 1);
    checks.n_total = res*res;
    checks.info_only_stability_was_computed = false;
end


%% ========================================================================
%  HELPERS
%  ========================================================================

function save_branch_file(out_path, s)
    save(out_path, '-struct', 's');
end

function v = get_or(s, name, default)
    if isfield(s, name), v = s.(name); else, v = default; end
end

function print_checks(label, c)
    % info_only_* fields are status, not pass/fail -> printed as INFO.
    fprintf('  [%s checks]\n', label);
    fn = fieldnames(c);
    for k = 1:numel(fn)
        val = c.(fn{k});
        is_info = startsWith(fn{k}, 'info_only_');
        if islogical(val) && isscalar(val)
            if is_info, tag = ternary(val, 'INFO (true)', 'INFO (false)');
            else,       tag = ternary(val, 'PASS', 'FAIL');
            end
            fprintf('    %-42s %s\n', fn{k}, tag);
        else
            fprintf('    %-42s %s\n', fn{k}, mat2str(val));
        end
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
