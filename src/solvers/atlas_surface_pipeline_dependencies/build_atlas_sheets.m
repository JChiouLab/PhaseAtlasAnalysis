function report = build_atlas_sheets(hpss_path, hss_path, mpss_path, out_dir)
    % BUILD_ATLAS_SHEETS  Packet 3: stitch left/right branch pairs into the nine
    % HSS/HPSS/MPSS atlas sheets.
    %
    % Pure aggregation over the three branch-library files -- NO new BVPs are
    % solved, NO branch file is modified (read-only via load()), and NO missing
    % cell is fabricated. A sheet cell's flux/stability fields are populated
    % ONLY where sheet_exists; NaN/false/"" elsewhere (Invariant 2).
    %
    % Mirror convention (verified 2026-07-10, exact to machine zero against the
    % actual on-disk grids): the shared grid is built everywhere in this codebase
    % via `[PAK_DIFF,GEF_DIFF] = meshgrid(deln_pak, deln_G)` with deln_pak,
    % deln_G symmetric linspaces about 0. So PAK_DIFF depends only on column j,
    % GEF_DIFF only on row i, and the mass-partition mirror point (-dP,-dG) of
    % grid cell (i,j) is EXACTLY grid cell (res+1-i, res+1-j) -- pure index
    % arithmetic, no interpolation, no lookup-by-value. Implemented as
    % X(end:-1:1, end:-1:1, :), which is its own exact involution.
    %
    % sheet_exists(i,j) = left.existence_mask(i,j) & right.existence_mask(mirror(i,j))
    % (atlas_plan.md, Existence and Stability Masks). Stability is annotation
    % only (Invariant 3): sheet_stable_provisional = sheet_exists &
    % half_stable_left & half_stable_right -- never a bare "sheet_stable" claim.
    %
    % Flux baseline (verified HPSS+HPSS case, plot_phase_planes.m calc_fluxes,
    % index-mirrored instead of griddedInterpolant-mirrored so no value is
    % interpolated):
    %   U_flux = -2*Dc*(p_left - p_right),  V_flux = -2*Dc*(vw_left - vw_right)
    % where p/vw are read directly from each branch's own flux_values (already
    % populated by Packet 1D/2 at the branch's own correct reference location --
    % trough for HPSS, uniform value for HSS, outer trough at the half-domain
    % interface for MPSS) and flux_left_source/flux_right_source are copied
    % through verbatim from the branch's own flux_source_location (provenance
    % preserved, not relabelled).
    %
    % Args (all optional; default to the canonical branch_*_latest.mat):
    %   hpss_path, hss_path, mpss_path, out_dir
    %
    % Output: data/outputs/atlas_sheets_<id>.mat (+ atlas_sheets_latest.mat).
    % Top-level fields: PAK_DIFF, GEF_DIFF, res, L_total, M_BG, M_P, par,
    % computation_id, source_HSS_path/HPSS_path/MPSS_path, source_HSS_id/
    % HPSS_id/MPSS_id, sheet_names (1x9 cellstr), sheets (1x1 struct with the
    % 9 fields '<LEFT>_<RIGHT>', each holding the Required Output Fields:
    % left_type, right_type, PAK_DIFF, GEF_DIFF, sheet_exists, half_stable_left,
    % half_stable_right, sheet_stable_provisional, U_flux, V_flux,
    % left_type_used, right_type_used, flux_left_source, flux_right_source,
    % flux_left_values, flux_right_values, flux_assumption, par, L_total, M_BG,
    % M_P, computation_id, source_left_id, source_right_id).

    if nargin < 1 || isempty(hpss_path), hpss_path = 'data/outputs/branch_HPSS_latest.mat'; end
    if nargin < 2 || isempty(hss_path),  hss_path  = 'data/outputs/branch_HSS_latest.mat';  end
    if nargin < 3 || isempty(mpss_path), mpss_path = 'data/outputs/branch_MPSS_latest.mat'; end
    if nargin < 4 || isempty(out_dir),   out_dir   = 'data/outputs'; end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    computation_id = datestr(now, 'yyyymmdd_HHMMSS');
    fprintf('=== build_atlas_sheets (Packet 3) ===\n');

    % ---- load branch libraries, READ-ONLY (load() only, never save() to these) ----
    types = {'HSS', 'HPSS', 'MPSS'};
    paths = struct('HSS', hss_path, 'HPSS', hpss_path, 'MPSS', mpss_path);
    B = struct();
    for k = 1:3
        t = types{k};
        if ~exist(paths.(t), 'file')
            error('build_atlas_sheets:missing_branch', 'Required branch file missing: %s (%s)', t, paths.(t));
        end
        B.(t) = load(paths.(t));
        fprintf('  loaded %-4s %s  (id=%s, res=%d, valid=%d/%d)\n', t, paths.(t), ...
            char(B.(t).computation_id), B.(t).res, sum(B.(t).existence_mask(:)), numel(B.(t).existence_mask));
    end

    % ---- grid + parameter compatibility (hard requirement -- sheets need an exact shared grid) ----
    checks = struct();
    checks.grid_shape_matches = isequal(size(B.HSS.PAK_DIFF), size(B.HPSS.PAK_DIFF)) && ...
                                isequal(size(B.HSS.PAK_DIFF), size(B.MPSS.PAK_DIFF));
    checks.grid_values_match = isequal(B.HSS.PAK_DIFF, B.HPSS.PAK_DIFF) && isequal(B.HSS.PAK_DIFF, B.MPSS.PAK_DIFF) && ...
                               isequal(B.HSS.GEF_DIFF, B.HPSS.GEF_DIFF) && isequal(B.HSS.GEF_DIFF, B.MPSS.GEF_DIFF);
    checks.res_matches = (B.HSS.res == B.HPSS.res) && (B.HSS.res == B.MPSS.res);
    checks.par_shared_params_match = ...
        (B.HSS.par.Dc == B.HPSS.par.Dc) && (B.HSS.par.Dc == B.MPSS.par.Dc) && ...
        (B.HSS.M_BG == B.HPSS.M_BG) && (B.HSS.M_BG == B.MPSS.M_BG) && ...
        (B.HSS.M_P  == B.HPSS.M_P)  && (B.HSS.M_P  == B.MPSS.M_P)  && ...
        (B.HSS.L_total == B.HPSS.L_total) && (B.HSS.L_total == B.MPSS.L_total);
    if ~(checks.grid_shape_matches && checks.grid_values_match && checks.res_matches && checks.par_shared_params_match)
        error('build_atlas_sheets:grid_mismatch', ...
            ['Branch grids/parameters are not compatible (grid_shape=%d, grid_values=%d, res=%d, ', ...
             'par=%d). Sheets require an exact shared grid; refusing to stitch mismatched branches.'], ...
            checks.grid_shape_matches, checks.grid_values_match, checks.res_matches, checks.par_shared_params_match);
    end
    fprintf('  grid compatibility: shape=OK values=OK res=OK shared-par(Dc,M_BG,M_P,L)=OK\n');

    PAK_DIFF = B.HSS.PAK_DIFF; GEF_DIFF = B.HSS.GEF_DIFF;
    res = B.HSS.res; L_total = B.HSS.L_total; M_BG = B.HSS.M_BG; M_P = B.HSS.M_P;
    par = B.HSS.par;

    % ---- Invariant-3 input sanity: is_stable must already be false wherever
    % existence is false in every source branch (informational; source data is
    % trusted but verified, not assumed) ----
    for k = 1:3
        t = types{k};
        fname = sprintf('source_invariant3_%s', t);
        checks.(fname) = all(~B.(t).is_stable(~B.(t).existence_mask));
    end

    % ---- build all nine pure sheets ----
    sheet_names = {'HSS_HSS','HSS_HPSS','HSS_MPSS', ...
                   'HPSS_HSS','HPSS_HPSS','HPSS_MPSS', ...
                   'MPSS_HSS','MPSS_HPSS','MPSS_MPSS'};
    sheets = struct();
    n_valid = struct();
    for k = 1:numel(sheet_names)
        parts = strsplit(sheet_names{k}, '_');
        Lt = parts{1}; Rt = parts{2};
        sh = build_one_sheet(B.(Lt), Lt, paths.(Lt), B.(Rt), Rt, paths.(Rt), ...
                             par, L_total, M_BG, M_P, res, computation_id);
        sheets.(sheet_names{k}) = sh;
        n_valid.(sheet_names{k}) = sum(sh.sheet_exists(:));
        fprintf('  sheet %-12s sheet_exists=%5d/%d\n', sheet_names{k}, n_valid.(sheet_names{k}), res*res);
    end
    checks.all_nine_sheets_present = numel(fieldnames(sheets)) == 9;
    checks.n_valid = n_valid;

    % ---- mirror-index self-consistency (should be exact by construction) ----
    PM = PAK_DIFF(end:-1:1, end:-1:1); GM = GEF_DIFF(end:-1:1, end:-1:1);
    checks.mirror_involution_exact = max(abs(PAK_DIFF(:) + PM(:))) == 0 && max(abs(GEF_DIFF(:) + GM(:))) == 0;

    % ---- swap-symmetry checks: sheet_AB.exists(i,j) == sheet_BA.exists(mirror(i,j)) ----
    % (atlas_plan.md: sheets 2/4, 3/7, 6/8 -- HSS+HPSS/HPSS+HSS, HSS+MPSS/MPSS+HSS, HPSS+MPSS/MPSS+HPSS)
    swap_pairs = {'HSS_HPSS','HPSS_HSS'; 'HSS_MPSS','MPSS_HSS'; 'HPSS_MPSS','MPSS_HPSS'};
    for r = 1:size(swap_pairs, 1)
        a = swap_pairs{r,1}; b = swap_pairs{r,2};
        exA = sheets.(a).sheet_exists;
        exB_mirrored = sheets.(b).sheet_exists(end:-1:1, end:-1:1);
        n_viol = sum(exA(:) ~= exB_mirrored(:));
        checks.(sprintf('swap_symmetry_%s_vs_%s_violations', a, b)) = n_viol;
    end

    % ---- self-symmetry checks for pure same-type sheets (1,5,9) ----
    for name = {'HSS_HSS','HPSS_HPSS','MPSS_MPSS'}
        ex = sheets.(name{1}).sheet_exists;
        ex_mirrored = ex(end:-1:1, end:-1:1);
        checks.(sprintf('self_symmetry_%s_violations', name{1})) = sum(ex(:) ~= ex_mirrored(:));
    end

    % ---- no-fabrication check across all nine sheets ----
    % Per-side raw data (flux_left_values/source, flux_right_values/source)
    % must be NaN/"" exactly where its OWN branch is absent (L_ex for the left
    % side, R_ex_mirror for the right side) -- NOT wherever the combined
    % sheet_exists is false. A branch that genuinely exists everywhere (e.g.
    % HSS) legitimately keeps real, non-fabricated data on its own side even
    % when paired with a missing partner; scrubbing it to NaN there would
    % destroy real information, not prevent fabrication. Only the SHEET-level
    % outputs (U_flux, V_flux, sheet_stable_provisional) are gated on the
    % combined sheet_exists.
    no_fab_ok = true; first_bad = '';
    for k = 1:numel(sheet_names)
        parts = strsplit(sheet_names{k}, '_');
        sh = sheets.(sheet_names{k});
        ex = sh.sheet_exists;
        L_ex = B.(parts{1}).existence_mask;
        R_ex_mirror = B.(parts{2}).existence_mask(end:-1:1, end:-1:1);
        ok_flux = all(isnan(sh.U_flux(~ex))) && all(isnan(sh.V_flux(~ex)));
        ok_values = all(isnan(sh.flux_left_values(repmat(~L_ex,1,1,5)))) && ...
                    all(isnan(sh.flux_right_values(repmat(~R_ex_mirror,1,1,5))));
        ok_source = all(sh.flux_left_source(~L_ex) == "") && all(sh.flux_right_source(~R_ex_mirror) == "");
        ok_stable = all(~sh.sheet_stable_provisional(~ex));
        if ~(ok_flux && ok_values && ok_source && ok_stable) && isempty(first_bad)
            first_bad = sheet_names{k};
        end
        no_fab_ok = no_fab_ok && ok_flux && ok_values && ok_source && ok_stable;
    end
    checks.no_fabrication_all_sheets = no_fab_ok;
    checks.no_fabrication_first_violation = first_bad;

    % ---- MPSS flux-label check: verify (not assume) the required label ----
    mpss_label_ok = true;
    for k = 1:numel(sheet_names)
        parts = strsplit(sheet_names{k}, '_');
        sh = sheets.(sheet_names{k});
        if strcmp(parts{1}, 'MPSS')
            m = sh.sheet_exists;
            mpss_label_ok = mpss_label_ok && all(sh.flux_left_source(m) == "outer_trough_at_half_domain_interface");
        end
        if strcmp(parts{2}, 'MPSS')
            m = sh.sheet_exists;
            mpss_label_ok = mpss_label_ok && all(sh.flux_right_source(m) == "outer_trough_at_half_domain_interface");
        end
    end
    checks.mpss_flux_label_correct = mpss_label_ok;

    print_checks('AtlasSheets', checks);

    % ---- provenance ----
    source_HSS_path = hss_path; source_HPSS_path = hpss_path; source_MPSS_path = mpss_path; %#ok<NASGU>
    source_HSS_id = B.HSS.computation_id; source_HPSS_id = B.HPSS.computation_id; source_MPSS_id = B.MPSS.computation_id; %#ok<NASGU>

    out_path = fullfile(out_dir, sprintf('atlas_sheets_%s.mat', computation_id));
    save(out_path, 'sheets', 'sheet_names', 'PAK_DIFF', 'GEF_DIFF', 'res', 'L_total', 'M_BG', 'M_P', ...
         'par', 'computation_id', 'source_HSS_path', 'source_HPSS_path', 'source_MPSS_path', ...
         'source_HSS_id', 'source_HPSS_id', 'source_MPSS_id');
    latest_path = fullfile(out_dir, 'atlas_sheets_latest.mat');
    copyfile(out_path, latest_path);
    fprintf('  wrote %s\n  wrote %s\n', out_path, latest_path);

    report = struct();
    report.out_path = out_path;
    report.computation_id = computation_id;
    report.sheet_names = sheet_names;
    report.checks = checks;
    fprintf('=== build_atlas_sheets complete ===\n');
end


% ==========================================================================
%  ONE SHEET
% ==========================================================================

function sh = build_one_sheet(BL, Lt, Lpath, BR, Rt, Rpath, par, L_total, M_BG, M_P, res, computation_id)
    L_ex = BL.existence_mask;
    R_ex_mirror = BR.existence_mask(end:-1:1, end:-1:1);
    sheet_exists = L_ex & R_ex_mirror;

    half_stable_left  = BL.is_stable;                          % raw pass-through; already 0 where ~L_ex
    half_stable_right = BR.is_stable(end:-1:1, end:-1:1);       % raw pass-through at the mirror cell
    sheet_stable_provisional = sheet_exists & half_stable_left & half_stable_right;

    flux_left_values  = BL.flux_values;                        % res x res x 5, NaN where ~L_ex
    flux_right_values = BR.flux_values(end:-1:1, end:-1:1, :);  % mirrored, NaN where ~R_ex_mirror
    flux_left_source  = BL.flux_source_location;                % "" where ~L_ex
    flux_right_source = BR.flux_source_location(end:-1:1, end:-1:1);

    p_left  = flux_left_values(:,:,5);
    vw_left = flux_left_values(:,:,3) + flux_left_values(:,:,4);
    p_right  = flux_right_values(:,:,5);
    vw_right = flux_right_values(:,:,3) + flux_right_values(:,:,4);

    Dc = par.Dc;
    U_flux = -2 * Dc * (p_left - p_right);
    V_flux = -2 * Dc * (vw_left - vw_right);
    % Defensive explicit mask (NaN already propagates from source NaNs by
    % construction -- ~sheet_exists = ~L_ex | ~R_ex_mirror, and each term is
    % individually NaN exactly there -- but assign explicitly so a future
    % change to flux_values' NaN convention cannot silently fabricate a value).
    U_flux(~sheet_exists) = NaN;
    V_flux(~sheet_exists) = NaN;

    left_type_used  = repmat(string(Lt), res, res);
    right_type_used = repmat(string(Rt), res, res);

    if strcmp(Lt, 'MPSS') || strcmp(Rt, 'MPSS')
        flux_assumption = 'MPSS outer-trough/interface -- UNVALIDATED, provisional';
    else
        flux_assumption = '';
    end

    sh = struct();
    sh.left_type = Lt; sh.right_type = Rt;
    sh.PAK_DIFF = BL.PAK_DIFF; sh.GEF_DIFF = BL.GEF_DIFF;
    sh.sheet_exists = sheet_exists;
    sh.half_stable_left = half_stable_left; sh.half_stable_right = half_stable_right;
    sh.sheet_stable_provisional = sheet_stable_provisional;
    sh.U_flux = U_flux; sh.V_flux = V_flux;
    sh.left_type_used = left_type_used; sh.right_type_used = right_type_used;
    sh.flux_left_source = flux_left_source; sh.flux_right_source = flux_right_source;
    sh.flux_left_values = flux_left_values; sh.flux_right_values = flux_right_values;
    sh.flux_assumption = flux_assumption;
    sh.par = par; sh.L_total = L_total; sh.M_BG = M_BG; sh.M_P = M_P;
    sh.computation_id = computation_id;
    sh.source_left_id = BL.computation_id; sh.source_right_id = BR.computation_id;
    sh.source_left_path = Lpath; sh.source_right_path = Rpath;
end


% ==========================================================================
%  REPORTING (mirrors build_mpss_branch.m / build_branch_library.m style)
% ==========================================================================

function print_checks(label, c)
    fprintf('  [%s checks]\n', label);
    fn = fieldnames(c);
    for k = 1:numel(fn)
        val = c.(fn{k});
        if isstruct(val)
            sfn = fieldnames(val);
            for m = 1:numel(sfn)
                fprintf('    %-46s %s\n', sprintf('%s.%s', fn{k}, sfn{m}), mat2str(val.(sfn{m})));
            end
        elseif islogical(val) && isscalar(val)
            fprintf('    %-46s %s\n', fn{k}, ternary(val, 'PASS', 'FAIL'));
        else
            fprintf('    %-46s %s\n', fn{k}, mat2str(val));
        end
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
