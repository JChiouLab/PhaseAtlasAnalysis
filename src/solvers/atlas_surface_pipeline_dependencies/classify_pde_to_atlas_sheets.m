function report = classify_pde_to_atlas_sheets(opts)
    % CLASSIFY_PDE_TO_ATLAS_SHEETS  Packet 5: classify every frame of a PDE
    % trajectory by which atlas sheet best reconstructs its full SPATIAL profile.
    %
    % Core rule (Invariant 7): classification uses a spatial-profile residual --
    % the relative L2 norm between the PDE frame's 5-species profile and a sheet's
    % reconstructed stitched profile -- NOT phase-plane coordinate distance. The
    % (Delta n_PAK, Delta n_GEF) coordinate is used ONLY to snap the frame to the
    % grid cell whose sheets share that frame's mass partition; the decision is
    % then made on profile residuals.
    %
    % For each saved frame:
    %   1. extract the full spatial profile y_pde(x, 1:5) on x_full;
    %   2. compute its mass asymmetry (Delta n_PAK, Delta n_GEF) (same trapz
    %      definition as get_pde_trace.m) and snap to the NEAREST grid cell (i,j);
    %   3. for each of the nine sheets that EXISTS at (i,j) (sheet_exists true --
    %      HSS+HSS always exists, so there is always >=1 candidate), reconstruct
    %      its stitched profile at (i,j) via reconstruct_stitched_sheet.m (the same
    %      helper run_atlas_pde_trace.m uses to build ICs -> t=0 residual ~ 0) and
    %      compute residuals;
    %   4. best_sheet = argmin summed residual; flag AMBIGUOUS if the 2nd-best is
    %      within a (provisional) relative gap of the best; flag OFF_ATLAS if even
    %      the best residual exceeds a (provisional) threshold.
    %
    % Masks respected (Invariant 2): only sheet_exists cells are candidates;
    % reconstruction uses only the raw stored per-cell profiles (no cross-cell
    % interpolation of flux fields, no fabrication of missing cells). Missing
    % sheets at a cell are NaN in residual_by_sheet.
    %
    % NOTE ON THRESHOLDS (task_routing.md Packet 5): ambiguity_gap and off_atlas
    % thresholds are PROVISIONAL convenience flags, NOT calibrated on/off-sheet
    % decisions. The continuous residual_by_sheet / residual_by_species arrays are
    % the primary output; Packet 6 (or the user) should calibrate any hard
    % threshold from them. Defaults are deliberately loose and clearly labelled.
    %
    % Options (struct, optional):
    %   trace_path  - PDE trace .mat (default data/outputs/atlas_pde_trace_latest.mat).
    %                 Must contain sol_pde (nt x nx x 5), x_full, t_span, par,
    %                 L_total, M_BG, M_P. (Files from run_atlas_pde_trace.m
    %                 qualify; a get_pde_trace.m file also works if it carries
    %                 those fields.)
    %   atlas_path  - default data/outputs/atlas_sheets_latest.mat
    %   hss_path/hpss_path/mpss_path - branch files (defaults to *_latest.mat)
    %   ambiguity_gap - PROVISIONAL: flag ambiguous if (r2 - r1)/r1 < this.
    %                   Default 0.10.
    %   off_atlas_tol - PROVISIONAL: flag off-atlas if best residual > this.
    %                   Default 0.15.
    %   out_dir     - default data/outputs
    %   verbose     - print a per-frame line. Default false.
    %
    % Output: data/outputs/pde_atlas_classification_<id>.mat (+ _latest.mat) and a
    % returned report struct. Saved classification fields (all length nt over time):
    %   t_span, dpak, dgef, best_cell (nt x 2), residual_by_sheet (nt x 9),
    %   residual_by_species (nt x 9 x 5), best_sheet_idx, best_sheet_name,
    %   best_residual, second_best_residual, ambiguous, off_atlas, n_candidates,
    %   sheet_names, residual_formula (string), and provenance.

    if nargin < 1, opts = struct(); end
    opts = def(opts, 'trace_path', 'data/outputs/atlas_pde_trace_latest.mat');
    opts = def(opts, 'atlas_path', 'data/outputs/atlas_sheets_latest.mat');
    opts = def(opts, 'hss_path',  'data/outputs/branch_HSS_latest.mat');
    opts = def(opts, 'hpss_path', 'data/outputs/branch_HPSS_latest.mat');
    opts = def(opts, 'mpss_path', 'data/outputs/branch_MPSS_latest.mat');
    opts = def(opts, 'ambiguity_gap', 0.10);
    opts = def(opts, 'off_atlas_tol', 0.15);
    opts = def(opts, 'out_dir', 'data/outputs');
    opts = def(opts, 'verbose', false);

    computation_id = datestr(now, 'yyyymmdd_HHMMSS');
    residual_formula = ['relative L2 on full-domain 5-species spatial profile: ', ...
        'resid = ||Y_pde - Y_sheet||_F / ||Y_pde||_F (summed); ', ...
        'per species k: ||y_pde_k - y_sheet_k||_2 / ||y_pde_k||_2. ', ...
        'Sheet profile reconstructed at the frame''s nearest grid cell via ', ...
        'reconstruct_stitched_sheet.m (raw stored per-cell profiles; no cross-cell ', ...
        'interpolation, no fabrication).'];

    % ---- load ----
    if ~exist(opts.trace_path, 'file'), error('classify_pde_to_atlas_sheets: trace not found: %s', opts.trace_path); end
    T = load(opts.trace_path);
    req = {'sol_pde','x_full','t_span','par','L_total','M_BG','M_P'};
    miss = req(~isfield(T, req));
    if ~isempty(miss)
        error('classify_pde_to_atlas_sheets: trace file missing fields: %s', strjoin(miss, ', '));
    end
    A = load(opts.atlas_path);
    branches.HSS  = load(opts.hss_path);
    branches.HPSS = load(opts.hpss_path);
    branches.MPSS = load(opts.mpss_path);

    % ---- compatibility (Invariant 11) ----
    if T.L_total ~= A.L_total || T.M_BG ~= A.M_BG || T.M_P ~= A.M_P
        error(['classify_pde_to_atlas_sheets: trace geometry (L=%.3f,M_BG=%.3f,M_P=%.3f) does not ', ...
               'match atlas (L=%.3f,M_BG=%.3f,M_P=%.3f).'], ...
               T.L_total, T.M_BG, T.M_P, A.L_total, A.M_BG, A.M_P);
    end
    for bt = {'HSS','HPSS','MPSS'}
        if ~isequal(branches.(bt{1}).PAK_DIFF, A.PAK_DIFF) || ~isequal(branches.(bt{1}).GEF_DIFF, A.GEF_DIFF)
            error('classify_pde_to_atlas_sheets: branch %s grid does not match atlas grid.', bt{1});
        end
    end

    x_full = T.x_full(:).';
    L_half = T.L_total / 2;
    idx_L = x_full <= L_half; idx_R = ~idx_L;
    xL = x_full(idx_L); xR = x_full(idx_R);
    PAK_DIFF = A.PAK_DIFF; GEF_DIFF = A.GEF_DIFF;
    sheet_names = A.sheet_names; nS = numel(sheet_names);
    nt = numel(T.t_span);

    fprintf('=== classify_pde_to_atlas_sheets ===\n');
    fprintf('  trace=%s (id=%s, %d frames, nx=%d)\n', opts.trace_path, ...
        getfield_default(T,'computation_id','?'), nt, numel(x_full));
    fprintf('  atlas id=%s ; IC sheet/cell (if present): %s / %s\n', char(A.computation_id), ...
        getfield_default(T,'ic_sheet','n/a'), mat2str(getfield_default(T,'ic_cell',[])));

    % ---- per-frame classification ----
    dpak = zeros(nt,1); dgef = zeros(nt,1);
    best_cell = zeros(nt,2);
    residual_by_sheet = nan(nt, nS);
    residual_by_species = nan(nt, nS, 5);
    best_sheet_idx = nan(nt,1);
    best_residual = nan(nt,1);
    second_best_residual = nan(nt,1);
    ambiguous = false(nt,1);
    off_atlas = false(nt,1);
    n_candidates = zeros(nt,1);

    for f = 1:nt
        y_pde = squeeze(T.sol_pde(f, :, :)).';           % 5 x nx
        % mass asymmetry (same definition as get_pde_trace.m / run_atlas_pde_trace.m)
        gefL = trapz(xL, sum(y_pde(1:4, idx_L),1)); gefR = trapz(xR, sum(y_pde(1:4, idx_R),1));
        pakL = trapz(xL, y_pde(2,idx_L)+y_pde(5,idx_L)); pakR = trapz(xR, y_pde(2,idx_R)+y_pde(5,idx_R));
        dgef(f) = gefL - gefR; dpak(f) = pakL - pakR;

        % snap to nearest grid cell
        d = (PAK_DIFF - dpak(f)).^2 + (GEF_DIFF - dgef(f)).^2;
        [~, kmin] = min(d(:)); [ci, cj] = ind2sub(size(d), kmin);
        best_cell(f,:) = [ci, cj];

        npde = norm(y_pde, 'fro');
        npde_s = sqrt(sum(y_pde.^2, 2));                  % per-species norm (5x1)
        ncand = 0;
        for s = 1:nS
            y_sheet = reconstruct_stitched_sheet(A, branches, sheet_names{s}, ci, cj, x_full);
            if isempty(y_sheet), continue; end            % sheet absent at this cell
            ncand = ncand + 1;
            residual_by_sheet(f, s) = norm(y_pde - y_sheet, 'fro') / max(npde, eps);
            for k = 1:5
                residual_by_species(f, s, k) = norm(y_pde(k,:) - y_sheet(k,:)) / max(npde_s(k), eps);
            end
        end
        n_candidates(f) = ncand;

        rr = residual_by_sheet(f, :);
        [rsort, order] = sort(rr);                        % NaNs sort to the end
        best_residual(f) = rsort(1);
        best_sheet_idx(f) = order(1);
        if nS >= 2 && ~isnan(rsort(2)), second_best_residual(f) = rsort(2); end
        % provisional flags
        if ~isnan(second_best_residual(f)) && best_residual(f) > 0
            ambiguous(f) = (second_best_residual(f) - best_residual(f)) / best_residual(f) < opts.ambiguity_gap;
        end
        off_atlas(f) = best_residual(f) > opts.off_atlas_tol;

        if opts.verbose
            fprintf('  t=%.3g cell(%d,%d) best=%s r=%.4f %s%s\n', T.t_span(f), ci, cj, ...
                sheet_names{best_sheet_idx(f)}, best_residual(f), ...
                tern(ambiguous(f),'[AMBIG]',''), tern(off_atlas(f),'[OFF-ATLAS]',''));
        end
    end
    best_sheet_name = sheet_names(best_sheet_idx);

    % ---- provenance + save ----
    t_span = T.t_span(:).';
    provenance = struct();
    provenance.trace_path = opts.trace_path;
    provenance.trace_id = getfield_default(T,'computation_id','?');
    provenance.atlas_path = opts.atlas_path; provenance.atlas_id = A.computation_id;
    provenance.branch_ids = struct('HSS',branches.HSS.computation_id, ...
        'HPSS',branches.HPSS.computation_id,'MPSS',branches.MPSS.computation_id);
    provenance.ic_sheet = getfield_default(T,'ic_sheet','n/a');
    provenance.ic_cell  = getfield_default(T,'ic_cell',[]);
    provenance.ambiguity_gap = opts.ambiguity_gap;
    provenance.off_atlas_tol = opts.off_atlas_tol;
    provenance.thresholds_are_provisional = true;

    if ~exist(opts.out_dir, 'dir'), mkdir(opts.out_dir); end
    out_path = fullfile(opts.out_dir, sprintf('pde_atlas_classification_%s.mat', computation_id));
    save(out_path, 't_span', 'dpak', 'dgef', 'best_cell', 'residual_by_sheet', ...
        'residual_by_species', 'best_sheet_idx', 'best_sheet_name', 'best_residual', ...
        'second_best_residual', 'ambiguous', 'off_atlas', 'n_candidates', ...
        'sheet_names', 'residual_formula', 'provenance', 'computation_id');
    copyfile(out_path, fullfile(opts.out_dir, 'pde_atlas_classification_latest.mat'));
    fprintf('  wrote %s\n', out_path);

    % ---- console summary ----
    fprintf('  frame-0 best sheet: %s (residual %.3e -- anchor: should be ~0 if IC is an atlas sheet)\n', ...
        best_sheet_name{1}, best_residual(1));
    tab = tabulate_names(best_sheet_name);
    fprintf('  best-sheet frame counts:\n');
    for k = 1:size(tab,1), fprintf('     %-12s %d\n', tab{k,1}, tab{k,2}); end
    fprintf('  ambiguous frames: %d/%d ; off-atlas frames: %d/%d (provisional thresholds)\n', ...
        sum(ambiguous), nt, sum(off_atlas), nt);

    report = struct('out_path', out_path, 'computation_id', computation_id, ...
        'nt', nt, 'frame0_residual', best_residual(1), 'frame0_best', best_sheet_name{1}, ...
        'n_ambiguous', sum(ambiguous), 'n_off_atlas', sum(off_atlas), ...
        'best_sheet_name', {best_sheet_name}, 'residual_by_sheet', residual_by_sheet);
    fprintf('=== classify_pde_to_atlas_sheets complete ===\n');
end


% ---- helpers ----
function opts = def(opts, name, val)
    if ~isfield(opts, name), opts.(name) = val; end
end
function v = getfield_default(s, name, d)
    if isfield(s, name), v = s.(name); else, v = d; end
end
function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
function tab = tabulate_names(names)
    u = unique(names, 'stable');
    tab = cell(numel(u), 2);
    for k = 1:numel(u)
        tab{k,1} = u{k}; tab{k,2} = sum(strcmp(names, u{k}));
    end
    % sort by count desc
    [~, ord] = sort(cell2mat(tab(:,2)), 'descend');
    tab = tab(ord, :);
end
