function report = visualize_atlas_surfaces(opts)
    % VISUALIZE_ATLAS_SURFACES  Surface visualization for the half-cell library files
    % (branch_<TYPE>_latest.mat; no-peak HSS, half-peak HPSS, full-peak MPSS).
    % Part of the uncoupled approximation; not used for any figure of the paper.
    %
    % Old-style call is supported:
    %   visualize_atlas_surfaces(struct('interpolate_missing', true, 'show_markers', false))
    %
    % What it plots (only for branch types whose file loads):
    %   * existence masks (always from the RAW mask -- never interpolated)
    %   * 3-D redistribution-flux surfaces v+w and p (from flux_values), mirrored
    %     across the origin, holes at failed cells
    %   * scalar fields amp, peak_x, mass_residual (2-D), masked to existence
    %   * nullcline overlays ONLY if atlas sheet / vector-field files already
    %     exist (data/outputs/sheet_*.mat); never constructed from the library
    %     surfaces alone. If missing, prints a message and continues.
    %
    % Options (struct, optional):
    %   interpolate_missing - display-only NaN fill of the numeric SURFACES via
    %                         scatteredInterpolant before plotting (out-of-hull
    %                         cells stay NaN). NEVER written back to any .mat, and
    %                         the existence-mask panels always show the raw mask.
    %                         Default false.
    %   interp_method       - 'natural' (default) | 'linear' | 'nearest'
    %   show_markers        - overlay 'x' on failed cells and 'o' on unstable
    %                         solved cells (only where stability was computed) and
    %                         '+' on half-peak cells of UNKNOWN shape class. Default true.
    %   plot_v_only         - plot v alone instead of v+w. Default false.
    %   fields              - cellstr of scalar fields to 2-D plot when present.
    %                         Default {'amp','peak_x','mass_residual'}.
    %   hss_path/hpss_path/mpss_path - override the default library (branch_*_latest) paths.
    %   out_dir             - where to look for sheet_*.mat (default data/outputs).
    %   save_dir            - if non-empty, export each figure as a PDF here
    %                         (for headless use). Default '' (interactive only).
    %   visible             - 'on' (default) | 'off' figure visibility.
    %
    % Returns (optional) a report struct: loaded, paths, grid_compatible,
    % comparable_types, n_figures, nullclines_drawn, sheet_message.

    if nargin < 1, opts = struct(); end
    opts = default_opt(opts, 'interpolate_missing', false);
    opts = default_opt(opts, 'interp_method', 'natural');
    opts = default_opt(opts, 'show_markers', true);
    opts = default_opt(opts, 'plot_v_only', false);
    opts = default_opt(opts, 'fields', {'amp','peak_x','mass_residual'});
    opts = default_opt(opts, 'hss_path',  'data/outputs/branch_HSS_latest.mat');
    opts = default_opt(opts, 'hpss_path', 'data/outputs/branch_HPSS_latest.mat');
    opts = default_opt(opts, 'mpss_path', 'data/outputs/branch_MPSS_latest.mat');
    opts = default_opt(opts, 'out_dir', 'data/outputs');
    opts = default_opt(opts, 'save_dir', '');
    opts = default_opt(opts, 'visible', 'on');

    if opts.plot_v_only, vw_title = 'v'; else, vw_title = 'v+w'; end

    types = {'HSS','HPSS','MPSS'};
    paths = struct('HSS', opts.hss_path, 'HPSS', opts.hpss_path, 'MPSS', opts.mpss_path);
    colors = struct('HSS', [0.20 0.60 0.80], 'HPSS', [0.85 0.45 0.10], 'MPSS', [0.60 0.15 0.55]);

    % ---- load ------------------------------------------------------------
    fprintf('\n=== visualize_atlas_surfaces: loading branch libraries ===\n');
    B = struct(); loaded = struct('HSS',false,'HPSS',false,'MPSS',false);
    for k = 1:3
        t = types{k};
        [B.(t), loaded.(t)] = load_branch(paths.(t), t);
    end
    if ~any(structfun(@(x)x, loaded))
        warning('visualize_atlas_surfaces:none_loaded', ...
            'No branch_*_latest.mat files could be loaded; nothing to plot.');
        report = emptyreport(paths);
        return;
    end

    % ---- provenance printout ------------------------------
    for k = 1:3
        t = types{k};
        if loaded.(t)
            d = B.(t);
            fprintf('  %-4s  %s\n', t, paths.(t));
            fprintf('        id=%s  res=%d  M_BG=%g  M_P=%g  L_total=%g  construction=%s  valid=%d/%d\n', ...
                char(d.computation_id), d.res, d.M_BG, d.M_P, d.L_total, char(d.construction), ...
                sum(d.existence_mask(:)), numel(d.existence_mask));
        else
            fprintf('  %-4s  MISSING (%s)\n', t, paths.(t));
        end
    end

    % ---- grid compatibility ---------------------------------------------
    fprintf('\n=== grid compatibility ===\n');
    present = types(structfun(@(x)x, loaded));
    grid_compatible = struct();
    for a = 1:numel(present)
        for b = a+1:numel(present)
            ta = present{a}; tb = present{b};
            [ok, why] = grids_match(B.(ta), B.(tb));
            grid_compatible.([ta '_' tb]) = ok;
            if ok, fprintf('  %s == %s : MATCH\n', ta, tb);
            else,  fprintf('  %s != %s : MISMATCH (%s)\n', ta, tb, why); end
        end
    end
    if numel(present) < 2, fprintf('  (only one branch loaded; nothing to cross-check)\n'); end

    % Reference grid = first loaded type; comparable = those matching it exactly.
    ref = present{1};
    comparable = {ref};
    for k = 1:numel(present)
        t = present{k};
        if strcmp(t, ref), continue; end
        key1 = [ref '_' t]; key2 = [t '_' ref];
        ok = (isfield(grid_compatible,key1) && grid_compatible.(key1)) || ...
             (isfield(grid_compatible,key2) && grid_compatible.(key2));
        if ok, comparable{end+1} = t; end %#ok<AGROW>
    end

    save_dir = opts.save_dir;
    if ~isempty(save_dir) && ~exist(save_dir,'dir'), mkdir(save_dir); end
    figs = gobjects(0);

    % ---- 1. existence masks (RAW, never interpolated) --------------------
    figs(end+1) = panel_mask(B, loaded, present, opts, colors); %#ok<AGROW>

    % ---- 2. scalar fields (amp, peak_x, mass_residual, ...) --------------
    for k = 1:numel(opts.fields)
        fn = opts.fields{k};
        f = panel_scalar(B, loaded, present, fn, opts);
        if isgraphics(f), figs(end+1) = f; end %#ok<AGROW>
    end

    % ---- 3. 3-D redistribution-flux surfaces (v+w and p), per branch -----
    for k = 1:numel(present)
        t = present{k};
        if ~isfield(B.(t), 'flux_values')
            fprintf('  [%s] no flux_values field -- skipping 3-D flux surface.\n', t);
            continue;
        end
        f = flux_surface_figure(B.(t), t, opts, colors.(t), vw_title);
        if isgraphics(f), figs(end+1) = f; end %#ok<AGROW>
    end

    % ---- 4. nullclines: only if sheet/vector-field files exist -----------
    [nullclines_drawn, sheet_message, nfig] = maybe_nullclines(opts, present, B, loaded);
    if isgraphics(nfig), figs(end+1) = nfig; end %#ok<AGROW>

    % ---- save (headless) -------------------------------------------------
    if ~isempty(save_dir)
        for i = 1:numel(figs)
            nm = get(figs(i), 'Name'); if isempty(nm), nm = sprintf('fig%d', i); end
            nm = regexprep(nm, '[^\w]+', '_');
            export_fig_pdf(figs(i), fullfile(save_dir, [nm '.pdf']));
        end
        fprintf('\nSaved %d figure(s) to %s\n', numel(figs), save_dir);
    end

    fprintf('\n=== visualize_atlas_surfaces done: %d figure(s) ===\n', numel(figs));

    if nargout > 0
        report = struct();
        report.loaded = loaded;
        report.paths = paths;
        report.grid_compatible = grid_compatible;
        report.comparable_types = {comparable};
        report.n_figures = numel(figs);
        report.nullclines_drawn = nullclines_drawn;
        report.sheet_message = sheet_message;
        report.save_dir = save_dir;
    end
end

% ======================================================================
%  Panels
% ======================================================================

function f = panel_mask(B, loaded, present, opts, colors)
    f = figure('Name','Atlas existence masks','Color','w','Visible',opts.visible, ...
               'Position',[60 60 420*max(1,numel(present)) 380]);
    for k = 1:numel(present)
        t = present{k}; d = B.(t);
        subplot(1, numel(present), k);
        surf(d.PAK_DIFF, d.GEF_DIFF, double(d.existence_mask), 'EdgeColor','none');
        view(2); axis tight; colormap(gca, [0.92 0.92 0.92; colors.(t)]); caxis([0 1]);
        xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}');
        title(sprintf('%s existence (%d/%d)', t, sum(d.existence_mask(:)), numel(d.existence_mask)));
    end
    sg(sprintf('Existence masks (raw; never interpolated)'));
end

function f = panel_scalar(B, loaded, present, fname, opts)
    % 1xN 2-D heatmaps of a scalar field, masked to existence, display-interp
    % optional. Returns [] if no loaded type has the field.
    f = [];
    avail = present(cellfun(@(t) isfield(B.(t), fname) && isnumeric(B.(t).(fname)) ...
                                 && isequal(size(B.(t).(fname)), size(B.(t).existence_mask)), present));
    if isempty(avail)
        fprintf('  field "%s" not present as a grid on any loaded branch -- skipped.\n', fname);
        return;
    end
    uselog = strcmp(fname, 'mass_residual');
    f = figure('Name',['Atlas ' fname],'Color','w','Visible',opts.visible, ...
               'Position',[80 80 420*numel(avail) 380]);
    for k = 1:numel(avail)
        t = avail{k}; d = B.(t);
        Z = double(d.(fname));
        Z(~d.existence_mask) = NaN;      % show only where the library state exists
        nfill = 0;
        if opts.interpolate_missing
            [Z, nfill] = fill_nans_display(d.PAK_DIFF, d.GEF_DIFF, Z, opts.interp_method);
        end
        if uselog, Zplot = log10(max(Z, realmin)); zl = 'log_{10}'; else, Zplot = Z; zl = ''; end
        subplot(1, numel(avail), k);
        surf(d.PAK_DIFF, d.GEF_DIFF, Zplot, 'EdgeColor','none'); view(2); axis tight; colorbar;
        xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}');
        ttl = sprintf('%s -- %s%s', t, zl, fname);
        if opts.interpolate_missing && nfill > 0, ttl = sprintf('%s (disp-fill %d)', ttl, nfill); end
        title(ttl, 'Interpreter','none');
    end
    sg(fname);
end

function f = flux_surface_figure(d, t, opts, col, vw_title)
    % 3-D flux surfaces v+w and p from flux_values (rows 1..5 = u,q,v,w,p),
    % mirrored across the origin, holes at failed cells, optional markers.
    f = figure('Name',['Atlas flux 3D ' t],'Color','w','Visible',opts.visible, ...
               'Position',[100 100 1200 500]);
    ex = d.existence_mask;
    if opts.plot_v_only
        VW = d.flux_values(:,:,3);
    else
        VW = d.flux_values(:,:,3) + d.flux_values(:,:,4);
    end
    P = d.flux_values(:,:,5);
    VW(~ex) = NaN; P(~ex) = NaN;
    nfillVW = 0; nfillP = 0;
    if opts.interpolate_missing
        [VW, nfillVW] = fill_nans_display(d.PAK_DIFF, d.GEF_DIFF, VW, opts.interp_method);
        [P,  nfillP ] = fill_nans_display(d.PAK_DIFF, d.GEF_DIFF, P,  opts.interp_method);
    end
    c_mirror = [0.6 0.6 0.6];

    subplot(1,2,1);
    surf(d.PAK_DIFF, d.GEF_DIFF, VW, 'FaceAlpha',0.6,'EdgeColor','none','FaceColor',col); hold on;
    surf(-d.PAK_DIFF, -d.GEF_DIFF, VW, 'FaceAlpha',0.25,'EdgeColor','none','FaceColor',c_mirror);
    if opts.show_markers, overlay_markers(d, VW); end
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}'); zlabel('conc'); grid on; view(3);
    ttl = sprintf('%s flux %s (%s)', t, vw_title, flux_loc_label(d));
    if opts.interpolate_missing && nfillVW>0, ttl = sprintf('%s [disp-fill %d]', ttl, nfillVW); end
    title(ttl, 'Interpreter','none');

    subplot(1,2,2);
    surf(d.PAK_DIFF, d.GEF_DIFF, P, 'FaceAlpha',0.6,'EdgeColor','none','FaceColor',col); hold on;
    surf(-d.PAK_DIFF, -d.GEF_DIFF, P, 'FaceAlpha',0.25,'EdgeColor','none','FaceColor',c_mirror);
    if opts.show_markers, overlay_markers(d, P); end
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}'); zlabel('conc'); grid on; view(3);
    ttl = sprintf('%s flux p (%s)', t, flux_loc_label(d));
    if opts.interpolate_missing && nfillP>0, ttl = sprintf('%s [disp-fill %d]', ttl, nfillP); end
    title(ttl, 'Interpreter','none');
end

function overlay_markers(d, Z)
    % 'o' unstable (only where stability computed), '+' half-peak UNKNOWN shape class.
    hold on;
    if isfield(d,'is_stable_computed') && isfield(d,'is_stable')
        comp = logical(d.is_stable_computed);
        unstable = comp & ~logical(d.is_stable) & d.existence_mask;
        if any(unstable(:))
            plot3(d.PAK_DIFF(unstable), d.GEF_DIFF(unstable), Z(unstable), 'o', ...
                  'MarkerFaceColor',[0.9 0.2 0.2],'MarkerEdgeColor','k','MarkerSize',5);
        end
    end
    if isfield(d,'mask_unknown') && any(d.mask_unknown(:))
        mu = d.mask_unknown;
        zt = zeros(nnz(mu),1);   % UNKNOWN cells have no surface value; mark at z=0 plane
        plot3(d.PAK_DIFF(mu), d.GEF_DIFF(mu), zt, '+', 'MarkerEdgeColor',[0.1 0.1 0.1], 'MarkerSize',5);
    end
end

% ======================================================================
%  Nullclines (guarded -- only from existing sheet/vector-field files)
% ======================================================================

function [drawn, msg, f] = maybe_nullclines(opts, present, B, loaded)
    drawn = false; f = gobjects(1); f = []; %#ok<NASGU>
    f = [];
    sheet_files = dir(fullfile(opts.out_dir, 'sheet_*.mat'));
    atlas_files = dir(fullfile(opts.out_dir, 'atlas_*.mat'));
    vfield_files = dir(fullfile(opts.out_dir, '*vector_field*.mat'));
    all_sheets = [sheet_files; atlas_files; vfield_files];
    fprintf('\n=== nullclines ===\n');
    if isempty(all_sheets)
        msg = sprintf(['No atlas sheet / vector-field files found in %s (looked for ', ...
            'sheet_*.mat, atlas_*.mat, *vector_field*.mat). Nullclines require ', ...
            'output of build_atlas_sheets (per-sheet U_flux/V_flux); not fabricating them from ', ...
            'branch surfaces alone. Continuing with surface plots only.'], opts.out_dir);
        fprintf('  %s\n', msg);
        return;
    end
    % Sheet files exist: draw U_flux=0 / V_flux=0 zero-contours if present.
    fprintf('  Found %d sheet/vector-field file(s); drawing nullcline overlays.\n', numel(all_sheets));
    f = figure('Name','Atlas nullclines','Color','w','Visible',opts.visible, ...
               'Position',[120 120 700 600]); hold on;
    ndrawn = 0;
    for i = 1:numel(all_sheets)
        sp = fullfile(all_sheets(i).folder, all_sheets(i).name);
        s = load(sp);
        if ~all(isfield(s, {'PAK_DIFF','GEF_DIFF','U_flux','V_flux'}))
            fprintf('  [skip] %s lacks PAK_DIFF/GEF_DIFF/U_flux/V_flux.\n', all_sheets(i).name);
            continue;
        end
        dp = s.PAK_DIFF(1,:); dg = s.GEF_DIFF(:,1);
        Uf = s.U_flux; Vf = s.V_flux;
        if isfield(s,'sheet_exists'), Uf(~s.sheet_exists) = NaN; Vf(~s.sheet_exists) = NaN; end
        contour(s.PAK_DIFF, s.GEF_DIFF, Uf, [0 0], 'LineWidth', 2, 'LineStyle','-');
        contour(s.PAK_DIFF, s.GEF_DIFF, Vf, [0 0], 'LineWidth', 2, 'LineStyle','--');
        ndrawn = ndrawn + 1;
    end
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}');
    title('Atlas nullclines (solid = PAK-flux=0, dashed = GEF-flux=0)');
    drawn = ndrawn > 0;
    msg = sprintf('Drew nullclines from %d sheet file(s).', ndrawn);
    fprintf('  %s\n', msg);
    if ~drawn && isgraphics(f), close(f); f = []; end
end

% ======================================================================
%  Small helpers
% ======================================================================

function opts = default_opt(opts, name, val)
    if ~isfield(opts, name), opts.(name) = val; end
end

function [d, ok] = load_branch(path, label)
    d = struct(); ok = false;
    if ~exist(path,'file'), return; end
    try
        d = load(path); ok = true;
    catch err
        warning('visualize_atlas_surfaces:load_failed','Failed to load %s (%s): %s', ...
                label, path, err.message);
    end
end

function [ok, why] = grids_match(da, db)
    why = '';
    if da.res ~= db.res, ok=false; why=sprintf('res %d vs %d',da.res,db.res); return; end
    if da.M_BG~=db.M_BG || da.M_P~=db.M_P || da.L_total~=db.L_total
        ok=false; why='M_BG/M_P/L_total differ'; return;
    end
    if ~isequal(da.PAK_DIFF,db.PAK_DIFF) || ~isequal(da.GEF_DIFF,db.GEF_DIFF)
        ok=false; why='PAK_DIFF/GEF_DIFF meshgrids differ'; return;
    end
    ok = true;
end

function lbl = flux_loc_label(d)
    if isfield(d,'flux_source_location')
        u = unique(d.flux_source_location(d.existence_mask));
        if isempty(u), lbl = 'n/a'; else, lbl = char(u(1)); end
    else
        lbl = 'unlabeled';
    end
end

function [Z, n_filled] = fill_nans_display(X, Y, Z, method)
    % DISPLAY-ONLY NaN fill (out-of-hull cells stay NaN). Never written back to
    % any .mat.
    nan_mask = isnan(Z); n_filled = 0;
    if ~any(nan_mask(:)) || all(nan_mask(:)), return; end
    valid = ~nan_mask;
    if nnz(valid) < 3, return; end
    F = scatteredInterpolant(X(valid), Y(valid), Z(valid), method, 'none');
    Zf = F(X(nan_mask), Y(nan_mask));
    Z(nan_mask) = Zf; n_filled = sum(~isnan(Zf));
end

function sg(txt)
    try, sgtitle(txt, 'Interpreter','none'); catch, end
end

function export_fig_pdf(fig, path)
    try
        arrayfun(@(ax) set(ax.Toolbar,'Visible','off'), findall(fig,'Type','axes'));
    catch
    end
    try
        exportgraphics(fig, path);
    catch
        try, saveas(fig, path); catch err
            warning('visualize_atlas_surfaces:save_failed','Could not save %s: %s', path, err.message);
        end
    end
end

function report = emptyreport(paths)
    report = struct('loaded', struct('HSS',false,'HPSS',false,'MPSS',false), ...
        'paths', paths, 'grid_compatible', struct(), 'comparable_types', {{}}, ...
        'n_figures', 0, 'nullclines_drawn', false, 'sheet_message', '', 'save_dir', '');
end
