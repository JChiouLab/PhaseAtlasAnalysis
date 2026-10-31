function report = plot_atlas_sheets(opts)
    % PLOT_ATLAS_SHEETS  Packet 4: per-sheet flux field + nullcline plots for the
    % nine HSS/HPSS/MPSS atlas sheets built by build_atlas_sheets.m (Packet 3).
    %
    % Extends plot_phase_planes_multibranch.m's conventions (opts struct,
    % NaN-safe contour guard, quiver subsampling, display-only NaN fill) to the
    % atlas-sheet schema rather than forking its HSS-multibranch-specific logic
    % (different data source: 9 named sheets, not N HSS roots). Does NOT modify
    % plot_phase_planes_multibranch.m or touch its two canonical output files
    % (HSS_PhasePlane_Multibranch.pdf, BVP_PhasePlane_Multibranch.pdf).
    %
    % Reuses the ALREADY-COMPUTED per-sheet U_flux/V_flux stored in
    % atlas_sheets_latest.mat (Packet 3's exact discrete-mirror flux, the same
    % verified baseline formula -2*Dc*(left-right) as plot_phase_planes.m's
    % calc_fluxes) -- this plotter does NOT recompute flux from branch data.
    %
    % Options (struct, optional):
    %   opts.atlas_path      - path to the atlas sheets file.
    %                          Default 'data/outputs/atlas_sheets_latest.mat'.
    %   opts.sheets           - cellstr subset of sheet names to plot, e.g.
    %                          {'HPSS_HPSS','HPSS_MPSS'}. Default: all 9.
    %   opts.show_markers    - if true, overlay 'o' on cells that exist but are
    %                          NOT sheet_stable_provisional (annotation only,
    %                          per Invariant 3 -- never gates existence).
    %                          Default true.
    %   opts.arrow_stride    - integer >=1, quiver arrow subsampling. Default 1.
    %   opts.plot_v_only     - if true, derive the GEF-flux nullcline from v
    %                          alone instead of v+w (re-derived from the sheet's
    %                          own stored flux_left_values/flux_right_values,
    %                          same raw data, different linear combination --
    %                          matches the existing opts.plot_v_only convention
    %                          in plot_phase_planes.m / _multibranch.m). U_flux
    %                          (p-based) is unaffected. Default false (uses the
    %                          sheet's stored V_flux, which is v+w).
    %   opts.save_dir        - output directory for per-sheet PDFs + overview +
    %                          combined nullcline plot. Default 'data/plots'.
    %                          Filenames are stable/repeatable
    %                          (Atlas_Sheet_<L>_<R>.pdf, Atlas_Sheets_Overview.pdf,
    %                          Atlas_Nullclines_Combined.pdf); re-running
    %                          overwrites only THIS function's own prior output,
    %                          never the multibranch plotter's files.
    %   opts.visible         - 'on'|'off' figure visibility. Default 'off'
    %                          (headless-safe; matches _multibranch.m).
    %
    % Returns report: atlas_path, computation_id, sheet_names plotted,
    % n_exists per sheet, figures_written (paths), any_empty_sheet (cellstr of
    % sheets with sheet_exists all false -- reported, not an error).

    if nargin < 1, opts = struct(); end
    if ~isfield(opts, 'atlas_path'),   opts.atlas_path = 'data/outputs/atlas_sheets_latest.mat'; end
    if ~isfield(opts, 'sheets'),        opts.sheets = {}; end   % empty = all
    if ~isfield(opts, 'show_markers'), opts.show_markers = true; end
    if ~isfield(opts, 'arrow_stride'), opts.arrow_stride = 1; end
    if ~isfield(opts, 'plot_v_only'),  opts.plot_v_only = false; end
    if ~isfield(opts, 'save_dir'),     opts.save_dir = 'data/plots'; end
    if ~isfield(opts, 'visible'),      opts.visible = 'off'; end
    if opts.arrow_stride < 1 || opts.arrow_stride ~= round(opts.arrow_stride)
        error('opts.arrow_stride must be a positive integer.');
    end

    if ~exist(opts.atlas_path, 'file')
        error(['plot_atlas_sheets: atlas sheets file not found: %s.\n', ...
               'Run build_atlas_sheets first (Packet 3).'], opts.atlas_path);
    end
    A = load(opts.atlas_path);
    fprintf('=== plot_atlas_sheets: loaded %s ===\n', opts.atlas_path);
    fprintf('  computation_id=%s  res=%d  sources: HSS=%s HPSS=%s MPSS=%s\n', ...
        char(A.computation_id), A.res, char(A.source_HSS_id), char(A.source_HPSS_id), char(A.source_MPSS_id));

    if isempty(opts.sheets)
        names_to_plot = A.sheet_names;
    else
        names_to_plot = opts.sheets;
        bad = setdiff(names_to_plot, A.sheet_names);
        if ~isempty(bad)
            error('plot_atlas_sheets: unknown sheet name(s): %s', strjoin(bad, ', '));
        end
    end

    if ~exist(opts.save_dir, 'dir'), mkdir(opts.save_dir); end

    Dc = A.par.Dc;
    sheet_colors = distinct_sheet_colors();

    figures_written = {};
    n_exists = struct();
    empty_sheets = {};

    % ---- per-sheet individual pages ----
    for k = 1:numel(names_to_plot)
        name = names_to_plot{k};
        sh = A.sheets.(name);
        n_exists.(name) = sum(sh.sheet_exists(:));
        if n_exists.(name) == 0, empty_sheets{end+1} = name; end %#ok<AGROW>

        col = sheet_colors.(name);
        fpath = fullfile(opts.save_dir, sprintf('Atlas_Sheet_%s.pdf', name));
        plot_one_sheet(sh, name, Dc, col, opts, fpath);
        figures_written{end+1} = fpath; %#ok<AGROW>
        fprintf('  [%s] sheet_exists=%d/%d  flux_assumption="%s"  -> %s\n', ...
            name, n_exists.(name), numel(sh.sheet_exists), sh.flux_assumption, fpath);

        % Provenance spot-check (Invariant 6/11): print one representative
        % existing cell's raw flux inputs so the number is reproducible by hand.
        idx = find(sh.sheet_exists, 1);
        if ~isempty(idx)
            [i, j] = ind2sub(size(sh.sheet_exists), idx);
            fprintf(['      spot-check (%d,%d): p_left=%.4f p_right=%.4f  vw_left=%.4f vw_right=%.4f  ', ...
                     'flux_left_source=%s flux_right_source=%s  U_flux=%.4f V_flux=%.4f\n'], ...
                i, j, sh.flux_left_values(i,j,5), sh.flux_right_values(i,j,5), ...
                sh.flux_left_values(i,j,3)+sh.flux_left_values(i,j,4), ...
                sh.flux_right_values(i,j,3)+sh.flux_right_values(i,j,4), ...
                sh.flux_left_source(i,j), sh.flux_right_source(i,j), sh.U_flux(i,j), sh.V_flux(i,j));
        end
    end

    % ---- 3x3 overview grid (rows=left type, cols=right type, atlas_plan.md order) ----
    overview_path = fullfile(opts.save_dir, 'Atlas_Sheets_Overview.pdf');
    plot_overview_grid(A, sheet_colors, opts, overview_path);
    figures_written{end+1} = overview_path;
    fprintf('  overview grid -> %s\n', overview_path);

    % ---- combined nullcline overlay: all plotted sheets on ONE set of axes ----
    combined_path = fullfile(opts.save_dir, 'Atlas_Nullclines_Combined.pdf');
    plot_all_nullclines(A, names_to_plot, sheet_colors, opts, combined_path);
    figures_written{end+1} = combined_path;
    fprintf('  combined nullclines -> %s\n', combined_path);

    if ~isempty(empty_sheets)
        fprintf('  NOTE: sheet(s) with zero existing cells (plotted as an empty frame, no crash): %s\n', ...
            strjoin(empty_sheets, ', '));
    end

    report = struct();
    report.atlas_path = opts.atlas_path;
    report.computation_id = A.computation_id;
    report.sheet_names = names_to_plot;
    report.n_exists = n_exists;
    report.figures_written = figures_written;
    report.empty_sheets = empty_sheets;
    fprintf('=== plot_atlas_sheets complete: %d figure(s) written ===\n', numel(figures_written));
end


% ==========================================================================
%  ONE SHEET
% ==========================================================================

function plot_one_sheet(sh, name, Dc, col, opts, fpath)
    ex = sh.sheet_exists;
    [U_flux, V_flux] = sheet_flux(sh, Dc, opts.plot_v_only);

    fig = figure('Name', sprintf('Atlas Sheet: %s', name), 'Color', 'w', ...
                 'Position', [80 80 900 700], 'Visible', opts.visible);
    hold on; box on;

    % 1. existence-region shading (same contourf-backdrop technique as the
    % bistable-wedge shading in plot_phase_planes_multibranch.m)
    if any(ex(:))
        contourf(sh.PAK_DIFF, sh.GEF_DIFF, double(ex), [0.5 0.5], ...
                 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.12, ...
                 'HandleVisibility', 'off');
    end

    % 2. quiver (normalized), NaN-safe
    legend_handles = gobjects(0); legend_labels = {};
    if any(~isnan(U_flux(:)))
        mag = sqrt(U_flux.^2 + V_flux.^2);
        mag(mag == 0 | isnan(mag)) = 1;
        [Xq, Yq, Uq, Vq] = subsample_quiver(sh.PAK_DIFF, sh.GEF_DIFF, ...
                                            U_flux ./ mag, V_flux ./ mag, opts.arrow_stride);
        quiver(Xq, Yq, Uq, Vq, 0.5, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    end

    % 3. nullclines (solid = PAK-flux=0, dashed = GEF-flux=0), NaN-guarded
    h_pak = safe_contour_h(sh.PAK_DIFF, sh.GEF_DIFF, U_flux, [0 0], col, 2, '-');
    h_gef = safe_contour_h(sh.PAK_DIFF, sh.GEF_DIFF, V_flux, [0 0], col, 2, '--');
    h_pak_ok = ~isempty(h_pak) && isgraphics(h_pak);
    h_gef_ok = ~isempty(h_gef) && isgraphics(h_gef);
    if h_pak_ok
        legend_handles(end+1) = h_pak; legend_labels{end+1} = 'PAK-flux=0 (U_{flux})'; %#ok<AGROW>
    end
    if h_gef_ok
        legend_handles(end+1) = h_gef; legend_labels{end+1} = 'GEF-flux=0 (V_{flux})'; %#ok<AGROW>
    end

    % 4. stability annotation markers (Invariant 3: never gates existence,
    % annotation only -- 'o' on cells that exist but are not provisionally stable)
    if opts.show_markers
        unstable = ex & ~sh.sheet_stable_provisional;
        if any(unstable(:))
            h_uns = plot(sh.PAK_DIFF(unstable), sh.GEF_DIFF(unstable), 'o', ...
                 'MarkerEdgeColor', [0.9 0.2 0.2], 'MarkerSize', 4, 'LineWidth', 1.2);
            legend_handles(end+1) = h_uns; legend_labels{end+1} = 'exists, not provisionally stable'; %#ok<AGROW>
        end
    end

    title(sprintf('Atlas Sheet: %s + %s   (%d/%d cells exist)', ...
        sh.left_type, sh.right_type, sum(ex(:)), numel(ex)), 'Interpreter', 'none');
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}');
    grid on; axis tight;
    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'eastoutside');
    end

    cap = 'Solid = PAK-flux nullcline (U_{flux}=0)   |   Dashed = GEF-flux nullcline (V_{flux}=0)';
    if ~isempty(sh.flux_assumption)
        cap = sprintf('%s\nflux_assumption: %s', cap, sh.flux_assumption);
    end
    text(0.02, -0.14, cap, 'Units', 'normalized', 'FontSize', 8, 'Color', [0.3 0.3 0.3], 'Interpreter', 'none');

    try
        arrayfun(@(ax) set(ax.Toolbar, 'Visible', 'off'), findall(fig, 'Type', 'axes'));
    catch
    end
    try
        exportgraphics(fig, fpath, 'ContentType', 'vector');
    catch
        saveas(fig, fpath);
    end
    if strcmp(opts.visible, 'off'), close(fig); end
end


% ==========================================================================
%  OVERVIEW GRID (3x3, rows=left type, cols=right type)
% ==========================================================================

function plot_overview_grid(A, sheet_colors, opts, fpath)
    types = {'HSS', 'HPSS', 'MPSS'};
    fig = figure('Name', 'Atlas Sheets Overview', 'Color', 'w', ...
                 'Position', [60 60 1400 1300], 'Visible', opts.visible);
    Dc = A.par.Dc;
    for r = 1:3
        for c = 1:3
            name = sprintf('%s_%s', types{r}, types{c});
            sh = A.sheets.(name);
            ex = sh.sheet_exists;
            [U_flux, V_flux] = sheet_flux(sh, Dc, opts.plot_v_only);
            subplot(3, 3, (r-1)*3 + c); hold on; box on;
            if any(ex(:))
                contourf(sh.PAK_DIFF, sh.GEF_DIFF, double(ex), [0.5 0.5], ...
                         'FaceColor', sheet_colors.(name), 'EdgeColor', 'none', ...
                         'FaceAlpha', 0.15, 'HandleVisibility', 'off');
            end
            safe_contour_h(sh.PAK_DIFF, sh.GEF_DIFF, U_flux, [0 0], sheet_colors.(name), 1.2, '-');
            safe_contour_h(sh.PAK_DIFF, sh.GEF_DIFF, V_flux, [0 0], sheet_colors.(name), 1.2, '--');
            title(sprintf('%s + %s (%d)', types{r}, types{c}, sum(ex(:))), ...
                'Interpreter', 'none', 'FontSize', 9);
            if r == 3, xlabel('\Delta n_{PAK}'); end
            if c == 1, ylabel('\Delta n_{GEF}'); end
            grid on; axis tight;
        end
    end
    try
        sgtitle('Atlas Sheets Overview (rows=left type, cols=right type; solid=PAK-flux=0, dashed=GEF-flux=0)', ...
            'Interpreter', 'none');
    catch
    end
    try
        arrayfun(@(ax) set(ax.Toolbar, 'Visible', 'off'), findall(fig, 'Type', 'axes'));
    catch
    end
    try
        exportgraphics(fig, fpath, 'ContentType', 'vector');
    catch
        saveas(fig, fpath);
    end
    if strcmp(opts.visible, 'off'), close(fig); end
end


% ==========================================================================
%  COMBINED NULLCLINE OVERLAY (all sheets, one set of 2-D axes, no surfaces)
% ==========================================================================

function plot_all_nullclines(A, names_to_plot, sheet_colors, opts, fpath)
    Dc = A.par.Dc;
    fig = figure('Name', 'Atlas Nullclines Combined', 'Color', 'w', ...
                 'Position', [80 80 950 750], 'Visible', opts.visible);
    hold on; box on;
    legend_handles = gobjects(0); legend_labels = {};
    for k = 1:numel(names_to_plot)
        name = names_to_plot{k};
        sh = A.sheets.(name);
        [U_flux, V_flux] = sheet_flux(sh, Dc, opts.plot_v_only);
        col = sheet_colors.(name);
        h_pak = safe_contour_h(sh.PAK_DIFF, sh.GEF_DIFF, U_flux, [0 0], col, 1.8, '-');
        h_gef = safe_contour_h(sh.PAK_DIFF, sh.GEF_DIFF, V_flux, [0 0], col, 1.8, '--');
        h_pak_ok = ~isempty(h_pak) && isgraphics(h_pak);
        h_gef_ok = ~isempty(h_gef) && isgraphics(h_gef);
        % One legend entry per sheet (whichever handle drew first), matching
        % plot_phase_planes_multibranch.m's one-entry-per-branch convention.
        if h_pak_ok || h_gef_ok
            if h_pak_ok, h_for_legend = h_pak; else, h_for_legend = h_gef; end
            legend_handles(end+1) = h_for_legend; %#ok<AGROW>
            legend_labels{end+1} = strrep(name, '_', '+'); %#ok<AGROW>
        end
    end
    title('Atlas Nullclines -- All Sheets Combined');
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}');
    grid on; axis tight;
    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'eastoutside');
    end
    text(0.02, -0.12, ...
         'Solid = PAK-flux nullcline (U_{flux}=0)   |   Dashed = GEF-flux nullcline (V_{flux}=0)', ...
         'Units', 'normalized', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

    try
        arrayfun(@(ax) set(ax.Toolbar, 'Visible', 'off'), findall(fig, 'Type', 'axes'));
    catch
    end
    try
        exportgraphics(fig, fpath, 'ContentType', 'vector');
    catch
        saveas(fig, fpath);
    end
    if strcmp(opts.visible, 'off'), close(fig); end
end


% ==========================================================================
%  FLUX SELECTION (stored U_flux/V_flux, or v-only re-derivation from the
%  same stored raw per-species values -- never recomputed from branch files)
% ==========================================================================

function [U_flux, V_flux] = sheet_flux(sh, Dc, plot_v_only)
    U_flux = sh.U_flux;
    if ~plot_v_only
        V_flux = sh.V_flux;   % stored value is already v+w (Packet 3 default)
        return;
    end
    p_left  = sh.flux_left_values(:,:,3);   % v alone (row 3 of [u,q,v,w,p])
    p_right = sh.flux_right_values(:,:,3);
    V_flux = -2 * Dc * (p_left - p_right);
    V_flux(~sh.sheet_exists) = NaN;
end

function colors = distinct_sheet_colors()
    palette = { [0,167,157]/255, [230,97,1]/255, [197,27,138]/255, ...
                [33,102,172]/255, [102,37,6]/255, [77,175,74]/255, ...
                [152,78,163]/255, [255,127,0]/255, [166,86,40]/255 };
    names = {'HSS_HSS','HSS_HPSS','HSS_MPSS','HPSS_HSS','HPSS_HPSS', ...
             'HPSS_MPSS','MPSS_HSS','MPSS_HPSS','MPSS_MPSS'};
    colors = struct();
    for k = 1:numel(names)
        colors.(names{k}) = palette{k};
    end
end


% ==========================================================================
%  NaN-SAFE PLOTTING PRIMITIVES (reused verbatim from
%  plot_phase_planes_multibranch.m's proven pattern -- copied, not forked,
%  since these are local/unexported helpers; guards the exact
%  gobjects(0)/isgraphics crash fixed earlier this session)
% ==========================================================================

function h = safe_contour_h(X, Y, Z, lvls, col, lw, ls)
    h = gobjects(0);
    if all(isnan(Z(:))), return; end
    [~, h] = contour(X, Y, Z, lvls, 'Color', col, 'LineWidth', lw, 'LineStyle', ls);
end

function [Xs, Ys, Us, Vs] = subsample_quiver(X, Y, U, V, n)
    if n <= 1
        Xs = X; Ys = Y; Us = U; Vs = V;
        return;
    end
    rs = 1:n:size(X, 1);
    cs = 1:n:size(X, 2);
    Xs = X(rs, cs); Ys = Y(rs, cs);
    Us = U(rs, cs); Vs = V(rs, cs);
end
