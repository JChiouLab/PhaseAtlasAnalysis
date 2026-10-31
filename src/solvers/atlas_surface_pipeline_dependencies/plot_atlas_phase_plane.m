function report = plot_atlas_phase_plane(opts)
    % PLOT_ATLAS_PHASE_PLANE  Atlas equivalent of plot_phase_planes_multibranch.m:
    % one phase-plane page PER SELECTED SHEET, each showing that sheet's own
    % flux field/nullclines (reused from Packet 3/4, never recomputed) with the
    % FULL PDE trajectory overlaid -- SOLID where the trajectory was classified
    % ON THAT sheet (best-fit & on-atlas), DASHED for the rest of the
    % simulation -- so each page reads as "here is where in this sheet's phase
    % plane the trajectory sat, in the context of its whole path."
    %
    % Which sheets get a page (selection precedence, highest first):
    %   1. opts.sheets non-empty  -> exactly that explicit list.
    %   2. opts.n_top set         -> the N best-fitting sheets, RANKED ascending
    %      by the minimum spatial-profile residual each sheet ever achieves
    %      anywhere on the trajectory (the best that sheet ever fits the trace).
    %      This deliberately CAN include sheets the trajectory was never
    %      classified onto -- a top-N sheet that never agrees simply shows an
    %      all-dashed trace (no solid segment) plus its best-fit residual in the
    %      title, which is the point of asking for "top N" rather than "dwelled".
    %      Sheets that never even EXIST at any cell the trajectory visits (all-
    %      NaN residual) rank last. n_top is clamped to the 9 sheets; Inf = all.
    %   3. neither                -> auto: every sheet the trajectory was
    %      genuinely on-atlas for >=1 frame (dwelled), in visiting order.
    %      (off-atlas frames never count as "dwelling" even though
    %      best_sheet_idx still names a least-bad match there -- Packet 7.)
    %
    % Requires a classification (class_path) -- there is no meaningful per-sheet
    % split or ranking without per-frame per-sheet residuals. Modifies neither
    % the sheets nor the classification; writes only its own output files.
    %
    % Options (struct, optional):
    %   atlas_path   - default data/outputs/atlas_sheets_latest.mat
    %   trace_path   - PDE trace .mat (default data/outputs/atlas_pde_trace_latest.mat).
    %                  Only used to report provenance; the (Delta n_PAK,Delta
    %                  n_GEF) curve itself is read from the classification.
    %   class_path   - classification .mat (default
    %                  data/outputs/pde_atlas_classification_latest.mat).
    %                  REQUIRED (errors if missing).
    %   n_top        - integer >=1 (or Inf): plot the N best-fitting sheets by
    %                  min-residual rank (see precedence #2 above). Default []
    %                  (empty = use the dwelled-sheet auto default, #3).
    %                  Ignored if opts.sheets is given.
    %   sheets       - cellstr explicit subset of sheet names to plot a page for
    %                  (precedence #1). Default {} (empty).
    %   show_markers - mark the trace's start with a solid grey dot (same grey
    %                  as the "rest of simulation" segments). The end is not
    %                  marked. Default true.
    %   arrow_stride - quiver subsampling, integer >=1. Default 1.
    %   plot_v_only  - as in plot_atlas_sheets.m (re-derive GEF-flux nullcline
    %                  from v alone). Default false.
    %   save_dir     - default data/plots.
    %   visible      - 'on'|'off'. Default 'off'.
    %
    % Output: data/plots/Atlas_PhasePlane_<sheet_name>.pdf, one per selected
    % sheet (distinct from Packet 4's Atlas_Sheet_<name>.pdf, which has no
    % trace).
    %
    % Returns report: atlas_path, trace_path, class_path, n_frames,
    % selection_mode (string), plotted_sheets (cellstr, in plotted order),
    % frames_on (struct: sheet name -> #frames classified on it), best_residual
    % (struct: sheet name -> its min residual to the trace, Inf if never
    % existent), n_off_atlas_frames, figures_written (paths).

    if nargin < 1, opts = struct(); end
    opts = def(opts, 'atlas_path', 'data/outputs/atlas_sheets_latest.mat');
    opts = def(opts, 'trace_path', 'data/outputs/atlas_pde_trace_latest.mat');
    opts = def(opts, 'class_path', 'data/outputs/pde_atlas_classification_latest.mat');
    opts = def(opts, 'n_top', []);
    opts = def(opts, 'sheets', {});
    opts = def(opts, 'show_markers', true);
    opts = def(opts, 'arrow_stride', 1);
    opts = def(opts, 'plot_v_only', false);
    opts = def(opts, 'save_dir', 'data/plots');
    opts = def(opts, 'visible', 'off');
    if opts.arrow_stride < 1 || opts.arrow_stride ~= round(opts.arrow_stride)
        error('plot_atlas_phase_plane: opts.arrow_stride must be a positive integer.');
    end
    if ~isempty(opts.n_top)
        if ~(isscalar(opts.n_top) && opts.n_top >= 1 && (isinf(opts.n_top) || opts.n_top == round(opts.n_top)))
            error('plot_atlas_phase_plane: opts.n_top must be a positive integer or Inf.');
        end
    end

    if ~exist(opts.atlas_path, 'file')
        error(['plot_atlas_phase_plane: atlas sheets file not found: %s.\n', ...
               'Run build_atlas_sheets first (Packet 3).'], opts.atlas_path);
    end
    if ~exist(opts.class_path, 'file')
        error(['plot_atlas_phase_plane: classification file not found: %s.\n', ...
               'Run classify_pde_to_atlas_sheets first (Packet 5) -- a per-sheet ', ...
               'trace split requires knowing which sheet each frame belongs to.'], opts.class_path);
    end
    A = load(opts.atlas_path);
    C = load(opts.class_path);
    fprintf('=== plot_atlas_phase_plane: loaded %s ===\n', opts.atlas_path);
    fprintf('  classification: %d frames (id=%s)\n', numel(C.dpak), char(C.computation_id));
    if exist(opts.trace_path, 'file')
        fprintf('  trace provenance: %s\n', opts.trace_path);
    end

    pak = C.dpak(:).'; gef = C.dgef(:).';
    best_idx = C.best_sheet_idx(:).';
    off_atlas = C.off_atlas(:).';
    sheet_names_c = C.sheet_names;
    n = numel(pak);
    nS = numel(sheet_names_c);

    % ---- rank ALL sheets by how well they ever fit the trajectory: the
    % minimum spatial-profile residual each sheet achieves over the trace
    % (ascending = best-fitting first). residual_by_sheet is NaN at frames
    % where a sheet does not exist at the trajectory's current cell; a sheet
    % that never exists anywhere on the path has an all-NaN column and is
    % ranked last (Inf). This ranking drives the opt.n_top selection. ----
    best_resid_per_sheet = inf(1, nS);
    for s = 1:nS
        rs = C.residual_by_sheet(:, s);
        rs = rs(~isnan(rs));
        if ~isempty(rs), best_resid_per_sheet(s) = min(rs); end
    end
    [~, rank_order] = sort(best_resid_per_sheet, 'ascend');  % Inf sorts last

    % dwelled sheets (auto default): genuinely on-atlas, not the least-bad
    % off-atlas match (Packet 7: e.g. HSS+HSS is the argmin for most of a
    % moving-peak trajectory's off-atlas tail without ever being "on" it).
    on_mask = ~isnan(best_idx) & ~off_atlas;
    dwelled_order = unique(best_idx(on_mask), 'stable');

    % ---- selection precedence: explicit sheets > n_top > dwelled ----
    if ~isempty(opts.sheets)
        sheets_to_plot = opts.sheets;
        bad = setdiff(sheets_to_plot, A.sheet_names);
        if ~isempty(bad)
            error('plot_atlas_phase_plane: unknown sheet name(s): %s', strjoin(bad, ', '));
        end
        selection_mode = 'explicit sheet list';
    elseif ~isempty(opts.n_top)
        ntop = min(opts.n_top, nS);
        sheets_to_plot = sheet_names_c(rank_order(1:ntop));
        selection_mode = sprintf('top-%d by best-fit residual', ntop);
    else
        sheets_to_plot = sheet_names_c(dwelled_order);
        selection_mode = 'dwelled (on-atlas) sheets';
    end
    fprintf('  selection: %s -> %d sheet(s)\n', selection_mode, numel(sheets_to_plot));

    if isempty(sheets_to_plot)
        % only reachable in the dwelled-default mode when the trajectory was
        % never on-atlas anywhere (n_top always yields >=1 sheet).
        fprintf(['  trajectory was never on-atlas for any sheet -- nothing to plot ', ...
                 '(pass opts.n_top to force the top-ranked sheets anyway).\n']);
        report = struct('atlas_path', opts.atlas_path, 'trace_path', opts.trace_path, ...
            'class_path', opts.class_path, 'n_frames', n, 'selection_mode', selection_mode, ...
            'plotted_sheets', {{}}, 'frames_on', struct(), 'best_residual', struct(), ...
            'n_off_atlas_frames', sum(off_atlas), 'figures_written', {{}});
        return;
    end

    if ~exist(opts.save_dir, 'dir'), mkdir(opts.save_dir); end
    Dc = A.par.Dc;
    sheet_colors = distinct_sheet_colors();
    figures_written = {};
    frames_on = struct();
    best_residual = struct();

    for k = 1:numel(sheets_to_plot)
        name = sheets_to_plot{k};
        s_idx_this = find(strcmp(sheet_names_c, name), 1);
        on_this = (best_idx == s_idx_this) & ~off_atlas;
        frames_on.(name) = sum(on_this);
        bfit = best_resid_per_sheet(s_idx_this);
        best_residual.(name) = bfit;

        fpath = fullfile(opts.save_dir, sprintf('Atlas_PhasePlane_%s.pdf', name));
        plot_one_sheet_with_trace(A.sheets.(name), name, Dc, sheet_colors.(name), ...
            pak, gef, on_this, bfit, opts, fpath);
        figures_written{end+1} = fpath; %#ok<AGROW>
        fprintf('  [%s] on-sheet frames=%d/%d  best-fit residual=%s -> %s\n', ...
            name, frames_on.(name), n, resid_str(bfit), fpath);
    end

    report = struct('atlas_path', opts.atlas_path, 'trace_path', opts.trace_path, ...
        'class_path', opts.class_path, 'n_frames', n, 'selection_mode', selection_mode, ...
        'plotted_sheets', {sheets_to_plot}, 'frames_on', frames_on, ...
        'best_residual', best_residual, 'n_off_atlas_frames', sum(off_atlas), ...
        'figures_written', {figures_written});
    fprintf('=== plot_atlas_phase_plane complete: %d figure(s) written ===\n', numel(figures_written));
end


% ==========================================================================
%  ONE SHEET + trace (solid = on this sheet, dashed = rest of trajectory)
% ==========================================================================

function plot_one_sheet_with_trace(sh, name, Dc, col, pak, gef, on_this, bfit, opts, fpath)
    ex = sh.sheet_exists;
    [U_flux, V_flux] = sheet_flux(sh, Dc, opts.plot_v_only);

    fig = figure('Name', sprintf('Atlas Phase Plane: %s', name), 'Color', 'w', ...
                 'Position', [80 80 900 700], 'Visible', opts.visible);
    hold on; box on;

    % 1. existence-region shading (same technique as plot_atlas_sheets.m)
    if any(ex(:))
        contourf(sh.PAK_DIFF, sh.GEF_DIFF, double(ex), [0.5 0.5], ...
                 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.12, ...
                 'HandleVisibility', 'off');
    end

    % 2. quiver (this sheet's own flux, normalized), NaN-safe
    legend_handles = gobjects(0); legend_labels = {};
    if any(~isnan(U_flux(:)))
        mag = sqrt(U_flux.^2 + V_flux.^2);
        mag(mag == 0 | isnan(mag)) = 1;
        [Xq, Yq, Uq, Vq] = subsample_quiver(sh.PAK_DIFF, sh.GEF_DIFF, ...
                                            U_flux ./ mag, V_flux ./ mag, opts.arrow_stride);
        quiver(Xq, Yq, Uq, Vq, 0.5, 'Color', [0.75 0.75 0.75], 'HandleVisibility', 'off');
    end

    % 3. nullclines (solid=PAK-flux=0, dashed=GEF-flux=0), NaN-guarded
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

    % 4. trajectory: the FULL trace (every frame, regardless of on/off-sheet)
    % is drawn once as a continuous grey dotted background line, then the
    % on-this-sheet segments are overlaid ON TOP of it as a thick, 50%-opaque
    % black line -- so the underlying full path is always visible, with the
    % matched portion highlighted rather than replaced. Overlay segments are
    % grouped into contiguous runs of constant on_this (rather than one
    % plot() call per frame pair) for rendering efficiency; adjacent runs
    % share their boundary frame so the overlay itself stays visually
    % connected across any brief off-sheet gaps.
    n = numel(pak);
    rest_color = [0.4 0.4 0.4];   % shared by the full background trace and
                                   % the start marker below, so they read as
                                   % the same element

    h_rest = plot(pak, gef, ':', 'Color', rest_color, 'LineWidth', 1.0, 'HandleVisibility', 'off');

    run_starts = [1, find(diff(on_this) ~= 0) + 1];
    run_ends   = [run_starts(2:end) - 1, n];
    h_on = gobjects(0);
    for r = 1:numel(run_starts)
        a = run_starts(r); b = min(run_ends(r) + 1, n);
        if on_this(a)
            seg_pak = pak(a:b); seg_gef = gef(a:b);
            h_on = plot(seg_pak, seg_gef, '-', 'Color', [0 0 0 0.5], 'LineWidth', 2.5, 'HandleVisibility', 'off');
        end
    end
    if ~isempty(h_on)
        legend_handles(end+1) = h_on; %#ok<AGROW>
        legend_labels{end+1} = sprintf('trace: on this sheet (%d frames)', sum(on_this)); %#ok<AGROW>
    end
    legend_handles(end+1) = h_rest; %#ok<AGROW>
    legend_labels{end+1} = sprintf('trace: full simulation (%d frames)', n); %#ok<AGROW>

    if opts.show_markers
        % Start only, marked as a solid dot in the same grey used for
        % "rest of simulation" (off-this-sheet) segments -- not the sheet's
        % own colour, since the trace here is being read for its shape/path,
        % not to imply the start belongs to this sheet. No end marker.
        plot(pak(1), gef(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', rest_color, ...
             'MarkerEdgeColor', rest_color, 'HandleVisibility', 'off');
    end

    title(sprintf('Atlas Phase Plane: %s   (on this sheet %d/%d frames; best-fit residual %s)', ...
        strrep(name, '_', '+'), sum(on_this), n, resid_str(bfit)), 'Interpreter', 'none');
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}');
    grid on; axis tight;
    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'eastoutside', 'FontSize', 8);
    end
    cap = ['Grey dotted = full simulation trace | solid black (50% opaque), overlaid on top = ', ...
           'segment(s) on THIS sheet | grey dot = trace start (end not marked)'];
    if sum(on_this) == 0
        cap = [cap, sprintf('\n(this sheet was never the best-fit classification; shown by rank -- best it ever fit the trace is residual %s)', resid_str(bfit))];
    end
    text(0.02, -0.12, cap, 'Units', 'normalized', 'FontSize', 8, ...
         'Color', [0.3 0.3 0.3], 'Interpreter', 'none');

    try, arrayfun(@(ax) set(ax.Toolbar, 'Visible', 'off'), findall(fig, 'Type', 'axes')); catch, end
    try, exportgraphics(fig, fpath, 'ContentType', 'vector'); catch, saveas(fig, fpath); end
    if strcmp(opts.visible, 'off'), close(fig); end
end


% ==========================================================================
%  helpers (copied from plot_atlas_sheets.m's proven pattern -- these are
%  local/unexported subfunctions there, so they are duplicated here rather
%  than cross-called, consistent with this codebase's existing convention,
%  e.g. plot_atlas_surfaces_combined.m)
% ==========================================================================

function s = resid_str(bfit)
    % Human-readable best-fit residual: 'n/a (never exists on path)' for Inf.
    if isinf(bfit), s = 'n/a (never exists on path)'; else, s = sprintf('%.3f', bfit); end
end

function [U_flux, V_flux] = sheet_flux(sh, Dc, plot_v_only)
    U_flux = sh.U_flux;
    if ~plot_v_only
        V_flux = sh.V_flux;
        return;
    end
    p_left  = sh.flux_left_values(:,:,3);
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

function opts = def(opts, name, val)
    if ~isfield(opts, name), opts.(name) = val; end
end
