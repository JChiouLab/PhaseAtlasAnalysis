function report = plot_atlas_surfaces_combined(opts)
    % PLOT_ATLAS_SURFACES_COMBINED  One combined 3-D plot: the flux surfaces of the three
    % half-cell types (no-peak HSS, half-peak HPSS, full-peak MPSS) overlaid on shared axes -- the same
    % per-type mirrored-surf style as visualize_atlas_surfaces.m, layered into one
    % figure instead of three -- with all nine sheet nullclines draped on
    % top at their correct height.
    %
    % Note on the paper: the flux here is the uncoupled approximation of build_atlas_sheets.m
    % (-2*Dc*(trough difference)), and the surfaces and nullclines are not used for any figure of the
    % paper (Fig. 3B is built from the coupled quasi-steady states).
    %
    % Two subplots: v+w and p (v alone if opts.plot_v_only).
    % Each subplot shows all three surfaces (semi-transparent, mirrored
    % across the origin as in visualize_atlas_surfaces.m), plus every
    % sheet's nullcline curve for that field (PAK-flux=0 solid on the p subplot,
    % GEF-flux=0 dashed on the v+w subplot), reusing the nullcline data already
    % computed in atlas_sheets_latest.mat -- nothing
    % is recomputed from the libraries except the surfaces themselves.
    %
    % Nullcline height (unambiguous): a sheet
    % <L>+<R>'s nullcline is exactly where left library L's own value at (i,j)
    % equals right library R's value at the mirror cell. So on the nullcline,
    % L's own surface height and R's mirrored surface height are identical by
    % definition; draping the curve on L's own displayed surface (which is
    % always one of the three surfaces already on the plot) is exact,
    % with no need to choose or blend a "which side" height.
    %
    % Options (struct, optional):
    %   opts.hss_path/hpss_path/mpss_path - library file paths (defaults to the
    %                          branch_*_latest.mat files).
    %   opts.atlas_path       - atlas sheets file (default
    %                          data/outputs/atlas_sheets_latest.mat).
    %   opts.interpolate_missing - display-only NaN fill of the SURFACES (not
    %                          the nullclines, which always use raw data).
    %                          Default false.
    %   opts.interp_method    - 'natural' (default) | 'linear' | 'nearest'
    %   opts.plot_v_only      - plot v alone instead of v+w. Default false.
    %   opts.show_nullclines  - if false, plot only the combined surfaces (no
    %                          curves). Default true.
    %   opts.sheets           - cellstr subset of sheet nullclines to draw.
    %                          Default: all 9.
    %   opts.save_path        - output PDF path. Default
    %                          'data/plots/Atlas_Surfaces_Combined.pdf'.
    %   opts.visible          - 'on'|'off'. Default 'off'.
    %
    % Returns report: paths loaded, sheets drawn, save_path.

    if nargin < 1, opts = struct(); end
    opts = default_opt(opts, 'hss_path',  'data/outputs/branch_HSS_latest.mat');
    opts = default_opt(opts, 'hpss_path', 'data/outputs/branch_HPSS_latest.mat');
    opts = default_opt(opts, 'mpss_path', 'data/outputs/branch_MPSS_latest.mat');
    opts = default_opt(opts, 'atlas_path', 'data/outputs/atlas_sheets_latest.mat');
    opts = default_opt(opts, 'interpolate_missing', false);
    opts = default_opt(opts, 'interp_method', 'natural');
    opts = default_opt(opts, 'plot_v_only', false);
    opts = default_opt(opts, 'show_nullclines', true);
    opts = default_opt(opts, 'sheets', {});
    opts = default_opt(opts, 'save_path', 'data/plots/Atlas_Surfaces_Combined.pdf');
    opts = default_opt(opts, 'visible', 'off');

    if opts.plot_v_only, vw_title = 'v'; else, vw_title = 'v+w'; end

    % ---- load libraries (read-only) ----
    types = {'HSS', 'HPSS', 'MPSS'};
    bpaths = struct('HSS', opts.hss_path, 'HPSS', opts.hpss_path, 'MPSS', opts.mpss_path);
    surf_colors = struct('HSS', [0.20 0.60 0.80], 'HPSS', [0.85 0.45 0.10], 'MPSS', [0.60 0.15 0.55]);
    c_mirror = [0.6 0.6 0.6];

    fprintf('=== plot_atlas_surfaces_combined ===\n');
    B = struct();
    for k = 1:3
        t = types{k};
        if ~exist(bpaths.(t), 'file')
            error('plot_atlas_surfaces_combined:missing_branch', 'Required branch file missing: %s', bpaths.(t));
        end
        B.(t) = load(bpaths.(t));
        fprintf('  loaded %-4s %s  (id=%s, valid=%d/%d)\n', t, bpaths.(t), ...
            char(B.(t).computation_id), sum(B.(t).existence_mask(:)), numel(B.(t).existence_mask));
    end
    % grid compatibility (hard requirement -- overlaying mismatched grids is meaningless)
    if ~(isequal(B.HSS.PAK_DIFF, B.HPSS.PAK_DIFF) && isequal(B.HSS.PAK_DIFF, B.MPSS.PAK_DIFF) && ...
         isequal(B.HSS.GEF_DIFF, B.HPSS.GEF_DIFF) && isequal(B.HSS.GEF_DIFF, B.MPSS.GEF_DIFF))
        error('plot_atlas_surfaces_combined:grid_mismatch', 'Branch grids are not identical; refusing to overlay.');
    end
    PAK_DIFF = B.HSS.PAK_DIFF; GEF_DIFF = B.HSS.GEF_DIFF;

    % ---- load atlas sheets (read-only), for nullclines ----
    A = [];
    if opts.show_nullclines
        if ~exist(opts.atlas_path, 'file')
            warning('plot_atlas_surfaces_combined:no_sheets', ...
                'opts.show_nullclines=true but %s not found; plotting surfaces only.', opts.atlas_path);
            opts.show_nullclines = false;
        else
            A = load(opts.atlas_path);
            if ~isequal(A.PAK_DIFF, PAK_DIFF) || ~isequal(A.GEF_DIFF, GEF_DIFF)
                error('plot_atlas_surfaces_combined:grid_mismatch', ...
                    'Atlas sheets grid does not match branch grids; refusing to overlay.');
            end
            fprintf('  loaded SHEETS %s  (id=%s)\n', opts.atlas_path, char(A.computation_id));
            if isempty(opts.sheets), sheet_list = A.sheet_names; else, sheet_list = opts.sheets; end
        end
    end
    if ~opts.show_nullclines, sheet_list = {}; end

    % ---- per-type surface data (raw, for nullcline z-lookup) ----
    VW_raw = struct(); P_raw = struct();
    for k = 1:3
        t = types{k};
        if opts.plot_v_only, VW_raw.(t) = B.(t).flux_values(:,:,3);
        else,                 VW_raw.(t) = B.(t).flux_values(:,:,3) + B.(t).flux_values(:,:,4);
        end
        VW_raw.(t)(~B.(t).existence_mask) = NaN;
        P_raw.(t) = B.(t).flux_values(:,:,5);
        P_raw.(t)(~B.(t).existence_mask) = NaN;
    end

    outdir = fileparts(opts.save_path);
    if ~isempty(outdir) && ~exist(outdir, 'dir'), mkdir(outdir); end

    fig = figure('Name', 'Atlas Surfaces Combined', 'Color', 'w', ...
                 'Position', [60 60 1500 650], 'Visible', opts.visible);

    % ==================== subplot 1: v+w ====================
    subplot(1, 2, 1); hold on; box on;
    legend_h = gobjects(0); legend_lbl = {};
    for k = 1:3
        t = types{k};
        Z = VW_raw.(t);
        if opts.interpolate_missing
            Z = fill_nans_display(PAK_DIFF, GEF_DIFF, Z, opts.interp_method);
        end
        hS = surf(PAK_DIFF, GEF_DIFF, Z, 'FaceAlpha', 0.55, 'EdgeColor', 'none', 'FaceColor', surf_colors.(t));
        surf(-PAK_DIFF, -GEF_DIFF, Z, 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'FaceColor', c_mirror, 'HandleVisibility', 'off');
        legend_h(end+1) = hS; legend_lbl{end+1} = t; %#ok<AGROW>
    end
    if opts.show_nullclines
        legend_h_nc = gobjects(0); legend_lbl_nc = {};
        for k = 1:numel(sheet_list)
            name = sheet_list{k};
            sh = A.sheets.(name);
            col = sheet_color(name);
            h = drape_nullcline(sh.V_flux, PAK_DIFF, GEF_DIFF, VW_raw.(sh.left_type), col, '--');
            if ~isempty(h)
                legend_h_nc(end+1) = h; legend_lbl_nc{end+1} = strrep(name, '_', '+'); %#ok<AGROW>
            end
        end
        legend_h = [legend_h, legend_h_nc]; legend_lbl = [legend_lbl, legend_lbl_nc];
    end
    title(sprintf('Combined %s Flux Surfaces + GEF-flux=0 Nullclines', vw_title));
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}'); zlabel('conc');
    grid on; view(3);
    if ~isempty(legend_h), legend(legend_h, legend_lbl, 'Location', 'eastoutside', 'FontSize', 7); end

    % ==================== subplot 2: p ====================
    subplot(1, 2, 2); hold on; box on;
    legend_h2 = gobjects(0); legend_lbl2 = {};
    for k = 1:3
        t = types{k};
        Z = P_raw.(t);
        if opts.interpolate_missing
            Z = fill_nans_display(PAK_DIFF, GEF_DIFF, Z, opts.interp_method);
        end
        hS = surf(PAK_DIFF, GEF_DIFF, Z, 'FaceAlpha', 0.55, 'EdgeColor', 'none', 'FaceColor', surf_colors.(t));
        surf(-PAK_DIFF, -GEF_DIFF, Z, 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'FaceColor', c_mirror, 'HandleVisibility', 'off');
        legend_h2(end+1) = hS; legend_lbl2{end+1} = t; %#ok<AGROW>
    end
    if opts.show_nullclines
        legend_h_nc2 = gobjects(0); legend_lbl_nc2 = {};
        for k = 1:numel(sheet_list)
            name = sheet_list{k};
            sh = A.sheets.(name);
            col = sheet_color(name);
            h = drape_nullcline(sh.U_flux, PAK_DIFF, GEF_DIFF, P_raw.(sh.left_type), col, '-');
            if ~isempty(h)
                legend_h_nc2(end+1) = h; legend_lbl_nc2{end+1} = strrep(name, '_', '+'); %#ok<AGROW>
            end
        end
        legend_h2 = [legend_h2, legend_h_nc2]; legend_lbl2 = [legend_lbl2, legend_lbl_nc2];
    end
    title('Combined p Flux Surfaces + PAK-flux=0 Nullclines');
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}'); zlabel('conc');
    grid on; view(3);
    if ~isempty(legend_h2), legend(legend_h2, legend_lbl2, 'Location', 'eastoutside', 'FontSize', 7); end

    try
        arrayfun(@(ax) set(ax.Toolbar, 'Visible', 'off'), findall(fig, 'Type', 'axes'));
    catch
    end
    try
        exportgraphics(fig, opts.save_path, 'ContentType', 'vector');
    catch
        saveas(fig, opts.save_path);
    end
    fprintf('  wrote %s\n', opts.save_path);
    if strcmp(opts.visible, 'off'), close(fig); end

    report = struct();
    report.hss_path = opts.hss_path; report.hpss_path = opts.hpss_path; report.mpss_path = opts.mpss_path;
    report.atlas_path = opts.atlas_path;
    report.sheets_drawn = sheet_list;
    report.save_path = opts.save_path;
    fprintf('=== plot_atlas_surfaces_combined complete ===\n');
end


% ==========================================================================
%  NULLCLINE DRAPING (contour on the sheet's own flux field, z-height from
%  the LEFT library's raw surface -- exact on the curve, see the header)
% ==========================================================================

function h = drape_nullcline(flux_field, X, Y, Zsurf, col, ls)
    h = [];
    C = safe_contourc(X(1,:), Y(:,1)', flux_field, [0 0]);
    if isempty(C), return; end
    k = 1;
    while k < size(C, 2)
        n = C(2, k);
        x_pts = C(1, k+1:k+n);
        y_pts = C(2, k+1:k+n);
        z_pts = interp2_safe(X, Y, Zsurf, x_pts, y_pts);
        hh = plot3(x_pts, y_pts, z_pts, 'Color', col, 'LineWidth', 2.2, 'LineStyle', ls);
        if isempty(h), h = hh; end   % first segment's handle represents the whole curve in the legend
        k = k + n + 1;
    end
end

function C = safe_contourc(x, y, Z, lvls)
    % NATIVE NaN -- contourc handles NaN cells correctly on its own (skips
    % them, does not error) and must NOT be given a sentinel substitute.
    %
    % Do not replace NaN by an out-of-range sentinel (as done for a single smooth outer NaN
    % boundary). The existence regions here have scattered internal NaN holes, and
    % the sentinel is always a large positive number, so every negative-valued real cell
    % adjacent to any internal hole would spuriously "cross zero" against the sentinel. That
    % produced about 2-3x extra contour points per sheet, with z-heights unrelated to the true
    % flux field and a worst-case left-vs-right-mirror height
    % mismatch of 6e-2 (the true identity is exact to machine precision). Plain contourc on the
    % raw NaN data reproduces plot_atlas_sheets.m's 2-D nullcline behaviour and removes the
    % artifacts (worst-case identity mismatch 1.1e-16).
    if all(isnan(Z(:))), C = []; return; end
    if numel(unique(Z(~isnan(Z)))) < 2, C = []; return; end
    try
        C = contourc(x, y, Z, lvls);
    catch
        C = [];
    end
end

function z = interp2_safe(X, Y, Z, xq, yq)
    % Bilinear interp2, with a nearest-neighbour fallback for query points
    % that fall in a grid cell touching a NaN corner (bilinear needs ALL 4
    % corners; a query point 98% of the way to a valid corner still returns
    % NaN otherwise). These points are the genuine, correctly-placed contour
    % endpoints right at a true existence boundary (median observed distance
    % to the nearest valid cell: ~3.6% of one grid step) -- nearest-neighbour
    % is an accurate approximation there (it always returns a
    % real stored value). Bilinear precision is unaffected everywhere
    % else; this fallback only engages where linear interpolation returns NaN.
    z = interp2(X, Y, Z, xq, yq, 'linear');
    bad = isnan(z);
    if any(bad)
        z(bad) = interp2(X, Y, Z, xq(bad), yq(bad), 'nearest');
    end
end


% ==========================================================================
%  SMALL HELPERS
% ==========================================================================

function opts = default_opt(opts, name, val)
    if ~isfield(opts, name), opts.(name) = val; end
end

function col = sheet_color(name)
    palette = { [0,167,157]/255, [230,97,1]/255, [197,27,138]/255, ...
                [33,102,172]/255, [102,37,6]/255, [77,175,74]/255, ...
                [152,78,163]/255, [255,127,0]/255, [166,86,40]/255 };
    names = {'HSS_HSS','HSS_HPSS','HSS_MPSS','HPSS_HSS','HPSS_HPSS', ...
             'HPSS_MPSS','MPSS_HSS','MPSS_HPSS','MPSS_MPSS'};
    col = palette{strcmp(names, name)};
end

function [Z, n_filled] = fill_nans_display(X, Y, Z, method)
    nan_mask = isnan(Z); n_filled = 0;
    if ~any(nan_mask(:)) || all(nan_mask(:)), return; end
    valid = ~nan_mask;
    if nnz(valid) < 3, return; end
    F = scatteredInterpolant(X(valid), Y(valid), Z(valid), method, 'none');
    Zf = F(X(nan_mask), Y(nan_mask));
    Z(nan_mask) = Zf; n_filled = sum(~isnan(Zf));
end
