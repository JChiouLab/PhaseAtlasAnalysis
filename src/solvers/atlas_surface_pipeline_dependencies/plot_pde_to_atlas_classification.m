function report = plot_pde_to_atlas_classification(opts)
    % PLOT_PDE_TO_ATLAS_CLASSIFICATION  Plotting helper for the output of
    % classify_pde_to_atlas_sheets.m. Part of the uncoupled approximation (not used for any
    % figure of the paper).
    %
    % Produces (one figure, saved to data/plots/):
    %   (1) best sheet over time (step plot, one lane per sheet), with ambiguous
    %       and off-atlas frames marked;
    %   (2) residual traces over time for the top sheets (lowest median residual),
    %       plus the best-per-frame residual;
    % and, optionally, a SEPARATE figure:
    %   (3) PDE-frame-vs-best-atlas-profile overlays (all 5 species) at selected
    %       times -- reconstructs the best sheet's profile the same way the
    %       classifier did (shared reconstruct_stitched_sheet.m).
    %
    % Options (struct, optional):
    %   class_path  - classification .mat (default
    %                 data/outputs/pde_atlas_classification_latest.mat).
    %   trace_path  - PDE trace .mat, needed only for the profile-overlay figure
    %                 (default data/outputs/atlas_pde_trace_latest.mat).
    %   atlas_path/hss_path/hpss_path/mpss_path - for the overlay reconstruction (library files).
    %   overlay_times - vector of times at which to draw frame-vs-atlas overlays,
    %                 or [] to skip the overlay figure. Default: [] unless the
    %                 trace is available, then 4 log-spaced-ish sample frames.
    %   n_top       - number of top sheets to show residual traces for. Default 4.
    %   save_dir    - default 'data/plots'.
    %   visible     - 'on'|'off'. Default 'off'.
    %
    % Returns report: figures_written (paths).

    if nargin < 1, opts = struct(); end
    opts = def(opts, 'class_path', 'data/outputs/pde_atlas_classification_latest.mat');
    opts = def(opts, 'trace_path', 'data/outputs/atlas_pde_trace_latest.mat');
    opts = def(opts, 'atlas_path', 'data/outputs/atlas_sheets_latest.mat');
    opts = def(opts, 'hss_path',  'data/outputs/branch_HSS_latest.mat');
    opts = def(opts, 'hpss_path', 'data/outputs/branch_HPSS_latest.mat');
    opts = def(opts, 'mpss_path', 'data/outputs/branch_MPSS_latest.mat');
    opts = def(opts, 'overlay_times', 'auto');
    opts = def(opts, 'n_top', Inf);   % Inf/>=9 => plot ALL sheets; else top-N by median residual
    opts = def(opts, 'time_log', 'auto');   % true|false|'auto' (auto=log if the trace is log-spaced)
    opts = def(opts, 'save_dir', 'data/plots');
    opts = def(opts, 'visible', 'off');

    if ~exist(opts.class_path, 'file'), error('plot_pde_to_atlas_classification: not found: %s', opts.class_path); end
    C = load(opts.class_path);
    if ~exist(opts.save_dir, 'dir'), mkdir(opts.save_dir); end
    t = C.t_span(:); nt = numel(t); nS = numel(C.sheet_names);
    figures_written = {};

    % log time axis? 'auto' = yes if the frames look log-spaced (ratio of
    % successive spacings roughly constant and spanning >1 decade).
    if ischar(opts.time_log)
        % 'auto': use a log time axis when the frames span more than ~2 decades
        % of positive time (typical for log-spaced traces), so early transients
        % are not crushed against t=0.
        tp = t(t > 0);
        use_logt = numel(t) > 5 && ~isempty(tp) && (t(end) / min(tp) > 100);
    else
        use_logt = opts.time_log;
    end

    % ============ Figure 1: best sheet over time + residual traces ============
    fig1 = figure('Name','PDE-to-Atlas classification','Color','w', ...
                  'Position',[80 80 1200 800],'Visible',opts.visible);

    % --- (1) best sheet over time ---
    subplot(2,1,1); hold on; box on;
    for f = 1:nt
        yv = C.best_sheet_idx(f);
        col = [0.2 0.4 0.8];
        if C.off_atlas(f), col = [0.85 0.1 0.1];
        elseif C.ambiguous(f), col = [0.9 0.6 0.1]; end
        plot(t(f), yv, '.', 'Color', col, 'MarkerSize', 12);
    end
    % legend proxies
    h1 = plot(nan,nan,'.','Color',[0.2 0.4 0.8],'MarkerSize',14);
    h2 = plot(nan,nan,'.','Color',[0.9 0.6 0.1],'MarkerSize',14);
    h3 = plot(nan,nan,'.','Color',[0.85 0.1 0.1],'MarkerSize',14);
    set(gca,'YTick',1:nS,'YTickLabel',strrep(C.sheet_names,'_','+'));
    ylim([0.5 nS+0.5]); xlabel(tern(use_logt,'time (log)','time')); ylabel('best sheet');
    title('Best atlas sheet per frame (by spatial-profile residual)');
    legend([h1 h2 h3], {'confident','ambiguous','off-atlas'}, 'Location','eastoutside');
    grid on; apply_logt(use_logt, t);

    % --- (2) residual traces over time ---
    % By default show ALL sheets (opts.n_top >= nS or Inf); otherwise the top
    % n_top by median residual. Each sheet is drawn in its canonical atlas colour
    % (consistent with plot_atlas_sheets.m and the combined nullcline plot) with a
    % '.'-marker so isolated existing frames of a sparse sheet still show; NaN
    % (sheet absent at that frame's cell) breaks the line -- existence-aware.
    subplot(2,1,2); hold on; box on;
    med = zeros(nS,1);
    for s = 1:nS, med(s) = median(C.residual_by_sheet(~isnan(C.residual_by_sheet(:,s)), s)); end
    med(isnan(med)) = inf;
    if opts.n_top >= nS
        showS = 1:nS;                        % all sheets, canonical order
    else
        [~, ord] = sort(med); showS = sort(ord(1:opts.n_top)).';
    end
    lh = gobjects(0); ll = {};
    for s = showS
        lh(end+1) = plot(t, C.residual_by_sheet(:,s), '.-', 'Color', sheet_palette(C.sheet_names{s}), ...
            'LineWidth', 1.4, 'MarkerSize', 8); %#ok<AGROW>
        ll{end+1} = strrep(C.sheet_names{s}, '_', '+'); %#ok<AGROW>
    end
    lh(end+1) = plot(t, C.best_residual, 'k--', 'LineWidth', 1.6); ll{end+1} = 'best-per-frame';
    if isfield(C,'provenance') && isfield(C.provenance,'off_atlas_tol')
        yline(C.provenance.off_atlas_tol, ':', 'off-atlas tol (provisional)', 'Color',[0.5 0.5 0.5]);
    end
    xlabel(tern(use_logt,'time (log)','time')); ylabel('relative L2 residual');
    if opts.n_top >= nS
        title('Spatial-profile residual of ALL sheets over time (breaks = sheet absent at that frame''s cell)');
    else
        title(sprintf('Spatial-profile residual over time (top %d sheets)', opts.n_top));
    end
    legend(lh, ll, 'Location','eastoutside', 'FontSize', 8); grid on; apply_logt(use_logt, t);

    f1path = fullfile(opts.save_dir, 'PDE_Atlas_Classification.pdf');
    export_fig(fig1, f1path); figures_written{end+1} = f1path;
    if strcmp(opts.visible,'off'), close(fig1); end

    % ============ Figure 2 (optional): frame-vs-atlas profile overlays ============
    overlay_done = false;
    if exist(opts.trace_path, 'file') && ~(ischar(opts.overlay_times) && isempty(opts.overlay_times))
        T = load(opts.trace_path);
        if isfield(T,'sol_pde')
            A = load(opts.atlas_path);
            branches.HSS = load(opts.hss_path); branches.HPSS = load(opts.hpss_path); branches.MPSS = load(opts.mpss_path);
            % choose overlay frames
            if ischar(opts.overlay_times) && strcmpi(opts.overlay_times,'auto')
                fsel = unique(round(linspace(1, nt, 4)));
            else
                fsel = arrayfun(@(tt) nearest_idx(t, tt), opts.overlay_times(:)');
                fsel = unique(fsel);
            end
            species = {'u','q','v','w','p'};
            fig2 = figure('Name','PDE frame vs best atlas profile','Color','w', ...
                          'Position',[60 60 1500 300*numel(fsel)],'Visible',opts.visible);
            for r = 1:numel(fsel)
                f = fsel(r);
                y_pde = squeeze(T.sol_pde(f,:,:)).';
                ci = C.best_cell(f,1); cj = C.best_cell(f,2);
                y_sheet = reconstruct_stitched_sheet(A, branches, C.best_sheet_name{f}, ci, cj, T.x_full);
                for k = 1:5
                    subplot(numel(fsel), 5, (r-1)*5 + k); hold on; box on;
                    plot(T.x_full, y_pde(k,:), 'b-', 'LineWidth', 1.6);
                    if ~isempty(y_sheet)
                        plot(T.x_full, y_sheet(k,:), 'r--', 'LineWidth', 1.4);
                    end
                    if k==1
                        ylabel(sprintf('t=%.3g\n%s r=%.3f', t(f), strrep(C.best_sheet_name{f},'_','+'), C.best_residual(f)));
                    end
                    if r==1, title(species{k}); end
                    if r==numel(fsel), xlabel('x'); end
                    grid on;
                end
            end
            try, sgtitle('PDE frame (blue) vs best atlas sheet reconstruction (red dashed)'); catch, end
            f2path = fullfile(opts.save_dir, 'PDE_Atlas_Frame_Overlays.pdf');
            export_fig(fig2, f2path); figures_written{end+1} = f2path;
            if strcmp(opts.visible,'off'), close(fig2); end
            overlay_done = true;
        end
    end

    report = struct('figures_written', {figures_written}, 'overlay_done', overlay_done);
    fprintf('plot_pde_to_atlas_classification: %d figure(s) written%s.\n', numel(figures_written), ...
        tern(overlay_done, '', ' (no overlay -- trace not available)'));
end


% ---- helpers ----
function opts = def(opts, name, val)
    if ~isfield(opts, name), opts.(name) = val; end
end
function s = tern(c, a, b), if c, s = a; else, s = b; end, end
function k = nearest_idx(v, x)
    [~, k] = min(abs(v - x));
end
function apply_logt(use_logt, t)
    % Log time axis; t=0 (the IC frame) cannot appear on a log axis, so start the
    % x-limit at the first positive time.
    if ~use_logt, return; end
    tp = t(t > 0);
    if isempty(tp), return; end
    set(gca, 'XScale', 'log'); xlim([tp(1), t(end)]);
end
function col = sheet_palette(name)
    % Canonical per-sheet colours, identical to plot_atlas_sheets.m's
    % distinct_sheet_colors so a sheet has the SAME colour across every atlas plot.
    palette = { [0,167,157]/255, [230,97,1]/255, [197,27,138]/255, ...
                [33,102,172]/255, [102,37,6]/255, [77,175,74]/255, ...
                [152,78,163]/255, [255,127,0]/255, [166,86,40]/255 };
    names = {'HSS_HSS','HSS_HPSS','HSS_MPSS','HPSS_HSS','HPSS_HPSS', ...
             'HPSS_MPSS','MPSS_HSS','MPSS_HPSS','MPSS_MPSS'};
    idx = find(strcmp(names, name), 1);
    if isempty(idx), col = [0.3 0.3 0.3]; else, col = palette{idx}; end
end
function export_fig(fig, path)
    try, arrayfun(@(ax) set(ax.Toolbar,'Visible','off'), findall(fig,'Type','axes')); catch, end
    try, exportgraphics(fig, path, 'ContentType','vector'); catch, saveas(fig, path); end
end
