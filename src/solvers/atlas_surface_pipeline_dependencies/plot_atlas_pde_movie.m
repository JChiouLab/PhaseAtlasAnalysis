function report = plot_atlas_pde_movie(opts)
    % PLOT_ATLAS_PDE_MOVIE  Render a full-domain PDE trajectory (the pdepe trace of
    % run_atlas_pde_trace.m) into an .mp4 movie: two stacked subplots -- u & q on
    % top, v & w & p on bottom, both vs x -- one video frame per stored PDE
    % time frame (or every opts.frame_stride-th one). Part of the uncoupled approximation
    % (not used for any figure of the paper).
    %
    % Reads the raw sol_pde (nt x nx x 5) already saved by run_atlas_pde_trace.m; does not run
    % any PDE. Read-only with respect to the trace file; writes only the movie file.
    %
    % Video is raster-only (every VideoWriter profile is a raster codec), so the movie is rendered at
    % opts.resolution_scale x supersampling (default 3): the figure canvas and every line width /
    % font size are scaled together, so the extra pixels add detail. getframe() on an invisible
    % figure returns the requested pixel size, so this works headless.
    %
    % Fixed species colour code (RGB / hex):
    %   u  [255 87  55]/255   #FF5737
    %   q  [128 130 133]/255  #808285
    %   v  [107 31  255]/255  #6B1FFF
    %   w  [188 190 192]/255  #BCBEC0
    %   p  [0   167 157]/255  #00A79D
    %
    % Options (struct, optional):
    %   trace_path       - PDE trace .mat (default
    %                      data/outputs/atlas_pde_trace_latest.mat). Must
    %                      contain sol_pde (nt x nx x 5), x_full, t_span.
    %   save_dir         - default data/plots.
    %   file_name        - default 'Atlas_PDE_Movie.mp4' (fixed name, overwritten each time, so
    %                      repeated runs do not accumulate large video files; the trace it came from
    %                      is printed to the console and returned in report).
    %   resolution_scale - raster supersampling factor, >=1. Default 3 (tripled
    %                      resolution vs. the 900x750 base canvas -> 2704x2256). Figure size, line
    %                      widths and font sizes all scale together. Width/height are rounded to the
    %                      nearest multiple of 16: MATLAB's MPEG-4 writer can corrupt
    %                      frames (diagonal streaking, wrong colours) at
    %                      large sizes that are not multiples of 16 (H.264 macroblock
    %                      alignment); merely even dimensions are not sufficient.
    %   fps              - playback frames per second. Default 15.
    %   frame_stride     - use every Nth stored PDE frame (bounds render
    %                      time and movie length independent of the trace's own
    %                      time resolution). Default 1 (every frame). A
    %                      console warning (not an error) fires above 3000 movie frames,
    %                      since that many frames at 3x resolution is a real
    %                      run-time and file-size cost.
    %   quality          - VideoWriter Quality, 0-100. Default 90.
    %   visible          - always treated as 'off' whatever is passed (kept as an option
    %                      for consistency with the other plotting functions; passing 'on'
    %                      prints a note that it was ignored). Forcing 'off' alone does not
    %                      prevent getframe() from returning a HiDPI-scaled frame size on some
    %                      machines; see the per-frame size normalization below.
    %
    % Output: data/plots/Atlas_PDE_Movie.mp4 (or opts.file_name).
    %
    % Returns report: trace_path, movie_path, n_frames_written, n_frames_total,
    % frame_stride, fps, resolution_scale, frame_size_px ([height width]),
    % duration_sec.

    if nargin < 1, opts = struct(); end
    opts = def(opts, 'trace_path', 'data/outputs/atlas_pde_trace_latest.mat');
    opts = def(opts, 'save_dir', 'data/plots');
    opts = def(opts, 'file_name', 'Atlas_PDE_Movie.mp4');
    opts = def(opts, 'resolution_scale', 3);
    opts = def(opts, 'fps', 15);
    opts = def(opts, 'frame_stride', 1);
    opts = def(opts, 'quality', 90);
    opts = def(opts, 'visible', 'off');
    if opts.frame_stride < 1 || opts.frame_stride ~= round(opts.frame_stride)
        error('plot_atlas_pde_movie: opts.frame_stride must be a positive integer.');
    end
    if opts.fps <= 0
        error('plot_atlas_pde_movie: opts.fps must be positive.');
    end
    if opts.resolution_scale < 1
        error('plot_atlas_pde_movie: opts.resolution_scale must be >= 1.');
    end

    if ~exist(opts.trace_path, 'file')
        error(['plot_atlas_pde_movie: trace file not found: %s.\n', ...
               'Run run_atlas_pde_trace first (pipeline Stage 8).'], opts.trace_path);
    end
    T = load(opts.trace_path);
    req = {'sol_pde', 'x_full', 't_span'};
    miss = req(~isfield(T, req));
    if ~isempty(miss)
        error('plot_atlas_pde_movie: trace file missing fields: %s', strjoin(miss, ', '));
    end

    x_full = T.x_full(:).';
    t_span = T.t_span(:).';
    nt_full = numel(t_span);
    frame_idx = 1:opts.frame_stride:nt_full;
    if frame_idx(end) ~= nt_full, frame_idx(end+1) = nt_full; end  % always show the final state
    nt = numel(frame_idx);

    % ---- supersampled canvas size: base 900x750 x resolution_scale, rounded
    % to the nearest multiple of 16 (see round_mb16) ----
    rs = opts.resolution_scale;
    fig_w = round_mb16(900 * rs); fig_h = round_mb16(750 * rs);
    base_lw = 2; base_font_axes = 10; base_font_legend = 9; base_font_title = 13;

    fprintf('=== plot_atlas_pde_movie ===\n');
    fprintf('  trace=%s (id=%s)\n', opts.trace_path, getfield_default(T, 'computation_id', '?'));
    fprintf('  %d stored PDE frames -> %d movie frames (frame_stride=%d), fps=%g -> ~%.1fs runtime\n', ...
        nt_full, nt, opts.frame_stride, opts.fps, nt / opts.fps);
    fprintf('  resolution_scale=%g -> frame size %dx%d px (base 900x750)\n', rs, fig_w, fig_h);
    if strcmpi(opts.visible, 'on')
        fprintf(['  NOTE: opts.visible=''on'' is ignored here -- the rendering figure is always created\n', ...
                 '  invisibly. A genuinely on-screen figure''s getframe() output can be scaled by the\n', ...
                 '  display''s HiDPI/Retina backing-store factor, and that factor is not guaranteed\n', ...
                 '  constant frame-to-frame -- this can crash VideoWriter with\n', ...
                 '  "Frame must be W by H" partway through a Visible=on render. Headless figures do not\n', ...
                 '  have this problem (every validated render in this project uses one).\n']);
    end
    if nt > 3000
        warning(['plot_atlas_pde_movie: %d movie frames at %dx%d is a lot -- rendering may take a ', ...
                 'while and the file will be large. Consider opts.frame_stride > 1 or a lower ', ...
                 'opts.resolution_scale.'], nt, fig_w, fig_h);
    end

    % ---- fixed species colours (not derived from the sheet
    % palette used elsewhere in the atlas tooling -- this is a raw PDE-state
    % movie, not a phase-plane/sheet plot) ----
    col.u = [255 87  55] / 255;
    col.q = [128 130 133] / 255;
    col.v = [107 31  255] / 255;
    col.w = [188 190 192] / 255;
    col.p = [0   167 157] / 255;

    % ---- fixed axis limits over the WHOLE (unsubsampled) trace so the movie
    % doesn't jitter frame to frame, computed once up front, not per frame ----
    top_lim = pad_lims([T.sol_pde(:,:,1); T.sol_pde(:,:,2)]);
    bot_lim = pad_lims([T.sol_pde(:,:,3); T.sol_pde(:,:,4); T.sol_pde(:,:,5)]);

    if ~exist(opts.save_dir, 'dir'), mkdir(opts.save_dir); end
    movie_path = fullfile(opts.save_dir, opts.file_name);

    fig = figure('Color', 'w', 'Position', [100 100 fig_w fig_h], 'Visible', 'off');

    ax1 = subplot(2, 1, 1); hold(ax1, 'on'); box(ax1, 'on');
    h_u = plot(ax1, x_full, squeeze(T.sol_pde(1, :, 1)), '-', 'Color', col.u, 'LineWidth', base_lw * rs);
    h_q = plot(ax1, x_full, squeeze(T.sol_pde(1, :, 2)), '-', 'Color', col.q, 'LineWidth', base_lw * rs);
    legend(ax1, {'u', 'q'}, 'Location', 'eastoutside', 'FontSize', base_font_legend * rs);
    ylabel(ax1, 'concentration'); ylim(ax1, top_lim); xlim(ax1, [x_full(1), x_full(end)]);
    grid(ax1, 'on'); set(ax1, 'FontSize', base_font_axes * rs, 'LineWidth', rs);

    ax2 = subplot(2, 1, 2); hold(ax2, 'on'); box(ax2, 'on');
    h_v = plot(ax2, x_full, squeeze(T.sol_pde(1, :, 3)), '-', 'Color', col.v, 'LineWidth', base_lw * rs);
    h_w = plot(ax2, x_full, squeeze(T.sol_pde(1, :, 4)), '-', 'Color', col.w, 'LineWidth', base_lw * rs);
    h_p = plot(ax2, x_full, squeeze(T.sol_pde(1, :, 5)), '-', 'Color', col.p, 'LineWidth', base_lw * rs);
    legend(ax2, {'v', 'w', 'p'}, 'Location', 'eastoutside', 'FontSize', base_font_legend * rs);
    xlabel(ax2, 'x'); ylabel(ax2, 'concentration'); ylim(ax2, bot_lim); xlim(ax2, [x_full(1), x_full(end)]);
    grid(ax2, 'on'); set(ax2, 'FontSize', base_font_axes * rs, 'LineWidth', rs);

    st = sgtitle('', 'FontSize', base_font_title * rs);

    try
        vw = VideoWriter(movie_path, 'MPEG-4'); %#ok<TNMLP>
        vw.Quality = opts.quality;
    catch
        movie_path = regexprep(movie_path, '\.mp4$', '.avi');
        vw = VideoWriter(movie_path, 'Motion JPEG AVI');
        vw.Quality = opts.quality;
    end
    vw.FrameRate = opts.fps;
    open(vw);

    % ---- capture loop ----
    % getframe() is not trusted to return a consistent pixel size frame to frame. Forcing
    % Visible='off' (above) does not always prevent it: on some machines getframe output is
    % scaled inconsistently (a HiDPI/Retina factor, e.g. exactly 1.5x, appears on some frames
    % but not others of the same headless render), which makes VideoWriter fail
    % ("Frame must be W by H"). Every captured frame is therefore normalized to the promised
    % [fig_h fig_w] size (via resize_frame, bilinear interpolation, no toolbox required) before it
    % reaches writeVideo. A one-time console note fires the first time a resize is needed.
    frame_size_px = [fig_h, fig_w, 3];
    resized_note_shown = false;
    for k = 1:nt
        f = frame_idx(k);
        set(h_u, 'YData', squeeze(T.sol_pde(f, :, 1)));
        set(h_q, 'YData', squeeze(T.sol_pde(f, :, 2)));
        set(h_v, 'YData', squeeze(T.sol_pde(f, :, 3)));
        set(h_w, 'YData', squeeze(T.sol_pde(f, :, 4)));
        set(h_p, 'YData', squeeze(T.sol_pde(f, :, 5)));
        set(st, 'String', sprintf('Atlas PDE trajectory: t = %.4g  (frame %d/%d)', ...
            t_span(f), f, nt_full));
        drawnow;
        fr = getframe(fig);
        if size(fr.cdata, 1) ~= fig_h || size(fr.cdata, 2) ~= fig_w
            if ~resized_note_shown
                fprintf(['  NOTE: getframe() returned %dx%d instead of the requested %dx%d on frame %d --\n', ...
                         '  normalizing every frame to %dx%d to match (the display is scaling\n', ...
                         '  captures inconsistently; the output video size is unaffected).\n'], ...
                    size(fr.cdata,2), size(fr.cdata,1), fig_w, fig_h, k, fig_w, fig_h);
                resized_note_shown = true;
            end
            fr.cdata = resize_frame(fr.cdata, fig_h, fig_w);
        end
        writeVideo(vw, fr);
    end
    close(vw);
    close(fig);  % always -- the figure is always invisible (see above), so there is
                 % nothing for opts.visible='on' to have kept open for viewing

    fprintf('  wrote %s (%d frames, %dx%d px, %.1fs @ %gfps)\n', movie_path, nt, ...
        frame_size_px(2), frame_size_px(1), nt / opts.fps, opts.fps);
    report = struct('trace_path', opts.trace_path, 'movie_path', movie_path, ...
        'n_frames_written', nt, 'n_frames_total', nt_full, ...
        'frame_stride', opts.frame_stride, 'fps', opts.fps, ...
        'resolution_scale', rs, 'frame_size_px', frame_size_px(1:2), ...
        'duration_sec', nt / opts.fps);
    fprintf('=== plot_atlas_pde_movie complete ===\n');
end


% ==========================================================================
%  helpers
% ==========================================================================

function lims = pad_lims(data)
    % [min,max] over all values with a little padding so lines don't sit on
    % the axis edge; guards a degenerate (all-equal) range.
    lo = min(data(:)); hi = max(data(:));
    if ~isfinite(lo) || ~isfinite(hi)
        lims = [0 1]; return;
    end
    if hi <= lo
        pad = max(abs(lo), 1) * 0.1 + eps;
        lims = [lo - pad, hi + pad];
        return;
    end
    pad = (hi - lo) * 0.08;
    lims = [lo - pad, hi + pad];
end

function opts = def(opts, name, val)
    if ~isfield(opts, name), opts.(name) = val; end
end

function n = round_mb16(x)
    % Nearest multiple of 16. MATLAB's MPEG-4 writer can silently corrupt
    % frames (diagonal streaking, wrong colours) at large sizes unless
    % width/height are multiples of 16 (H.264 macroblock alignment); a merely
    % even dimension does not satisfy it. Applied unconditionally, since the size
    % at which the corrupting code path activates is not known.
    n = 16 * round(x / 16);
end

function out = resize_frame(img, target_h, target_w)
    % Bilinear resize of an HxWx3 uint8 image to target_h x target_w, using
    % only interp2 (base MATLAB, no Image Processing Toolbox dependency --
    % imresize is not guaranteed available). Used to normalize a getframe()
    % capture that came back an unexpected size (see the capture loop above).
    [h, w, c] = size(img);
    if h == target_h && w == target_w
        out = img; return;
    end
    [Xq, Yq] = meshgrid(linspace(1, w, target_w), linspace(1, h, target_h));
    out = zeros(target_h, target_w, c, 'uint8');
    for ch = 1:c
        out(:, :, ch) = uint8(interp2(double(img(:, :, ch)), Xq, Yq, 'linear'));
    end
end

function v = getfield_default(s, name, d)
    if isfield(s, name), v = s.(name); else, v = d; end
end
