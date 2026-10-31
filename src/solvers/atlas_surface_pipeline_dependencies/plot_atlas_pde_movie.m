function report = plot_atlas_pde_movie(opts)
    % PLOT_ATLAS_PDE_MOVIE  Render a full-domain PDE trajectory (Stage 8's
    % atlas-seeded trace) into an .mp4 movie: two stacked subplots -- u & q on
    % top, v & w & p on bottom, both vs x -- one video frame per stored PDE
    % time frame (or every opts.frame_stride-th one).
    %
    % Reads the RAW sol_pde (nt x nx x 5) already saved by run_atlas_pde_trace.m
    % (Packet 5 / Stage 8); does not run or re-run any PDE. Read-only w.r.t. the
    % trace file; writes only the movie file.
    %
    % ON VECTOR FRAMES: not possible, checked directly against
    % VideoWriter.getProfiles() -- every available profile (MPEG-4, Motion JPEG
    % AVI/2000, Archival, Uncompressed/Grayscale/Indexed AVI) is a raster codec.
    % No standard, player-compatible video container (MP4/H.264, AVI, MOV, WebM)
    % stores vector frames; video is raster-by-design. So instead this renders
    % at opts.resolution_scale x supersampling (default 3, i.e. tripled): the
    % figure canvas AND every line width / font size are scaled together by the
    % same factor, so the extra pixels add real sharpness/detail rather than
    % just shrinking everything relative to a bigger blank canvas. Confirmed
    % empirically that getframe() on an invisible figure returns cdata at
    % EXACTLY the requested Position pixel size, so this works headless too
    % (also confirmed the raw cdata is pixel-perfect even at large sizes, by
    % round-tripping it through a .mat file and rendering with Python/PIL,
    % bypassing this machine's separately-broken imwrite/print/exportgraphics-
    % to-raster path -- see round_mb16 below for the one real capture bug that
    % WAS found, in the video ENCODER rather than getframe).
    %
    % Fixed species colour code (RGB / hex, user-specified):
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
    %   file_name        - default 'Atlas_PDE_Movie.mp4' (fixed name, like Stage
    %                      10's per-sheet PDFs -- overwritten each time this
    %                      runs, so repeated pipeline runs don't accumulate
    %                      large orphaned video files; provenance -- which trace
    %                      it came from -- is printed to the console and
    %                      returned in report, not baked into the filename).
    %   resolution_scale - raster supersampling factor, >=1. Default 3 (tripled
    %                      resolution vs. the original 900x750 base canvas ->
    %                      2704x2256). Figure size, line widths, and font sizes
    %                      all scale together so the composition looks the
    %                      same, just sharper -- NOT bigger lines on a bigger
    %                      blank canvas. Width/height are rounded to the
    %                      nearest multiple of 16 (confirmed empirically that
    %                      MATLAB's MPEG-4 writer on this machine corrupts
    %                      frames -- diagonal streaking, wrong colours -- at
    %                      large non-16-aligned sizes; classic H.264 macroblock
    %                      alignment. Merely-even dimensions are NOT sufficient).
    %   fps              - playback frames per second. Default 15.
    %   frame_stride     - use every Nth stored PDE frame (bounds render
    %                      time/movie length independent of the trace's own
    %                      time resolution -- a trace can easily have thousands
    %                      of stored frames). Default 1 (every frame). A
    %                      console warning (not an error or a silent
    %                      truncation) fires above 3000 movie frames, since
    %                      that many frames at 3x resolution is a real
    %                      wall-clock and file-size cost.
    %   quality          - VideoWriter Quality, 0-100. Default 90.
    %   visible          - ALWAYS TREATED AS 'off' regardless of what is passed
    %                      (kept as an option only for interface consistency
    %                      with the other Stage 8-10 functions; passing 'on'
    %                      prints a note explaining it was ignored). NOTE:
    %                      forcing this to 'off' alone was NOT sufficient to
    %                      prevent getframe() from returning an inconsistently
    %                      HiDPI-scaled frame size on at least one real
    %                      machine (still crashed 2026-07-13 even with
    %                      Visible='off') -- see the per-frame size
    %                      normalization below, which is the actual fix.
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
               'Run run_atlas_pde_trace first (Packet 5 / pipeline Stage 8).'], opts.trace_path);
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
    % to the nearest multiple of 16 (see round_mb16 -- MPEG-4 corruption fix) ----
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
                 '  constant frame-to-frame -- confirmed 2026-07-12: this crashed VideoWriter with\n', ...
                 '  "Frame must be W by H" partway through a Visible=on render. Headless figures do not\n', ...
                 '  have this problem (every validated render in this project uses one).\n']);
    end
    if nt > 3000
        warning(['plot_atlas_pde_movie: %d movie frames at %dx%d is a lot -- rendering may take a ', ...
                 'while and the file will be large. Consider opts.frame_stride > 1 or a lower ', ...
                 'opts.resolution_scale.'], nt, fig_w, fig_h);
    end

    % ---- fixed species colours (user-specified, not derived from the sheet
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
    % IMPORTANT: getframe() is not trusted to return a consistent pixel size
    % frame-to-frame. Forcing Visible='off' (above) was NOT sufficient by
    % itself to prevent this on at least one real machine -- confirmed
    % 2026-07-13: getframe output was inconsistently scaled (a HiDPI/Retina
    % backing-store factor, e.g. exactly 1.5x, showed up on SOME frames but
    % not others of the SAME headless render), which crashes VideoWriter
    % ("Frame must be W by H") the moment a mismatched frame appears. Rather
    % than continuing to chase the display-dependent root cause blind (this
    % session cannot attach to the reporting machine's display to reproduce
    % it), every captured frame is unconditionally normalized to the PROMISED
    % [fig_h fig_w] size (via resize_frame, plain bilinear interpolation, no
    % toolbox required) before it reaches writeVideo -- so a scaling
    % inconsistency, whatever its cause, can no longer crash the render. A
    % one-time console note fires the first time a resize is actually needed.
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
                         '  normalizing every frame to %dx%d to match (this machine''s display is scaling\n', ...
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
    % Nearest multiple of 16. Confirmed empirically (isolate_16align test,
    % 2026-07-12): MATLAB's MPEG-4 writer on this machine silently corrupts
    % frames (diagonal streaking, wrong colours) at large sizes unless
    % width/height are multiples of 16 -- a classic H.264 macroblock-alignment
    % requirement, evidently enforced by a hardware-accelerated encoder path
    % that a plain "even" dimension (the previous, insufficient guard) does not
    % satisfy. Applied unconditionally, not just for large resolution_scale,
    % since it's not certain at what size threshold the corrupting code path
    % activates.
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
