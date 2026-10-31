function report = run_atlas_pde_trace(opts)
    % RUN_ATLAS_PDE_TRACE  Run a full PDE simulation whose INITIAL CONDITION is an
    % exact atlas-sheet stitched steady state, so the trajectory is directly
    % comparable to the atlas by classify_pde_to_atlas_sheets.m (Packet 5).
    %
    % Analogous to get_pde_trace.m, but (a) the IC is built from the atlas branch
    % libraries via the SHARED reconstruct_stitched_sheet.m -- so at t=0 the PDE
    % frame equals the classifier's reconstruction of the same sheet cell to
    % machine precision -- and (b) it can stitch ANY of the nine sheet types
    % (HSS/HPSS/MPSS on either half), not just the old HPSS-only bvp file.
    %
    % The PDE model, mesh convention, no-flux BCs, and mass-asymmetry definition
    % are identical to get_pde_trace.m. It does NOT modify get_pde_trace.m and
    % writes to its own namespaced files (atlas_pde_trace_*.mat), never touching
    % pde_trace_latest.mat.
    %
    % Options (struct, all optional):
    %   sheet_name  - which sheet to build the IC from. Default 'HPSS_HPSS'.
    %   cell        - [i j] grid indices for the IC. If omitted, use (del_pak,
    %                 del_gef); if those omitted too, auto-pick the existing cell
    %                 of this sheet nearest the grid centre.
    %   del_pak,del_gef - target mass-asymmetry coords (snapped to nearest cell).
    %   T_end       - end time. Default 2000.
    %   n_frames    - number of saved time frames (BOUNDS memory: sol_pde is
    %                 n_frames x N_mesh x 5). Default 100.
    %   frame_spacing - 'linear' (default) or 'log' (t=0 plus log-spaced points,
    %                 dense at early times -- good for catching fast transients).
    %   t_span      - explicit strictly-increasing time vector starting at 0
    %                 (overrides T_end/n_frames/frame_spacing).
    %   N_mesh      - spatial mesh points over [0, L_total]. Default 200.
    %   atlas_path  - default data/outputs/atlas_sheets_latest.mat
    %   hss_path/hpss_path/mpss_path - branch files (defaults to *_latest.mat)
    %   out_dir     - default data/outputs
    %   save_plot   - write a summary PNG. Default true.
    %   visible     - 'on'|'off' figure visibility. Default 'off'.
    %   run_config  - opaque struct saved VERBATIM into the output .mat as the
    %                 field `run_config`. This function never interprets it; it
    %                 is a staleness key a driver can compare against on the next
    %                 run to decide whether the trace needs regenerating (see
    %                 run_atlas_surface_pipeline.m Stage 8). Default struct().
    %
    % Output: data/outputs/atlas_pde_trace_<id>.mat (+ atlas_pde_trace_latest.mat).
    % Saved fields: sol_pde (nt x nx x 5), x_full, t_span, nPAK_trace, nGEF_trace
    % (per-frame mass asymmetry), u0_stitched, par, L_total, M_BG, M_P, and IC/
    % provenance (ic_sheet, ic_cell, ic_mirror_cell, atlas_id, branch ids,
    % N_mesh, n_frames, frame_spacing, computation_id, run_config).

    if nargin < 1, opts = struct(); end
    opts = def(opts, 'sheet_name', 'HPSS_HPSS');
    opts = def(opts, 'cell', []);
    opts = def(opts, 'del_pak', []);
    opts = def(opts, 'del_gef', []);
    opts = def(opts, 'T_end', 2000);
    opts = def(opts, 'n_frames', 100);
    opts = def(opts, 'frame_spacing', 'linear');
    opts = def(opts, 't_span', []);
    opts = def(opts, 'N_mesh', 200);
    opts = def(opts, 'atlas_path', 'data/outputs/atlas_sheets_latest.mat');
    opts = def(opts, 'hss_path',  'data/outputs/branch_HSS_latest.mat');
    opts = def(opts, 'hpss_path', 'data/outputs/branch_HPSS_latest.mat');
    opts = def(opts, 'mpss_path', 'data/outputs/branch_MPSS_latest.mat');
    opts = def(opts, 'out_dir', 'data/outputs');
    opts = def(opts, 'save_plot', true);
    opts = def(opts, 'visible', 'off');
    opts = def(opts, 'run_config', struct());   % opaque caller-supplied staleness key,
                                                 % saved verbatim into the output .mat so a
                                                 % driver can detect when its inputs changed
                                                 % (see run_atlas_surface_pipeline.m Stage 8).

    computation_id = datestr(now, 'yyyymmdd_HHMMSS');

    % ---- load atlas + branches (read-only) ----
    A = load(opts.atlas_path);
    branches.HSS  = load(opts.hss_path);
    branches.HPSS = load(opts.hpss_path);
    branches.MPSS = load(opts.mpss_path);
    if ~any(strcmp(opts.sheet_name, A.sheet_names))
        error('run_atlas_pde_trace: unknown sheet_name "%s".', opts.sheet_name);
    end
    par = A.par;
    if ~isfield(par, 'D'), par.D = [par.Dm; par.Dm; par.Dc; par.Dc; par.Dc]; end
    L_total = A.L_total; L_half = L_total / 2;
    M_BG = A.M_BG; M_P = A.M_P; res = A.res;
    sh = A.sheets.(opts.sheet_name);

    % ---- choose IC cell ----
    if ~isempty(opts.cell)
        i0 = opts.cell(1); j0 = opts.cell(2);
    elseif ~isempty(opts.del_pak) && ~isempty(opts.del_gef)
        d = (A.PAK_DIFF - opts.del_pak).^2 + (A.GEF_DIFF - opts.del_gef).^2;
        [~, k] = min(d(:)); [i0, j0] = ind2sub(size(d), k);
    else
        i0 = round(res/2); j0 = round(res/2);
    end
    if ~sh.sheet_exists(i0, j0)
        % auto-pick the nearest existing cell for this sheet
        [ii, jj] = find(sh.sheet_exists);
        if isempty(ii)
            error('run_atlas_pde_trace: sheet "%s" has no existing cells.', opts.sheet_name);
        end
        dd = (ii - i0).^2 + (jj - j0).^2;
        [~, k] = min(dd); i1 = ii(k); j1 = jj(k);
        fprintf('run_atlas_pde_trace: sheet %s does not exist at (%d,%d); using nearest existing cell (%d,%d).\n', ...
            opts.sheet_name, i0, j0, i1, j1);
        i0 = i1; j0 = j1;
    end
    mi = res + 1 - i0; mj = res + 1 - j0;

    fprintf('=== run_atlas_pde_trace ===\n');
    fprintf('  IC sheet=%s  cell=(%d,%d)  mirror=(%d,%d)  dPAK=%.3f dGEF=%.3f\n', ...
        opts.sheet_name, i0, j0, mi, mj, A.PAK_DIFF(i0,j0), A.GEF_DIFF(i0,j0));

    % ---- build IC via the shared reconstruction (bit-identical to classifier) ----
    x_full = linspace(0, L_total, opts.N_mesh);
    u0_stitched = reconstruct_stitched_sheet(A, branches, opts.sheet_name, i0, j0, x_full);
    if isempty(u0_stitched)
        error('run_atlas_pde_trace: could not reconstruct IC for sheet %s at (%d,%d).', ...
            opts.sheet_name, i0, j0);
    end

    % ---- time span (bounded frames) ----
    if ~isempty(opts.t_span)
        t_span = opts.t_span(:).';
        if t_span(1) ~= 0 || any(diff(t_span) <= 0)
            error('run_atlas_pde_trace: t_span must start at 0 and be strictly increasing.');
        end
    elseif strcmpi(opts.frame_spacing, 'log')
        t0 = opts.T_end / 1e4;
        t_span = [0, logspace(log10(t0), log10(opts.T_end), opts.n_frames - 1)];
    else
        t_span = linspace(0, opts.T_end, opts.n_frames);
    end
    nt = numel(t_span);
    n_elems = nt * opts.N_mesh * 5;
    if n_elems > 5e6
        warning('run_atlas_pde_trace: sol_pde has %d elements (%.0f MB) -- consider fewer frames or mesh points.', ...
            n_elems, n_elems*8/1e6);
    end

    % ---- run PDE ----
    fprintf('  running pdepe: N_mesh=%d, n_frames=%d (%s), T_end=%.1f ...\n', ...
        opts.N_mesh, nt, opts.frame_spacing, t_span(end));
    sol_pde = pdepe(0, @(x,t,u,dudx) pdefun_local(x,t,u,dudx,par), ...
                       @(xi) interp1(x_full, u0_stitched', xi)', ...
                       @bcfun_local, x_full, t_span);

    % ---- per-frame mass asymmetry (same definition as get_pde_trace.m) ----
    idx_L = x_full <= L_half; idx_R = ~idx_L;
    xL = x_full(idx_L); xR = x_full(idx_R);
    nPAK_trace = zeros(nt, 1); nGEF_trace = zeros(nt, 1);
    for f = 1:nt
        yf = squeeze(sol_pde(f, :, :))';                     % 5 x nx
        gefL = trapz(xL, sum(yf(1:4, idx_L), 1)); gefR = trapz(xR, sum(yf(1:4, idx_R), 1));
        pakL = trapz(xL, yf(2, idx_L) + yf(5, idx_L)); pakR = trapz(xR, yf(2, idx_R) + yf(5, idx_R));
        nGEF_trace(f) = gefL - gefR;
        nPAK_trace(f) = pakL - pakR;
    end

    % ---- provenance ----
    ic_sheet = opts.sheet_name; ic_cell = [i0 j0]; ic_mirror_cell = [mi mj];
    atlas_id = A.computation_id;
    branch_ids = struct('HSS', branches.HSS.computation_id, ...
                        'HPSS', branches.HPSS.computation_id, 'MPSS', branches.MPSS.computation_id);
    N_mesh = opts.N_mesh; n_frames = nt; frame_spacing = opts.frame_spacing; %#ok<NASGU>
    run_config = opts.run_config; %#ok<NASGU>   % staleness key (empty unless a driver set it)

    if ~exist(opts.out_dir, 'dir'), mkdir(opts.out_dir); end
    out_path = fullfile(opts.out_dir, sprintf('atlas_pde_trace_%s.mat', computation_id));
    save(out_path, 'sol_pde', 'x_full', 't_span', 'nPAK_trace', 'nGEF_trace', ...
        'u0_stitched', 'par', 'L_total', 'M_BG', 'M_P', ...
        'ic_sheet', 'ic_cell', 'ic_mirror_cell', 'atlas_id', 'branch_ids', ...
        'N_mesh', 'n_frames', 'frame_spacing', 'computation_id', 'run_config');
    latest_path = fullfile(opts.out_dir, 'atlas_pde_trace_latest.mat');
    copyfile(out_path, latest_path);
    fprintf('  wrote %s\n  wrote %s\n', out_path, latest_path);

    if opts.save_plot
        summary_plot(x_full, t_span, sol_pde, u0_stitched, nPAK_trace, nGEF_trace, ...
                     opts.sheet_name, computation_id, opts.visible, opts.out_dir);
    end

    report = struct('out_path', out_path, 'computation_id', computation_id, ...
        'ic_sheet', ic_sheet, 'ic_cell', ic_cell, 'nt', nt, 'N_mesh', opts.N_mesh, ...
        'sol_pde_MB', n_elems*8/1e6);
    fprintf('=== run_atlas_pde_trace complete (sol_pde %.1f MB) ===\n', report.sol_pde_MB);
end


% ---- PDE model (identical to get_pde_trace.m) ----
function [c, f, s] = pdefun_local(~, ~, u, dudx, p)
    c = ones(5,1);
    if isfield(p, 'D'), f = p.D .* dudx; else, f = [p.Dm; p.Dm; p.Dc; p.Dc; p.Dc] .* dudx; end
    U=u(1); Q=u(2); V=u(3); W=u(4); P=u(5);
    Ru = (p.a*U^2+p.f)*V - p.b*U - p.c*U*P; Rq = p.c*U*P - p.d*Q;
    Rv = -(p.a*U^2+p.f)*V + p.b*U + p.e*W; Rw = p.d*Q - p.e*W; Rp = -p.c*U*P + p.d*Q;
    s = [Ru; Rq; Rv; Rw; Rp];
end
function [pl, ql, pr, qr] = bcfun_local(~,~,~,~,~)
    pl=zeros(5,1); ql=ones(5,1); pr=zeros(5,1); qr=ones(5,1);
end

function opts = def(opts, name, val)
    if ~isfield(opts, name), opts.(name) = val; end
end

function summary_plot(x_full, t_span, sol_pde, u0, nPAK, nGEF, sheet_name, id, vis, out_dir)
    fig = figure('Color','w','Position',[100 100 1200 380],'Visible',vis);
    subplot(1,3,1);
    plot(x_full, u0(1,:), 'k--', 'LineWidth', 1); hold on;
    plot(x_full, squeeze(sol_pde(end,:,1)), 'r-', 'LineWidth', 2);
    title(sprintf('u(x): IC (%s) vs final', sheet_name), 'Interpreter','none');
    xlabel('x'); ylabel('u'); legend('IC','final PDE'); grid on;
    subplot(1,3,2);
    imagesc(t_span, x_full, squeeze(sol_pde(:,:,1))'); set(gca,'YDir','normal');
    colormap('jet'); colorbar; title('u(x,t)'); xlabel('t'); ylabel('x');
    subplot(1,3,3);
    plot(nPAK, nGEF, 'b-', 'LineWidth', 2); hold on;
    plot(nPAK(1), nGEF(1), 'go', 'MarkerFaceColor','g');
    xlabel('\Delta n_{PAK}'); ylabel('\Delta n_{GEF}'); title('mass asymmetry'); grid on;
    % Diagnostic PNG: goes into the trace's own out_dir (the caller's IC-specific
    % folder), not a hardcoded flat data/plots. (Still wrapped in try/catch --
    % this install's PNG raster export is separately broken; see the movie tool.)
    if ~exist(out_dir,'dir'), mkdir(out_dir); end
    try, exportgraphics(fig, fullfile(out_dir, sprintf('Atlas_PDE_Trace_%s.png', id))); catch, end
    if strcmp(vis,'off'), close(fig); end
end
