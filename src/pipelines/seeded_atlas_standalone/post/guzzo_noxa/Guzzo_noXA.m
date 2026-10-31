function Guzzo_noXA
    % GUZZO_NOXA
    % Four-field reaction-diffusion model (RomR / MglB).
    % Reduced further from Guzzo_noX.m: MglA (A) is also eliminated,
    % leaving only RomR and MglB. A appears in B's equation only through
    % dBA*Am*Bm^2. Its quasi-steady state
    %   Am ~ kA*rho_A*Rm / (dA + dAB*Bm^2 + kA*Rm)
    % shows A is proportional to Rm and inhibited by Bm, so it is replaced by:
    %   Am  ->  Rm / (1 + (Bm/KA)^2)
    % The inhibition by B is essential: a purely linear substitution (Am -> Rm)
    % gave only 2 weakly damped hits in 1700 random parameter sets, whereas the
    % saturating form above gave 31 sustained limit cycles in 1200 sets.
    %
    % Field order (pdepe's u(1..4)):
    %   u1 = Rm  (membrane RomR)   u2 = Rc  (cytosolic RomR)
    %   u3 = Bm  (membrane MglB)   u4 = Bc  (cytosolic MglB)
    %
    % Two conserved masses (zero-flux BC + reaction terms cancel pairwise):
    %   ∫(Rm+Rc)dx = N_R,  ∫(Bm+Bc)dx = N_B
    %
    % Membrane equations (each cytosolic partner is the negative, with Dc):
    %   Aeff   = Rm / (1 + (Bm/KA)^2)
    %   dRm/dt = Dm*Rm_xx + (kR + kRB*Bm)*Rc - dR*Rm
    %   dBm/dt = Dm*Bm_xx + (kB + kBB*Bm)*Bc - dB*Bm/(1 + Bm/K) - dBA*Aeff*Bm^2

    clc; clear;

    % ---------------------------------------------------------
    % 1. PARAMETERS
    % ---------------------------------------------------------
    % Parameters come from a random search (1200 sets) screened for left-right
    % symmetry:
    %   half-period mirror condition B(x, t+T/2) = B(L-x, t), error 0.0038
    %   peak ratio 1.000 (both poles reach 1.057), dwell times 49/51
    % The result is a symmetric pole-to-pole oscillation with a period of about
    % 10.7 min (comparable to the ~10 min of the original eight-field model);
    % Bm ranges over 0.225-1.355.
    %
    % (A higher-contrast set with a ~36 min period is also available.)

    L      = 6;       % cell length (um)
    T_end  = 200;     % simulated time (min); about 19 periods
    N_mesh = 128;     % spatial resolution
    N_time = 601;     % number of output times

    % --- Diffusivities (um^2/min) ---
    Dm = 0.2;         % membrane
    Dc = 120;         % cytosol
    p.D = [Dm; Dc; Dm; Dc];

    % --- Attachment rates (min^-1) ---
    p.kR  = 0.028624;           % basal RomR attachment
    p.kRB = 1.0423;             % RomR attachment promoted by membrane MglB
    p.kB  = 8.1255;             % basal MglB attachment
    p.kBB = 3.5216;             % cooperative MglB self-recruitment

    % --- Detachment rates (min^-1) ---
    p.dR  = 0.094962;           % RomR detachment
    p.dB  = 5.4381;             % maximum rate of saturating MglB detachment
    p.dBA = 210.01;             % MglB removed by Aeff (~Aeff*Bm^2)
    p.K   = 0.11440;            % saturation constant of MglB detachment
    p.KA  = 0.14714;            % half-saturation constant for the inhibition of Aeff by Bm

    % --- Total masses (global constraints) ---
    rho_R = 1;     N_R = rho_R * L;
    rho_B = 1;     N_B = rho_B * L;
    p.rho = [rho_R; rho_B];

    % --- Perturbation amplitude (relative to the smaller HSS value of each pair) ---
    amp = 0.1;

    % --- Movie settings ---
    vid.save     = true;
    vid.filename = 'Guzzo_noXA';
    vid.fps      = 20;
    vid.stride   = 6;                   % 601 frames -> ~100 frames (a 5 s movie)
    vid.dpi      = 75;                  % print resolution; lowering it speeds recording up a lot

    % ---------------------------------------------------------
    % 2. SOLVE FOR THE HOMOGENEOUS STEADY STATE (HSS)
    % ---------------------------------------------------------
    fprintf('Solving HSS for densities: R=%.3f, B=%.3f\n', p.rho);

    options = optimset('Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);

    guesses = [0.5 0.5; 0.9 0.9; 0.2 0.8; ...
               0.8 0.2; 0.05 0.5; 0.1 0.1];
    best = []; best_res = inf;
    for k = 1:size(guesses,1)
        g0 = guesses(k,:) .* p.rho(:)';
        [sol_k, fval_k, flag_k] = fsolve(@(y) hss_equations(y, p), g0, options);
        res = norm(fval_k);
        ok  = all(sol_k(:) > 0) && all(sol_k(:) < p.rho(:)) && res < 1e-6;
        if flag_k > 0 && ok && res < best_res
            best = sol_k; best_res = res;
        end
    end

    if isempty(best)
        error('No physical HSS found. Adjust kB, dBA, KA or the totals rho_*.');
    end

    Rm_ss = best(1); Bm_ss = best(2);
    hss = [Rm_ss; rho_R - Rm_ss; ...
           Bm_ss; rho_B - Bm_ss];

    fprintf('  -> HSS: Rm=%.4f Rc=%.4f | Bm=%.4f Bc=%.4f  (res=%.1e)\n', ...
            hss, best_res);

    % ---------------------------------------------------------
    % 3. RUN THE PDE
    % ---------------------------------------------------------
    m = 0;
    x = linspace(0, L, N_mesh);
    t = linspace(0, T_end, N_time);

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    fprintf('  -> starting pdepe (%d mesh x %d fields = %d ODEs)...\n', ...
            N_mesh, numel(p.D), N_mesh*numel(p.D));
    tPDE = tic;
    sol = pdepe(m, @(x,t,u,dudx) pdefun(x,t,u,dudx,p), ...
                   @(x) icfun(x, hss, L, amp), ...
                   @bcfun, x, t, opts);

    fprintf('  -> pdepe finished in %.1f s\n', toc(tPDE));

    assignin('base', 't', t);
    assignin('base', 'sol', sol);
    assignin('base', 'x', x);

    % --- If integration fails early, pdepe returns fewer times than requested;
    %     truncate t so plots and movie stay aligned ---
    if size(sol,1) < numel(t)
        warning(['Integration incomplete: %d/%d time points (last t=%.4g). ' ...
                 'Plots and movie cover only this range.'], size(sol,1), numel(t), t(size(sol,1)));
        t = t(1:size(sol,1));
    end

    U = squeeze(sol(end,:,:))';             % 8 x N_mesh

    % --- Mass check ---
    fprintf('  -> Mass R: %.4f (target %.4f)\n', trapz(x, U(1,:) + U(2,:)), N_R);
    fprintf('  -> Mass B: %.4f (target %.4f)\n', trapz(x, U(3,:) + U(4,:)), N_B);

    % ---------------------------------------------------------
    % 4. PLOTS
    % ---------------------------------------------------------
    names  = {'Rm','Rc','Bm','Bc'};
    pairs  = {[1 2], [3 4]};
    titles = {'RomR','MglB'};

    figure('Name', sprintf('Guzzo noXA (R,B only), L=%g', L), 'Color', 'w', ...
           'Position', [100 100 900 650]);

    % Top row: final spatial profiles (membrane solid / cytosolic dashed)
    for k = 1:2
        subplot(2,2,k);
        idx = pairs{k};
        plot(x, U(idx(1),:), 'LineWidth', 2); hold on;
        plot(x, U(idx(2),:), '--', 'LineWidth', 1.5);
        xlabel('x (\mum)'); ylabel('Concentration');
        title([titles{k} sprintf('  (t=%g min)', t(end))]);
        legend(names{idx(1)}, names{idx(2)}, 'Location', 'best');
        grid on; xlim([0 L]);
    end

    % Bottom row: kymographs of the membrane fields
    memb = [1 3];
    for k = 1:2
        subplot(2,2,2+k);
        imagesc(t, x, sol(1:numel(t),:,memb(k))');
        set(gca, 'YDir', 'normal'); colormap('jet'); colorbar;
        xlim([t(1) t(end)]);
        xlabel('Time (min)'); ylabel('x (\mum)');
        title([names{memb(k)} '(x,t)']);
    end

    % ---------------------------------------------------------
    % 5. WRITE THE MOVIE
    % ---------------------------------------------------------
    if vid.save
        tVID = tic;
        write_movie(sol, x, t, vid, pairs, titles, names);
        fprintf('  -> movie finished in %.1f s\n', toc(tVID));
    end
end

% =========================================================
% A. HSS: reaction terms = 0, subject to mass conservation
%    Unknowns y = [Rm, Bm]; cytosolic values follow from the totals
% =========================================================
function F = hss_equations(y, p)
    Rm = y(1); Bm = y(2);
    Rc = p.rho(1) - Rm;
    Bc = p.rho(2) - Bm;

    [sR, sB] = reactions(Rm, Rc, Bm, Bc, p);
    F = [sR; sB];
end

% =========================================================
% B. Initial condition: HSS + zero-mean perturbation (mass exactly conserved)
% =========================================================
function u0 = icfun(x, hss, L, amp)
    pert = cos(pi * x / L);
    u0 = hss(:);
    % Scale by the smaller of the membrane/cytosolic pair so neither goes negative
    for k = 1:2
        m = 2*k-1; c = 2*k;
        d = amp * min(hss(m), hss(c)) * pert;
        u0(m) = hss(m) + d;
        u0(c) = hss(c) - d;
    end
end

% =========================================================
% C. Reaction terms (shared by the HSS solver and the PDE, so they cannot drift apart)
% =========================================================
function [sR, sB] = reactions(Rm, Rc, Bm, Bc, p)
    % RomR
    onR = (p.kR + p.kRB * Bm) * Rc;
    ofR = p.dR * Rm;

    % Stand-in for the eliminated MglA: proportional to Rm, inhibited by Bm
    Aeff = Rm / (1 + (Bm / p.KA)^2);

    % MglB: saturating detachment + Aeff-dependent removal
    onB = (p.kB + p.kBB * Bm) * Bc;
    ofB = p.dB * Bm / (1 + Bm / p.K) + p.dBA * Aeff * Bm^2;

    sR = onR - ofR;
    sB = onB - ofB;
end

% =========================================================
% D. PDE equations
% =========================================================
function [c, f, s] = pdefun(~, ~, u, dudx, p)
    c = ones(4, 1);
    f = p.D .* dudx;

    [sR, sB] = reactions(u(1), u(2), u(3), u(4), p);

    s = [ sR; -sR; ...
          sB; -sB ];
end

% =========================================================
% E. Boundary conditions: zero flux for every field
% =========================================================
function [pl, ql, pr, qr] = bcfun(~, ~, ~, ~, ~)
    pl = zeros(4,1); ql = ones(4,1);
    pr = zeros(4,1); qr = ones(4,1);
end

% =========================================================
% Write the time evolution to an mp4
% On error the partial file is deleted, so no broken movie is left behind
% =========================================================
function write_movie(sol, x, t, vid, pairs, titles, names)
    outfile = [vid.filename '.mp4'];

    vw = VideoWriter(vid.filename, 'MPEG-4');
    vw.FrameRate = vid.fps;
    vw.Quality   = 95;
    open(vw);

    % Capture with print, not getframe: getframe reads screen pixels and returns blank
    % frames when the figure is hidden or covered.
    fig = figure('Color', 'w', 'Units', 'pixels', ...
                 'Position', [100 100 960 352], ...
                 'Visible', 'off', 'PaperPositionMode', 'auto');

    ok = false; errMsg = '';
    try
        nT     = min(size(sol,1), numel(t));
        frames = 1:vid.stride:nT;

        % Fix the y range so frames do not rescale and jump around
        ylims = zeros(2, 2);
        for k = 1:2
            vals = sol(1:nT,:,pairs{k});
            lo = min(vals(:)); hi = max(vals(:));
            pad = 0.05 * max(hi - lo, eps);
            ylims(k,:) = [lo - pad, hi + pad];
        end

        % Build the lines once; each frame only updates YData
        hLine = gobjects(2,2);
        U0 = squeeze(sol(frames(1),:,:))';
        for k = 1:2
            ax  = subplot(1, 2, k, 'Parent', fig);
            idx = pairs{k};
            hLine(k,1) = plot(ax, x, U0(idx(1),:), 'LineWidth', 2); hold(ax, 'on');
            hLine(k,2) = plot(ax, x, U0(idx(2),:), '--', 'LineWidth', 1.5);
            xlabel(ax, 'x (\mum)'); ylabel(ax, 'Concentration');
            title(ax, titles{k});
            legend(ax, names{idx(1)}, names{idx(2)}, 'Location', 'best');
            grid(ax, 'on');
            xlim(ax, [x(1) x(end)]); ylim(ax, ylims(k,:));
        end
        hSup = sgtitle(fig, 't = 0 min', 'FontWeight', 'bold');

        sz = [];
        for n = frames
            U = squeeze(sol(n,:,:))';
            for k = 1:2
                idx = pairs{k};
                set(hLine(k,1), 'YData', U(idx(1),:));
                set(hLine(k,2), 'YData', U(idx(2),:));
            end
            set(hSup, 'String', sprintf('t = %.1f min', t(n)));
            drawnow;

            img = print(fig, '-RGBImage', sprintf('-r%d', vid.dpi));

            % Crop to a multiple of 16: H.264 macroblocks are 16x16, and sizes that
            % do not divide evenly corrupt the edge blocks (black diagonal bands).
            % Only the outer white margin is cropped away.
            if isempty(sz)
                sz = [floor(size(img,1)/16)*16, floor(size(img,2)/16)*16];
            end
            img = img(1:sz(1), 1:sz(2), :);

            writeVideo(vw, img);
        end
        ok = true;
    catch ME
        errMsg = ME.message;
    end

    close(vw);
    if ishandle(fig), close(fig); end

    if ok
        fprintf('  -> movie saved: %s (%d frames, %g fps)\n', ...
                outfile, numel(frames), vid.fps);
    else
        if exist(outfile, 'file'), delete(outfile); end
        warning('Recording failed; deleted the incomplete file %s\n%s', outfile, errMsg);
    end
end
