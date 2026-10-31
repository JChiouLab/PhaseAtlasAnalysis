function report = plot_hss_dispersion_atlas(opts)
% plot_hss_dispersion_atlas  --  HSS linear-stability atlas by dispersion relation.
%
% Linearizes the 5-species RD reaction about the per-cell homogeneous steady
% state, forms the dispersion relation lambda(k_n) = max Re eig(J - k_n^2 D)
% over the discrete no-flux modes k_n = n*pi/ell (n>=1), and maps, across the
% (Dn_PAK, Dn_GEF) grid, WHERE the HSS is linearly unstable and the ONSET
% (fastest-growing) mode number n*.
%
% Plot opacity encodes the onset as SPATIAL PERIODS on the half-domain = n*/2
% (convention 2, user-chosen 2026-07-17): 0.5 period = one half-peak (HPSS
% morphology at x=0), 1.0 = one full period, 1.5, 2.0, ...  The mode wavelength
% is Lambda = L/n*; on the reflected FULL domain the same onset is n* full peaks.
%
% IMPORTANT (Invariant-9-adjacent): this is the LINEAR ONSET — the pattern that
% first emerges from the uniform state — NOT the final coarsened peak count.
% Mass-conserved systems coarsen; do not read n*/2 as the attractor's peak number.
%
% Only n>=1 modes are scanned: the n=0 (homogeneous) modes are the two exactly-
% marginal mass-conservation eigenvalues (Ru+Rq+Rv+Rw==0, Rp+Rq==0 give J two
% left null-vectors), so restricting to n>=1 both respects mass conservation of
% the pattern-forming perturbations and sidesteps the k=0 marginal modes.
%
% opts fields:
%   hss_path   (required) branch_HSS_latest.mat (has *_hss_branches, grids, par)
%   save_dir   (required) folder for the PDF + .mat
%   domain     'half' (default, ell=L/2) or 'full' (ell=L)
%   n_max      max mode index scanned (default 80; asserts the argmax is interior)
%   tol        growth-rate instability threshold (default 1e-6)
%   visible    'off' (default) | 'on'
%
% Output: HSS_Dispersion_Stability.pdf + hss_dispersion_atlas.mat in save_dir,
% and a `report` struct with the computed grids and summary stats.

    if ~isfield(opts,'hss_path') || ~isfield(opts,'save_dir')
        error('plot_hss_dispersion_atlas: opts.hss_path and opts.save_dir required.');
    end
    domain  = getf(opts,'domain','half');
    n_max   = getf(opts,'n_max',80);
    tol     = getf(opts,'tol',1e-6);
    visible = getf(opts,'visible','off');
    if ~exist(opts.save_dir,'dir'), mkdir(opts.save_dir); end

    d   = load(opts.hss_path);
    par = d.par; L = d.L_total; res = d.res;
    switch lower(domain)
        case 'half', ell = L/2;
        case 'full', ell = L;
        otherwise, error('opts.domain must be ''half'' or ''full''.');
    end
    Ddiag = diag([par.Dm par.Dm par.Dc par.Dc par.Dc]);
    kn = (1:n_max) * pi / ell;

    % The branch-library SOL cell is the authoritative representative HSS.
    % Legacy tracked arrays can omit a discovered unique root after branch-ID
    % pruning; build_branch_library recovers those roots into SOL with
    % representative_branch_index=-1. Prefer SOL so dispersion coverage exactly
    % follows existence_mask. Fall back to the old arrays for older libraries.
    has_sol = isfield(d,'SOL') && isequal(size(d.SOL),[res res]);
    if ~has_sol
        rbi = d.representative_branch_index;
        Ub=d.U_hss_branches; Qb=d.Q_hss_branches; Vb=d.V_hss_branches;
        Wb=d.W_hss_branches; Pb=d.P_hss_branches;
    end

    unstable = false(res); nstar = nan(res); maxg = nan(res); osc = false(res);
    argmax_hit_edge = false;
    for i = 1:res
        for j = 1:res
            if ~d.existence_mask(i,j), continue; end
            if has_sol
                uu = d.SOL{i,j}(:);
            else
                k = rbi(i,j);
                if isnan(k) || k < 1, continue; end
                uu = [Ub(i,j,k); Qb(i,j,k); Vb(i,j,k); Wb(i,j,k); Pb(i,j,k)];
            end
            if any(~isfinite(uu)) || uu(1) <= 0, continue; end
            J = hss_jacobian(uu, par);
            g = -inf; ns = 0; is_osc = false;
            for m = 1:n_max
                ev = eig(J - kn(m)^2 * Ddiag);
                [gm, ix] = max(real(ev));
                if gm > g, g = gm; ns = m; is_osc = abs(imag(ev(ix))) > 1e-9; end
            end
            maxg(i,j) = g;
            if g > tol
                unstable(i,j) = true; nstar(i,j) = ns; osc(i,j) = is_osc;
                if ns >= n_max, argmax_hit_edge = true; end
            end
        end
    end
    if argmax_hit_edge
        warning(['plot_hss_dispersion_atlas: some onset modes hit n_max=%d — ', ...
                 'the fastest-growing mode may be truncated; raise opts.n_max.'], n_max);
    end

    % ---- summary ----
    ncomp = sum(~isnan(maxg(:)));
    nun   = sum(unstable(:));
    report = struct('unstable',unstable,'nstar',nstar,'maxgrowth',maxg, ...
        'osc',osc,'PAK_DIFF',d.PAK_DIFF,'GEF_DIFF',d.GEF_DIFF,'domain',domain, ...
        'ell',ell,'n_max',n_max,'tol',tol,'par',par,'L_total',L,'res',res, ...
        'M_BG',d.M_BG,'M_P',d.M_P,'source_hss',opts.hss_path);
    report.frac_unstable = nun / max(1,ncomp);
    fprintf(['HSS dispersion atlas (%s-domain, ell=%.3g): %d/%d HSS cells unstable ', ...
             '(%.1f%%). onset n*: '], domain, ell, nun, ncomp, 100*report.frac_unstable);
    if nun>0
        ns_all = nstar(unstable);
        for v = 1:max(ns_all(:))
            cnt = sum(ns_all==v); if cnt>0, fprintf('n*=%d:%d ', v, cnt); end
        end
        fprintf('| osc %d / stat %d', sum(osc(:)), nun-sum(osc(:)));
    end
    fprintf('\n');

    % Convention 2 (user-chosen 2026-07-17): report the onset as the number of
    % SPATIAL PERIODS across the half-domain [0, L/2] = n*/2. So 0.5 period = one
    % half-peak (HPSS morphology), 1.0 = one full period, 1.5, 2.0, ...  (The
    % underlying integer mode index n* is kept too; periods = n*/2.)
    report.onset_periods = nstar / 2;
    save(fullfile(opts.save_dir,'hss_dispersion_atlas.mat'), '-struct', 'report');

    % ---- plot ----
    X = d.PAK_DIFF; Y = d.GEF_DIFF;
    nmax_star = max(1, max(nstar(unstable)));
    f = figure('Visible', visible, 'Color', 'w', 'Position', [80 80 1180 520]);

    % Panel A: onset spatial periods on the half-domain (= n*/2), opacity
    % light -> dark with more periods.
    axA = subplot(1,2,1); hold(axA,'on'); box(axA,'on');
    C = nstar;                              % colour by integer mode; NaN -> white
    hp = pcolor(axA, X, Y, C); set(hp,'EdgeColor','none');
    cmap = blues_opacity(nmax_star);
    colormap(axA, cmap); caxis(axA, [0.5, nmax_star+0.5]);
    cb = colorbar(axA); cb.Ticks = 1:nmax_star;
    cb.TickLabels = arrayfun(@(k) sprintf('%.1f', k/2), 1:nmax_star, 'UniformOutput', false);
    cb.Label.String = 'onset spatial periods on half-domain  (0.5 = 1 half-peak = HPSS)';
    contour(axA, X, Y, double(~isnan(nstar)), [0.5 0.5], 'k-', 'LineWidth', 0.6); % unstable boundary
    xlabel(axA,'\Deltan_{PAK}'); ylabel(axA,'\Deltan_{GEF}');
    title(axA, sprintf(['HSS linear instability — onset spatial periods (half-domain, ' ...
        '0.5 = 1 half-peak)\n%.1f%% of cells unstable'], 100*report.frac_unstable));
    axis(axA,'tight');

    % Panel B: onset character (stable / stationary Turing / oscillatory Turing)
    axB = subplot(1,2,2); hold(axB,'on'); box(axB,'on');
    cat = nan(res);                         % NaN = stable (white)
    cat(unstable & ~osc) = 1;               % stationary Turing
    cat(unstable &  osc) = 2;               % oscillatory (Hopf-Turing)
    hp2 = pcolor(axB, X, Y, cat); set(hp2,'EdgeColor','none');
    colormap(axB, [0.20 0.45 0.70; 0.85 0.30 0.25]);  % stat=blue, osc=red
    caxis(axB, [0.5 2.5]);
    cb2 = colorbar(axB); cb2.Ticks = [1 2]; cb2.TickLabels = {'stationary','oscillatory'};
    xlabel(axB,'\Deltan_{PAK}'); ylabel(axB,'\Deltan_{GEF}');
    title(axB, sprintf('Onset character (osc = Hopf-Turing \\rightarrow shuttle)\nosc %d / stat %d', ...
        sum(osc(:)), nun-sum(osc(:))));
    axis(axB,'tight');

    sgtitle(sprintf(['HSS dispersion-relation stability atlas   (a=%g b=%g c=%g d=%g e=%g f=%g, ', ...
        'Dm=%g Dc=%g, L=%g, M_{BG}=%g, M_P=%g)   —   onset periods are LINEAR onset, not the final coarsened count'], ...
        par.a,par.b,par.c,par.d,par.e,par.f,par.Dm,par.Dc,L,d.M_BG,d.M_P), 'FontSize', 9);

    out_pdf = fullfile(opts.save_dir,'HSS_Dispersion_Stability.pdf');
    exportgraphics(f, out_pdf, 'ContentType', 'vector');
    close(f);
    fprintf('  wrote %s\n', out_pdf);
    report.pdf = out_pdf;
end

function J = hss_jacobian(uu, par)
    % Analytic reaction Jacobian of pde_eq_v7's source term at state uu.
    % Verified to satisfy the two mass-conservation null-vectors [1 1 1 1 0]
    % and [0 1 0 0 1] exactly (2026-07-17).
    u=uu(1); v=uu(3); p=uu(5);
    a=par.a; b=par.b; c=par.c; dd=par.d; e=par.e; f=par.f;
    J = [ 2*a*u*v-b-c*p,   0,        a*u^2+f,     0,   -c*u;
          c*p,            -dd,       0,           0,    c*u;
         -2*a*u*v+b,       0,       -(a*u^2+f),   e,    0;
          0,               dd,       0,          -e,    0;
         -c*p,             dd,       0,           0,   -c*u];
end

function cmap = blues_opacity(nlev)
    % Discrete light->dark blue: reads as increasing opacity with onset mode.
    nlev = max(1, nlev);
    t = linspace(0.18, 0.92, nlev)';        % darkness
    base = [0.12 0.30 0.62];
    cmap = (1-t).*[1 1 1] + t.*base;        % white -> base blue
end

function v = getf(s, name, def)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end
