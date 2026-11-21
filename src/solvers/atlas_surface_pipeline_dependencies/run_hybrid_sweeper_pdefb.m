function run_hybrid_sweeper_pdefb(p)
    % RUN_HYBRID_SWEEPER_PDEFB  Mass-axis sweep of the half-peak (HPSS) half-cell state.
    % -----------------------------------------------------------------
    % Goal: at fixed L (= p.L), build a (Delta n_PAK, Delta n_GEF) grid of
    % half-peak steady states of an isolated half cell (no flux at the ends; fission yeast model),
    % using the 14-variable BVP of the shared per-point core (5 species + 5 derivatives +
    % 2 uniform sources + 2 mass accumulators). Each solved point is used as the starting
    % profile for its neighbours (continuation).
    %
    % Note on the paper: these half-cell states are the starting profiles of the paper's coupled
    % pipeline; the rest of the uncoupled workflow built on them is an earlier approximation.
    %
    % Phases:
    %   Phase 0: seeding at the centre's mass density (or at another anchor bin if the centre
    %            has no half-peak state, see below).
    %            * Interactive (the default when p.seeding_mode is not set): run the PDE at
    %              several test L values, show a panel of profiles, ask the user to pick the
    %              cleanest one and how many half-waves it contains; extract the fundamental
    %              half-wave; creep UP or DOWN in L to L_half.
    %            * Auto: try a list of L candidates, accept the first with amplitude >= 0.5 in u
    %              (the no-peak rejection threshold).
    %            * Constructed guess (p.phase0_construct): solve the BVP directly from analytic
    %              edge-peak guesses, without PDE relaxation.
    %            All seeds and creep steps are checked by validate_halfwave, so a flat
    %            (no-peak) solution can never pass as a half-peak.
    %   Phase 1: Solve the anchor grid point with the seed.
    %   Phase 2: Column continuation from the anchor, two tiers at each step:
    %            (1) predictor-corrector BVP from the previous good column points,
    %            (2) PDE relaxation from an extrapolation of the last 1-2 good points,
    %                then a BVP at the relaxed state.
    %   Phase 3: Parallel row sweeps with the same predictor-corrector in PAK-mass space.
    %   Phase 4: Retry of every (i,j) the main sweep could not solve, from solved neighbours
    %            (parallel, multi-pass, per-cell time limit).
    %   Phase 5: Eigenvalue-based stability classification of every solved grid point
    %            (skipped if p.skip_stability).
    %
    % Failed points are not filled in. They are recorded in failure_mask and
    % left as NaN in the U/Q/V/W/P_trough arrays so downstream plots can mask them.
    %
    % Notes:
    %   * The BVP does not rescale the mass integrals; the sources s1, s2 absorb the residual.
    %   * Solutions are checked for shape (validate_halfwave) at every seed and accept point.

    computation_id = datestr(now, 'yyyymmdd_HHMMSS');

    % --- 1. Parameters (built-in default parameters, not the paper's set; see run_atlas_surface_pipeline.m / docs/MODELS.md) -------------------
    if nargin < 1
        p.a = 0.59; p.b = 0.13; p.c = 21.0;
        p.d = 1.3;  p.e = 33.0; p.f = 1.9;
        p.Dm = 0.001; p.Dc = 1.0;
        p.L = 4.0;
        p.M_BG = 20.0;
        p.M_P  = 8.0;
        p.res = 25;
    end
    % L-creep + PDE settings (override-able)
    if ~isfield(p, 'L_creep_step'),     p.L_creep_step    = 0.05; end
    if ~isfield(p, 'L_creep_start'),    p.L_creep_start   = 0.5;  end
    if ~isfield(p, 'pde_T_end'),        p.pde_T_end       = 2e7;  end
    if ~isfield(p, 'pde_T_chunk'),      p.pde_T_chunk     = 1e7;  end
    if ~isfield(p, 'pde_max_chunks'),   p.pde_max_chunks  = 20;   end
    if ~isfield(p, 'pde_fallback_T'),   p.pde_fallback_T  = 5e6;  end
    if ~isfield(p, 'skip_stability'),   p.skip_stability  = false; end
    % Oscillatory-regime PDE fallback (default OFF). When ON, the PDE-relaxation fallback (a) stops early once it
    % detects the trajectory is not relaxing to a steady state instead of
    % burning the whole chunk budget, and (b) offers several oscillation-phase
    % snapshots + their time-mean as bvp4c SEED candidates, accepting the first
    % that converges to a validated steady half-peak. Snapshots are seeds only;
    % morphology/mass acceptance is unchanged.
    if ~isfield(p, 'osc_fallback'),     p.osc_fallback    = false; end
    if ~isfield(p, 'osc_conv_tol'),     p.osc_conv_tol    = 1e-6;  end  % steady if inter-chunk max|du| < this
    if ~isfield(p, 'osc_min_chunks'),   p.osc_min_chunks  = 3;     end  % chunks before "not relaxing" may be declared
    if ~isfield(p, 'osc_rel'),          p.osc_rel         = 0.5;   end  % not relaxing if recent min|du| > osc_rel*first
    if ~isfield(p, 'osc_n_snap'),       p.osc_n_snap      = 6;     end  % phase snapshots per chunk
    % Phase 4 per-cell wall-clock budget (seconds). A single pathological
    % neighbour-seeded BVP call can otherwise stall bvp4c's mesh refinement
    % indefinitely and silently block the whole pass (so each cell runs as its own
    % parfeval task with a time limit).
    if ~isfield(p, 'phase4_cell_timeout'), p.phase4_cell_timeout = 30; end
    % Phase 0 seeding strategy:
    %   'interactive' (default) - PDE at several test L, user picks
    %                              the cleanest pattern, hybrid creeps to L_target.
    %   'auto'                  - small-L start, automatic creep up.
    if ~isfield(p, 'seeding_mode'),     p.seeding_mode    = 'interactive'; end
    if ~isfield(p, 'L_seed_test_values')
        % Default: 5 test L's bracketing L_half so the user can see how the
        % polarized branch develops with L.  Override if your bifurcation
        % point is far from these defaults.
        p.L_seed_test_values = (p.L/4) * [0.5, 0.75, 1.0, 1.5, 2.0];
    end
    % Phase 2 robustness (column-seeding fallback chain).
    %   subdivide_max_levels - bisection levels for the mass sub-step retry
    %                          (level k uses 2^k intermediate mass targets).
    %   centre_anchor_factor - multiplier on |i - mid_idx| for the number of
    %                          mass steps used by centre-anchored continuation
    %                          (final step size = grid_step / factor).
    if ~isfield(p, 'subdivide_max_levels'), p.subdivide_max_levels = 4; end
    if ~isfield(p, 'centre_anchor_factor'), p.centre_anchor_factor = 2; end

    par     = p;
    D_vals  = [p.Dm; p.Dm; p.Dc; p.Dc; p.Dc];
    L_half  = p.L / 2;
    L_total = p.L; M_BG = p.M_BG; M_P = p.M_P;   % needed early for checkpointing

    res        = p.res;
    n_GEF_bar  = p.M_BG / 2;
    n_PAK_bar  = p.M_P  / 2;

    deln_pak = linspace(-0.95, 0.95, res) * p.M_P;
    deln_G   = linspace(-0.95, 0.95, res) * p.M_BG;
    [PAK_DIFF, GEF_DIFF] = meshgrid(deln_pak, deln_G);

    % Output containers
    SOL_grid       = cell(res, res);
    failure_mask   = true(res, res);     % true == not solved (mark the cell, never fill it in)
    is_stable_grid = false(res, res);
    morphology_class_grid = repmat("NOT_FOUND", res, res);
    morphology_reason_grid = repmat("no_accepted_solution", res, res);
    morphology_diagnostics_grid = cell(res, res);
    quarantined_nonhpss_mask = false(res, res);
    U_trough = nan(res, res); Q_trough = nan(res, res);
    V_trough = nan(res, res); W_trough = nan(res, res); P_trough = nan(res, res);

    if isfield(p, 'out_dir') && ~isempty(p.out_dir), out_dir = p.out_dir; else, out_dir = 'data/outputs'; end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    filename = sprintf('Hybrid_Sweep_V7_%s.mat', computation_id);
    out_path = fullfile(out_dir, filename);

    % ====================================================================
    % PHASE 0/1  -  MULTI-BIN ANCHOR SEARCH
    % ====================================================================
    % Phase 0/1 does not fail outright when the centre bin has no seedable half-peak state,
    % e.g. when the centre mass is past a fold of the half-peak branch, so the static
    % half-peak does not exist there. Instead it tries up to p.phase01_max_bins candidate anchor
    % bins (default 9: a centre-first coarse 3x3 spread across the (Dn_PAK,
    % Dn_GEF) grid), running the Phase-0 seed + Phase-1 solve at each bin's own
    % mass density until one yields a validated half-peak. That bin becomes the anchor
    % (anchor_i, anchor_j) the column/row continuation marches out from; it need
    % not be the centre. Only if all candidates fail does the sweep error out.
    % (The quarter-domain sweeper re-roots the same way.)
    mid_idx = round(res/2);
    if ~isfield(p, 'phase01_max_bins'), p.phase01_max_bins = 9; end
    anchor_bins = candidate_anchor_bins(res, mid_idx, p.phase01_max_bins);
    fprintf('PHASE 0/1: anchor search over up to %d bin(s) (centre first)...\n', ...
            size(anchor_bins, 1));

    anchor_i = NaN; anchor_j = NaN; sol_anchor = []; seed_info = struct('n_steps', 0);
    for cb = 1:size(anchor_bins, 1)
        i0 = anchor_bins(cb, 1); j0 = anchor_bins(cb, 2);
        tg0 = n_GEF_bar + GEF_DIFF(i0, j0)/2;
        tp0 = n_PAK_bar + PAK_DIFF(i0, j0)/2;
        rho_G0 = tg0 / L_half; rho_P0 = tp0 / L_half;   % this bin's mass density
        fprintf('  [bin %d/%d] (i=%d,j=%d) tg=%.3f tp=%.3f: ', ...
                cb, size(anchor_bins, 1), i0, j0, tg0, tp0);

        [seed_sol, s_info] = phase0_lcreep_seed(par, D_vals, L_half, rho_G0, rho_P0, p);
        if isempty(seed_sol), fprintf('no Phase-0 seed.\n'); continue; end

        [sol0, ok] = run_bvp_v7(construct_bvp_guess_v7(seed_sol), L_half, par, D_vals, tg0, tp0);
        if ~ok
            sol0 = lcreep_fallback(L_half, par, D_vals, tg0, tp0, p);
            ok   = ~isempty(sol0);
        end
        if ok && validate_halfwave(sol0)
            anchor_i = i0; anchor_j = j0; sol_anchor = force_peak_at_zero(sol0); seed_info = s_info;
            fprintf('ANCHOR (u_max=%.4f, %d creep steps).\n', ...
                    max(sol0.y(1,:)), seed_info.n_steps);
            break;
        end
        fprintf('not seedable.\n');
    end
    if isempty(sol_anchor)
        error(['Phase 0/1 failed: none of the %d candidate anchor bins yielded a ', ...
               'validated HPSS at L_half=%.3f (the static half-peak may not exist ', ...
               'anywhere on this grid at these params).'], size(anchor_bins, 1), L_half);
    end
    if anchor_i ~= mid_idx || anchor_j ~= mid_idx
        fprintf('  NOTE: centre bin (%d,%d) was not seedable; anchored at (%d,%d) instead.\n', ...
                mid_idx, mid_idx, anchor_i, anchor_j);
    end
    SOL_grid{anchor_i, anchor_j}     = sol_anchor;
    failure_mask(anchor_i, anchor_j) = false;
    [U_trough(anchor_i,anchor_j), Q_trough(anchor_i,anchor_j), V_trough(anchor_i,anchor_j), ...
     W_trough(anchor_i,anchor_j), P_trough(anchor_i,anchor_j)] = extract_trough(sol_anchor);

    % ====================================================================
    % PHASE 2  -  COLUMN SEEDING (predictor-corrector + PDE fallback)
    % ====================================================================
    fprintf('PHASE 2 [pdefb]: Column seeding (anchor column j=%d from row i=%d), predictor + PDE fallback...\n', ...
            anchor_j, anchor_i);
    t_phase2 = tic;
    n_tier1 = 0; n_tier2 = 0;
    column_seeds = cell(res, 1);
    column_succ  = false(res, 1);
    column_seeds{anchor_i} = sol_anchor;
    column_succ(anchor_i)  = true;

    for dir_step = [1, -1]
        if dir_step > 0, i_range = anchor_i+1 : 1   : res;
        else,            i_range = anchor_i-1 : -1  : 1;
        end
        prev_sol      = sol_anchor;
        prev_prev_sol = [];
        for i = i_range
            tg = n_GEF_bar + GEF_DIFF(i, anchor_j)/2;
            tp = n_PAK_bar + PAK_DIFF(i, anchor_j)/2;

            % Tier 1: predictor-corrector BVP
            guess = predictor_step(prev_sol, prev_prev_sol);
            [bvp_sol, ok] = run_bvp_v7(guess, L_half, par, D_vals, tg, tp);
    % bvp4c routinely converges from a half-peak guess to a flat
            % (no-peak) profile. Accepting it would store a non-peak and also use it as the
            % starting profile of the next cell, propagating the collapse down the march
            % (at res = 50 about 1000/2500 cells were collapsed, and about 60% of the
            % cells next to a collapse recover a real half-peak from a validated
            % warm start). Reject a collapsed solve so tier 2 (PDE re-polarisation) runs.
            ok = ok && validate_halfwave(bvp_sol);
            tier_used = 1;

            % Tier 2: PDE relaxation from extrapolated polarized IC
            if ~ok
                bvp_sol = pde_fallback_extrapolated(prev_sol, prev_prev_sol, ...
                                                    tg, tp, L_half, par, D_vals, p);
                ok = ~isempty(bvp_sol);
                tier_used = 2;
            end

            if ok
                bvp_sol = force_peak_at_zero(bvp_sol);   % canonical orientation (peak at x=0)
                if tier_used == 1, n_tier1 = n_tier1 + 1; else, n_tier2 = n_tier2 + 1; end
                SOL_grid{i, anchor_j}    = bvp_sol;
                failure_mask(i, anchor_j) = false;
                column_seeds{i}          = bvp_sol;
                column_succ(i)           = true;
                [U_trough(i,anchor_j), Q_trough(i,anchor_j), V_trough(i,anchor_j), ...
                 W_trough(i,anchor_j), P_trough(i,anchor_j)] = extract_trough(bvp_sol);
                prev_prev_sol = prev_sol;
                prev_sol      = bvp_sol;
            else
                column_seeds{i} = prev_sol;
                column_succ(i)  = false;
                fprintf('  -> Column i=%d (tg=%.2f) FAILED (predictor + PDE fallback).\n', i, tg);
            end
        end
    end
    p2_wall = toc(t_phase2);
    fprintf('  -> Column complete. Solved %d/%d (tier1=%d, tier2-PDE=%d), failed=%d.  Wall: %.2fs.\n', ...
            sum(column_succ)-1, res-1, n_tier1, n_tier2, sum(~column_succ), p2_wall);

    % ====================================================================
    % PHASE 3  -  PARALLEL ROW SWEEPS WITH PREDICTOR-CORRECTOR
    % ====================================================================
    fprintf('PHASE 3: Parallel row sweeps (predictor-corrector)...\n');

    parfor i = 1:res
        row_SOL  = cell(1, res);
        row_fail = true(1, res);
        row_U = nan(1, res); row_Q = nan(1, res);
        row_V = nan(1, res); row_W = nan(1, res); row_P = nan(1, res);

        row_seed = column_seeds{i};
        seed_ok  = column_succ(i);

        if seed_ok
            row_SOL{anchor_j}  = row_seed;
            row_fail(anchor_j) = false;
            [row_U(anchor_j), row_Q(anchor_j), row_V(anchor_j), ...
             row_W(anchor_j), row_P(anchor_j)] = extract_trough(row_seed);

            for dir_step = [1, -1]
                if dir_step > 0, j_range = anchor_j+1 : 1  : res;
                else,            j_range = anchor_j-1 : -1 : 1;
                end
                prev_sol      = row_seed;
                prev_prev_sol = [];
                for j = j_range
                    tg = n_GEF_bar + GEF_DIFF(i, j)/2;
                    tp = n_PAK_bar + PAK_DIFF(i, j)/2;

                    guess = predictor_step(prev_sol, prev_prev_sol);
                    [bvp_sol, ok] = run_bvp_v7(guess, L_half, par, D_vals, tg, tp);
                    ok = ok && validate_halfwave(bvp_sol);   % do not accept or seed from a collapsed (flat/no-peak) BVP
                    if ~ok
                        bvp_sol = lcreep_fallback(L_half, par, D_vals, tg, tp, p);
                        ok = ~isempty(bvp_sol);
                    end
                    if ok
                        bvp_sol = force_peak_at_zero(bvp_sol);   % canonical orientation (peak at x=0)
                        row_SOL{j}  = bvp_sol;
                        row_fail(j) = false;
                        [row_U(j), row_Q(j), row_V(j), row_W(j), row_P(j)] = extract_trough(bvp_sol);
                        prev_prev_sol = prev_sol;
                        prev_sol      = bvp_sol;
                    end
                end
            end
        end

        % Sliced outputs
        SOL_grid(i, :)     = row_SOL;
        failure_mask(i, :) = row_fail;
        U_trough(i, :)     = row_U;
        Q_trough(i, :)     = row_Q;
        V_trough(i, :)     = row_V;
        W_trough(i, :)     = row_W;
        P_trough(i, :)     = row_P;
    end

    trough_grids = {U_trough, Q_trough, V_trough, W_trough, P_trough};
    [SOL_grid, failure_mask, trough_grids, morphology_class_grid, ...
        morphology_reason_grid, morphology_diagnostics_grid, ...
        quarantined_nonhpss_mask, n_quarantined] = ...
        quarantine_invalid_hpss(SOL_grid, failure_mask, trough_grids, ...
            morphology_class_grid, morphology_reason_grid, ...
            morphology_diagnostics_grid, quarantined_nonhpss_mask);
    U_trough = trough_grids{1}; Q_trough = trough_grids{2};
    V_trough = trough_grids{3}; W_trough = trough_grids{4};
    P_trough = trough_grids{5};

    n_failed_after_main = sum(failure_mask(:));
    fprintf(['  -> Parallel sweep done. Currently failed: %d / %d; ', ...
             'morphology-quarantined before recovery: %d.\n'], ...
            n_failed_after_main, res*res, n_quarantined);

    save_checkpoint_pdefb(out_dir, computation_id, 'post_phase3', PAK_DIFF, GEF_DIFF, ...
        U_trough, Q_trough, V_trough, W_trough, P_trough, SOL_grid, failure_mask, ...
        morphology_class_grid, morphology_reason_grid, morphology_diagnostics_grid, ...
        quarantined_nonhpss_mask, L_total, M_BG, M_P, par, p);

    % ====================================================================
    % PHASE 4  -  NEIGHBOUR-SEEDED BVP FALLBACK (parallel, multi-pass)
    % ====================================================================
    if n_failed_after_main > 0
        fprintf('PHASE 4 [pdefb]: neighbour-seeded BVP fallback on %d failed point(s)...\n', ...
                n_failed_after_main);

    % The cap is generous: at res = 50 a pass limit of 5 stopped recovery while
        % the rate was still declining (66, 50, 43, 36, 32 cells per pass, not yet zero). The loop
        % early-stops via `if n_rescued_pass == 0, break; end`, so a high cap costs
        % nothing once the recovery is exhausted. No shape or mass threshold depends on it.
        max_passes      = 15;
        total_attempts  = 0;
        total_rescued   = 0;
        passes_used     = 0;
        t_phase4 = tic;
        pool = gcp;
        for pass = 1:max_passes
            fail_idx_list = find(failure_mask);
            n_failed_now  = numel(fail_idx_list);
            if n_failed_now == 0, break; end

    % Pre-extract only the (up to 2) neighbour seeds each failed cell
            % actually needs. This is a cheap serial pass over cell-array
            % pointers (no BVP solves) and avoids broadcasting the entire
            % growing SOL_grid (every solved cell so far, up to res*res
            % bvp4c structs) to every worker on every pass, which would
            % cause memory pressure at large res.
            seed1_list = cell(n_failed_now, 1);
            seed2_list = cell(n_failed_now, 1);
            tg_list    = zeros(n_failed_now, 1);
            tp_list    = zeros(n_failed_now, 1);
            for k_idx = 1:n_failed_now
                idx_lin = fail_idx_list(k_idx);
                [i, j]  = ind2sub([res, res], idx_lin);
                tg_list(k_idx) = n_GEF_bar + GEF_DIFF(i, j)/2;
                tp_list(k_idx) = n_PAK_bar + PAK_DIFF(i, j)/2;
                [seed1_list{k_idx}, seed2_list{k_idx}] = ...
                    find_neighbour_seeds(SOL_grid, failure_mask, i, j);
            end

            new_sols   = cell(n_failed_now, 1);
            seed_avail = ~cellfun(@isempty, seed1_list);
            task_idx   = find(seed_avail);
            n_tasks    = numel(task_idx);

    % Submit each solvable cell as its own parfeval task so a single
            % pathological BVP call can be individually timed out and
            % cancelled instead of blocking the whole batch indefinitely
            % (parfor cannot be interrupted per iteration; one stuck worker holds up all of
            % them, invisibly).
            futures = parallel.FevalFuture.empty(0, n_tasks);
            for t = 1:n_tasks
                k = task_idx(t);
                futures(t) = parfeval(pool, @bvp_seed_extrapolated, 1, ...
                    seed1_list{k}, seed2_list{k}, tg_list(k), tp_list(k), ...
                    L_half, par, D_vals);
            end

            pending = true(1, n_tasks);
            n_timed_out = 0;
            while any(pending)
                live = find(pending);
                [rel_idx, value] = fetchNext(futures(live), 2);
                if isempty(rel_idx)
                    for li = live
                        % Snapshot StartDateTime ONCE and match datetime('now')
                        % to its TimeZone: a future can transition queued->running
                        % (unzoned NaT -> zoned timestamp) between reads, and
                        % mixing a zoned and an unzoned datetime in the subtraction
                        % errors. Skip futures that have not actually started.
                        sd = futures(li).StartDateTime;
                        if isnat(sd), continue; end
                        tz = sd.TimeZone;
                        if isempty(tz), now_dt = datetime('now');
                        else, now_dt = datetime('now','TimeZone',tz); end
                        if seconds(now_dt - sd) > p.phase4_cell_timeout
                            cancel(futures(li));
                            pending(li)  = false;
                            n_timed_out  = n_timed_out + 1;
                        end
                    end
                    continue;
                end
                li = live(rel_idx);
                pending(li) = false;
                new_sols{task_idx(li)} = value;
            end

            n_rescued_pass = 0;
            for k_idx = 1:n_failed_now
                if ~isempty(new_sols{k_idx})
                    idx_lin = fail_idx_list(k_idx);
                    [i, j]  = ind2sub([res, res], idx_lin);
                    rescued_sol = force_peak_at_zero(new_sols{k_idx});   % canonical orientation
                    [accepted, morphology] = validate_halfwave(rescued_sol);
                    morphology_class_grid(i,j) = morphology.class;
                    morphology_reason_grid(i,j) = morphology.reason;
                    morphology_diagnostics_grid{i,j} = morphology;
                    if accepted
                        SOL_grid{i, j}     = rescued_sol;
                        failure_mask(i, j) = false;
                        [U_trough(i,j), Q_trough(i,j), V_trough(i,j), ...
                         W_trough(i,j), P_trough(i,j)] = extract_trough(rescued_sol);
                        n_rescued_pass = n_rescued_pass + 1;
                    else
                        quarantined_nonhpss_mask(i,j) = true;
                    end
                end
            end

            % Audit the complete seed pool before another
            % recovery pass can draw neighbours from it.
            trough_grids = {U_trough, Q_trough, V_trough, W_trough, P_trough};
            [SOL_grid, failure_mask, trough_grids, morphology_class_grid, ...
                morphology_reason_grid, morphology_diagnostics_grid, ...
                quarantined_nonhpss_mask, n_quarantined_pass] = ...
                quarantine_invalid_hpss(SOL_grid, failure_mask, trough_grids, ...
                    morphology_class_grid, morphology_reason_grid, ...
                    morphology_diagnostics_grid, quarantined_nonhpss_mask);
            U_trough = trough_grids{1}; Q_trough = trough_grids{2};
            V_trough = trough_grids{3}; W_trough = trough_grids{4};
            P_trough = trough_grids{5};
            n_rescued_pass = max(0, n_rescued_pass - n_quarantined_pass);
            total_attempts = total_attempts + n_failed_now;
            total_rescued  = total_rescued  + n_rescued_pass;
            passes_used    = pass;
            fprintf(['  Pass %d: %d cells, %d had polarized seed, rescued %d, ', ...
                     'timed out %d. Remaining: %d.\n'], ...
                    pass, n_failed_now, sum(seed_avail), n_rescued_pass, ...
                    n_timed_out, sum(failure_mask(:)));

            save_checkpoint_pdefb(out_dir, computation_id, sprintf('post_phase4_pass%d', pass), ...
                PAK_DIFF, GEF_DIFF, U_trough, Q_trough, V_trough, W_trough, P_trough, ...
                SOL_grid, failure_mask, morphology_class_grid, morphology_reason_grid, ...
                morphology_diagnostics_grid, quarantined_nonhpss_mask, ...
                L_total, M_BG, M_P, par, p);

            if n_rescued_pass == 0, break; end
        end
        p4_wall = toc(t_phase4);
        fprintf(['  -> Phase 4 done. Rescued %d / %d attempts across %d passes. ', ...
                 'Remaining: %d. Wall: %.2fs.\n'], ...
                total_rescued, total_attempts, passes_used, ...
                sum(failure_mask(:)), p4_wall);
    end

    % ====================================================================
    % PHASE 5  -  STABILITY CLASSIFICATION
    % ====================================================================
    if p.skip_stability
        fprintf('PHASE 5: Skipped (p.skip_stability = true). is_stable_grid left all-false.\n');
        stability_notes = sprintf('Failed=%d/%d, Stability=skipped', ...
                                  sum(failure_mask(:)), res*res);
    else
        fprintf('PHASE 5: Stability classification...\n');
        fail_flat   = failure_mask(:);
        SOL_flat    = SOL_grid(:);
        stable_flat = ~fail_flat;   % solved points default to stable; failed points stay false
        parfor idx = 1:res*res
            if ~fail_flat(idx)
                try
                    stable_flat(idx) = check_stability_v7(SOL_flat{idx}, par, D_vals);
                catch
                    % keep default: stable
                end
            end
        end
        is_stable_grid = reshape(stable_flat, res, res);
        fprintf('  -> Stable: %d / %d solved points.\n', ...
                sum(is_stable_grid(:)), sum(~failure_mask(:)));
        stability_notes = sprintf('Failed=%d/%d, Stable=%d', ...
                                  sum(failure_mask(:)), res*res, sum(is_stable_grid(:)));
    end

    % ====================================================================
    % SAVE
    % ====================================================================
    assert(all(cellfun(@isempty, SOL_grid(failure_mask))), ...
           'Failed or quarantined HPSS cells must not retain seedable solutions.');
    assert(all(morphology_class_grid(~failure_mask) == "HPSS"), ...
           'Every stored HPSS solution must pass the morphology gate.');
    save(out_path, 'PAK_DIFF', 'GEF_DIFF', ...
         'U_trough', 'Q_trough', 'V_trough', 'W_trough', 'P_trough', ...
         'SOL_grid', 'failure_mask', 'is_stable_grid', ...
         'morphology_class_grid', 'morphology_reason_grid', ...
         'morphology_diagnostics_grid', 'quarantined_nonhpss_mask', ...
         'L_total', 'M_BG', 'M_P', 'par', 'computation_id', 'p');
    copyfile(out_path, fullfile(out_dir, 'bvp_sweep_full_results.mat'));
    log_metadata(computation_id, filename, 'Hybrid sweep (PDE-fallback variant)', ...
                 par, L_total, M_BG, M_P, M_BG/L_total, M_P/L_total, ...
                 stability_notes);
    fprintf('Hybrid Sweep (pdefb) complete: %s\n', out_path);
end


function save_checkpoint_pdefb(out_dir, computation_id, tag, PAK_DIFF, GEF_DIFF, ...
        U_trough, Q_trough, V_trough, W_trough, P_trough, ...
        SOL_grid, failure_mask, morphology_class_grid, morphology_reason_grid, ...
        morphology_diagnostics_grid, quarantined_nonhpss_mask, ...
        L_total, M_BG, M_P, par, p)
    % Best-effort checkpoint so a stalled or killed run doesn't lose
    % completed work. Overwrites one file per run (keyed by computation_id),
    % not versioned per pass, so this stays cheap even with many passes.
    ckpt_path = fullfile(out_dir, sprintf('checkpoint_%s.mat', computation_id));
    try
        save(ckpt_path, 'PAK_DIFF', 'GEF_DIFF', ...
             'U_trough', 'Q_trough', 'V_trough', 'W_trough', 'P_trough', ...
             'SOL_grid', 'failure_mask', 'morphology_class_grid', ...
             'morphology_reason_grid', 'morphology_diagnostics_grid', ...
             'quarantined_nonhpss_mask', 'L_total', 'M_BG', 'M_P', 'par', ...
             'computation_id', 'p', 'tag');
        fprintf('  [checkpoint saved: %s, tag=%s]\n', ckpt_path, tag);
    catch ME
        warning('Checkpoint save failed (tag=%s): %s', tag, ME.message);
    end
end


function [SOL_grid, failure_mask, trough_grids, class_grid, reason_grid, ...
        diagnostics_grid, quarantine_mask, n_quarantined] = ...
        quarantine_invalid_hpss(SOL_grid, failure_mask, trough_grids, ...
            class_grid, reason_grid, diagnostics_grid, quarantine_mask)
    % Hard gate before propagation. A converged non-half-peak solution remains documented, but
    % it is cleared from the seed pool and treated as a failed half-peak solve.
    [class_grid, reason_grid, diagnostics_grid, invalid_mask] = ...
        audit_hpss_grid(SOL_grid, failure_mask, class_grid, reason_grid, ...
                        diagnostics_grid);
    n_quarantined = sum(invalid_mask(:));
    if n_quarantined == 0, return; end

    quarantine_mask = quarantine_mask | invalid_mask;
    failure_mask(invalid_mask) = true;
    SOL_grid(invalid_mask) = {[]};
    for k = 1:numel(trough_grids)
        grid = trough_grids{k};
        grid(invalid_mask) = NaN;
        trough_grids{k} = grid;
    end
end


%% ========================================================================
%  PHASE 0 / FALLBACK HELPERS
%  ========================================================================

function bins = candidate_anchor_bins(res, mid_idx, max_bins)
    % Ordered candidate anchor bins for the Phase-0/1 search. A centre-first
    % coarse 3x3 spread across the (Dn_PAK, Dn_GEF) grid (rows/cols at ~1/4,
    % 1/2, 3/4 of each axis), the rest ordered by distance from the centre so
    % nearer bins are tried before corners. Capped at max_bins (default 9).
    % Deliberately a SPREAD, not the immediate neighbours: if the whole central
    % region is off-branch (e.g. the centre mass sits past a fold, so no static
    % HPSS exists near it), a bin at a genuinely different mass partition is
    % still tried. For small res the unique quartiles collapse and fewer bins
    % are returned (still centre-first).
    q = unique(min(res, max(1, round(res * [0.25 0.50 0.75]))));
    [JJ, II] = meshgrid(q, q);                      % II = rows (GEF), JJ = cols (PAK)
    cand = unique([mid_idx mid_idx; II(:) JJ(:)], 'rows', 'stable');
    d2 = (cand(:,1) - mid_idx).^2 + (cand(:,2) - mid_idx).^2;
    [~, ord] = sort(d2);                            % centre (d2=0) first
    cand = cand(ord, :);
    bins = cand(1:min(max_bins, size(cand,1)), :);
end

function sol = force_peak_at_zero(sol)
    % Normalize a half-domain half-peak (HPSS) solution to the canonical orientation: peak at x=0,
    % trough at x=L_half. bvp4c returns the mirror-image solution
    % (peak at x=L_half) about half the time -- an EXACT solution too, since the
    % reaction is x-independent and the no-flux BCs + mass constraints are
    % symmetric under x -> L_half - x -- but build_branch_library then classifies
    % it as "unknown (far-end peak)", NOT HPSS, and un-normalized seeds let the
    % continuation drift onto that mirror branch. Flip the 5 species if u(0) <
    % u(end) (a no-op for a genuine boundary peak at x=0, and for a symmetric
    % interior peak). Derivative/integral rows are zeroed; downstream reads the
    % species, and any re-seed rebuilds them via construct_bvp_guess_v7. Same
    % guard as the quarter/MPSS sweeper's force_peak_at_zero.
    if sol.y(1,1) < sol.y(1,end)
        yr = zeros(size(sol.y));
        yr([1,3,5,7,9],:) = fliplr(sol.y([1,3,5,7,9],:));
        sol.y = yr;
    end
end

function [seed_sol, info] = phase0_lcreep_seed(par, D, L_target, rho_G, rho_P, p)
    % Build a polarized seed at L_target by:
    %   1) running PDE at one or more test L values to find a clean half-peak
    %      ("interactive" lets the user pick visually; "auto" picks the first
    %       test L that converges to a polarized state with amplitude >= 0.5),
    %   2) extracting the fundamental half-wave at length L_fund,
    %   3) creeping in L (UP or DOWN) from L_fund to L_target at fixed
    %      density (rho_G, rho_P).
    info = struct('n_steps', 0, 'L_start', NaN, 'mode', '');
    mode = lower(p.seeding_mode);
    info.mode = mode;

    % --- Constructed-guess cold start (p.phase0_construct) -------------------
    % Solve the HPSS BVP directly at L_target from analytic edge-peak guesses,
    % with NO PDE relaxation. This is the robust path for oscillatory regimes
    % where the static half-peak is dynamically UNSTABLE (the stable attractor
    % is a moving/shuttling peak), so PDE relaxation flees it and the
    % PDE-relaxation seed reports "amplitude ~0". bvp4c does not care about temporal stability
    % and converges to the (existing) steady half-peak (HPSS) state from a decent guess. If the
    % branch has folded away at these masses no guess validates -> fall through.
    if ~isfield(p, 'phase0_construct') || p.phase0_construct
        cs = phase0_construct_seed(par, D, L_target, rho_G, rho_P, p);
        if ~isempty(cs) && validate_halfwave(cs)
            fprintf('  Phase 0: constructed-guess seed at L=%.3f (amp=%.2f); no relaxation/L-creep needed.\n', ...
                    L_target, max(cs.y(1,:)) - min(cs.y(1,:)));
            seed_sol = cs;
            info.L_start = L_target; info.n_steps = 0; info.mode = [mode '+construct'];
            return;
        end
        fprintf('  Phase 0: constructed-guess cold start found no HPSS; falling back to PDE relaxation.\n');
    end

    switch mode
        case 'interactive'
            [bvp_at_Lfund, L_fund] = phase0_interactive_seed(par, D, rho_G, rho_P, p);
        case 'auto'
            [bvp_at_Lfund, L_fund] = phase0_auto_seed(par, D, L_target, rho_G, rho_P, p);
        otherwise
            error('Unknown p.seeding_mode = "%s" (use "interactive" or "auto")', mode);
    end
    if isempty(bvp_at_Lfund)
        seed_sol = []; return;
    end
    info.L_start = L_fund;

    % Validate the L_fund seed BEFORE creeping; if amplitude < 0.5, reject.
    if ~validate_halfwave(bvp_at_Lfund)
        warning(['Phase 0: BVP at L_fund=%.3f has amplitude < 0.5 (looks like HSS). ', ...
                 'Aborting; pick a different test L or rerun.'], L_fund);
        seed_sol = []; return;
    end

    % Creep L_fund -> L_target (could be UP or DOWN).
    if abs(L_target - L_fund) < 1e-9
        seed_sol = bvp_at_Lfund;
    else
        seed_sol = lcreep_to_target(bvp_at_Lfund, L_fund, L_target, rho_G, rho_P, par, D, p);
    end
    if isempty(seed_sol), return; end
    info.n_steps = max(1, ceil(abs(L_target - L_fund) / max(p.L_creep_step, 1e-9)));

    % Re-validate at L_target.  If amplitude collapsed during the creep,
    % reject so the caller knows we crossed the bifurcation.
    if ~validate_halfwave(seed_sol)
        warning(['Phase 0: amplitude collapsed during L-creep from %.3f to %.3f. ', ...
                 'Aborting (likely crossed the bifurcation).'], L_fund, L_target);
        seed_sol = [];
    end
end


function seed_sol = phase0_construct_seed(par, D, L_target, rho_G, rho_P, p)
    % Cold-start HPSS seed by CONSTRUCTING analytic edge-peak guesses (peak at
    % x=0, trough at x=L_target) and solving the 14-var BVP directly -- no PDE
    % relaxation. Sweeps a small grid of peak amplitude x width, keeps the
    % highest-amplitude solution that converges AND passes validate_halfwave.
    % Returns [] if none validate (e.g. the static branch has folded away at
    % these masses). Tested to recover the amplitude-10 half-peak state at the
    % M_P = 6 centre, where PDE relaxation (subcritical + unstable half-peak) fails.
    seed_sol = [];
    nG = rho_G * L_target; nP = rho_P * L_target;
    try
        [uh, qh, vh, wh, ph] = get_hss_local(nG, nP, L_target, par);
    catch
        return;
    end
    if any(~isfinite([uh qh vh wh ph])), return; end

    if isfield(p, 'phase0_construct_amps'),   amps = p.phase0_construct_amps;
    else, amps = [4, 8, 12, 16, 20]; end
    if isfield(p, 'phase0_construct_widths'), widths = p.phase0_construct_widths;
    else, widths = [0.1, 0.25, 0.5, 0.9]; end

    xh = linspace(0, L_target, 200);
    best = []; best_amp = -inf;
    for A = amps
        for wpk = widths
            g = exp(-(xh / wpk).^2);                 % edge peak at x=0
            ss.x = xh;
            yf = zeros(14, numel(xh));
            yf(1,:) = uh + A*g;                      % u  (membrane GEF)
            yf(3,:) = qh + 0.5*A*g;                  % q  (complex tracks u)
            yf(5,:) = max(1e-3, vh - 0.3*A*g);       % v  (cytosol depletes)
            yf(7,:) = max(1e-3, wh - 0.3*A*g);       % w
            yf(9,:) = max(1e-3, ph - 0.1*A*g);       % p
            ss.y = yf;
            [bvp, ok] = run_bvp_v7(construct_bvp_guess_v7(ss), L_target, par, D, nG, nP);
            if ok && validate_halfwave(bvp)
                a = max(bvp.y(1,:)) - min(bvp.y(1,:));
                if a > best_amp, best_amp = a; best = bvp; end
            end
        end
        % Early exit once a clearly-polarized peak is in hand (avoid scanning
        % the whole grid every run); the amp>=2 gate keeps only real half-peaks.
        if best_amp >= 2, break; end
    end
    seed_sol = best;
end


%% ========================================================================
%  PHASE 2 PDE-EXTRAPOLATION FALLBACK
%  ========================================================================

function sol_out = bvp_seed_extrapolated(prev_sol, prev_prev_sol, ...
                                          tg_target, tp_target, L_half, par, D)
    % Phase-4 fast fallback: build an extrapolated polarized guess from the
    % neighbour solution(s) and feed it directly to bvp4c with a tight NMax
    % so it can't run away on stiff inputs.
    sol_out = [];
    if isempty(prev_sol), return; end

    x_ic = prev_sol.x;
    y_ic = prev_sol.y;
    if ~isempty(prev_prev_sol)
        y_pp = zeros(size(y_ic));
        for k = 1:size(y_ic, 1)
            y_pp(k,:) = interp1(prev_prev_sol.x, prev_prev_sol.y(k,:), x_ic, ...
                                'pchip', 'extrap');
        end
        y_ic = y_ic + (y_ic - y_pp);
    end
    y_ic([1,3,5,7,9],:) = max(0, y_ic([1,3,5,7,9],:));

    guess_in.x = x_ic;
    guess_in.y = y_ic;
    guess = construct_bvp_guess_v7(guess_in);

    % Tighter NMax for Phase 4 to bound wall-time per cell.
    options = bvpset('RelTol', 1e-4, 'NMax', 1500);
    ws = warning('off', 'all');
    cleanup = onCleanup(@() warning(ws));
    try
        bvp_sol = bvp4c(@(x,y) ode_sys_v7(x,y,par,D), ...
                        @(ya,yb) bc_aug_v7(ya,yb,tg_target,tp_target), ...
                        guess, options);
    catch
        return;
    end
    if isempty(bvp_sol), return; end
    if length(bvp_sol.x) >= 1500, return; end
    if any(isnan(bvp_sol.y(:))) || any(isinf(bvp_sol.y(:))), return; end
    if min(bvp_sol.y([1,3,5,7,9],:), [], 'all') < -1e-3, return; end
    if ~validate_halfwave(bvp_sol), return; end
    sol_out = bvp_sol;
end


function sol_out = pde_fallback_extrapolated(prev_sol, prev_prev_sol, ...
                                             tg_target, tp_target, ...
                                             L_half, par, D, p)
    % Single-tier Phase-2 fallback.
    %
    %   Step 1: build an IC by extrapolating the polarized profile of the
    %           previous successful column point(s) one mass step further
    %           in the same direction (linear in y), with a positivity
    %           guard.  If only one previous point exists, just reuse it.
    %   Step 2: relax this IC under the reaction-diffusion PDE on [0, L_half]
    %           with reflecting BCs (no symmetry-breaking cos perturbation;
    %           the IC is already polarized).
    %   Step 3: feed the relaxed PDE state to bvp4c at the target mass.
    %
    % Returns [] on any failure or if the BVP collapses to flat (HSS).
    sol_out = [];
    if isempty(prev_sol), return; end

    % --- Step 1: extrapolated IC on prev_sol's mesh ------------------------
    x_ic = prev_sol.x;
    y_ic = prev_sol.y;
    if ~isempty(prev_prev_sol)
        y_pp = zeros(size(y_ic));
        for k = 1:size(y_ic, 1)
            y_pp(k,:) = interp1(prev_prev_sol.x, prev_prev_sol.y(k,:), x_ic, ...
                                'pchip', 'extrap');
        end
        y_ic = y_ic + (y_ic - y_pp);
    end
    y_ic([1,3,5,7,9],:) = max(0, y_ic([1,3,5,7,9],:));     % positivity guard

    % --- Step 2: PDE relaxation on a fresh uniform mesh --------------------
    N = 100;
    x_pde = linspace(0, L_half, N);
    y0 = zeros(5, N);
    rows = [1, 3, 5, 7, 9];
    for k = 1:5
        y0(k,:) = interp1(x_ic, y_ic(rows(k),:), x_pde, 'pchip');
    end
    y0 = max(0, y0);

    max_ch = max(2, ceil(p.pde_max_chunks/2));

    if p.osc_fallback
        % --- Steps 2-3, oscillation-aware: bounded relaxation returning seed
        %     candidates, each tried through the SAME BVP core + validation ---
        cands = pde_relax_candidates(y0, x_pde, L_half, par, D, ...
                                     p.pde_fallback_T, p.pde_T_chunk, max_ch, p);
        sol_out = seed_bvp_from_candidates(cands, L_half, par, D, ...
                                           tg_target, tp_target, @(b) true);
        return;
    end

    % --- Step 2: PDE relaxation to a single steady state -------------------
    pde_sol = run_pde_polarised(y0, x_pde, L_half, par, D, ...
                                p.pde_fallback_T, p.pde_T_chunk, max_ch);
    if isempty(pde_sol), return; end

    % --- Step 3: BVP at the target mass, seeded by the relaxed PDE state --
    [bvp_sol, ok] = run_bvp_v7(construct_bvp_guess_v7(pde_sol), ...
                                L_half, par, D, tg_target, tp_target);
    if ~ok, return; end
    if ~validate_halfwave(bvp_sol), return; end
    sol_out = bvp_sol;
end


function cands = pde_relax_candidates(y0, x, ~, par, D, T_end, T_chunk, max_chunks, p)
    % Bounded, oscillation-aware relaxation. Returns an ORDERED cell array of
    % candidate seed structs (.x, 14-row .y in the construct_bvp_guess_v7
    % layout) to try as bvp4c seeds -- NEVER solutions. If the relaxation
    % converges to a steady state it returns just that state (identical seed to
    % the old single-shot path). If it does NOT relax (sustained oscillation),
    % it stops early and returns the time-mean of the last chunk's snapshots
    % first (best guess for the unstable fixed point near a Hopf), then each
    % phase snapshot. Nothing here decides acceptance -- see
    % seed_bvp_from_candidates.
    cands = {};
    if T_end <= 0, return; end
    if T_chunk > T_end, T_chunk = T_end; end
    if T_chunk <= 0, T_chunk = T_end/2; end
    if max_chunks < 1, max_chunks = 1; end
    if size(y0, 2) == 1, y0 = repmat(y0, 1, length(x)); end
    nsnap = max(2, p.osc_n_snap);

    curr_y = y0;
    deltas = [];
    last_frames = {};
    converged = false;
    ws = warning('off', 'all'); cleanup = onCleanup(@() warning(ws)); %#ok<NASGU>
    try
        for k = 1:max_chunks
            tspan = linspace(0, T_chunk, nsnap + 1);   % t=0 plus nsnap sub-frames
            sol = pdepe(0, ...
                @(xq,t,uu,dudx) pde_eq_v7(xq,t,uu,dudx,par,D), ...
                @(xq) pde_ic_local(xq, curr_y, x), ...
                @pde_bc_v7, x, tspan);
            last_frames = cell(1, nsnap);
            for s = 1:nsnap
                last_frames{s} = squeeze(sol(s+1,:,:))';   % 5 x N
            end
            frame_end = last_frames{end};
            d = max(abs(frame_end(1,:) - curr_y(1,:)));
            deltas(end+1) = d; %#ok<AGROW>
            curr_y = frame_end;
            if d < p.osc_conv_tol
                converged = true; break;                    % genuine steady state
            end
            if k >= p.osc_min_chunks && ...
               min(deltas(max(1,end-p.osc_min_chunks+1):end)) > p.osc_rel * deltas(1)
                break;                                       % not relaxing -> stop early (bounded)
            end
        end
    catch
        cands = {}; return;
    end

    if converged
        cands = { state5_to_seed(curr_y, x) };
        return;
    end
    if isempty(last_frames), return; end
    ymean = zeros(size(last_frames{1}));
    for s = 1:numel(last_frames), ymean = ymean + last_frames{s}; end
    ymean = max(0, ymean / numel(last_frames));
    cands{end+1} = state5_to_seed(ymean, x);                 % time-mean first
    for s = 1:numel(last_frames)
        cands{end+1} = state5_to_seed(max(0, last_frames{s}), x); %#ok<AGROW>
    end
end


function ss = state5_to_seed(y5, x)
    % 5xN species -> seed struct with 14-row .y in the construct_bvp_guess_v7
    % layout (species in rows 1,3,5,7,9; derivative/integral rows filled there).
    ss.x = x;
    yf = zeros(14, size(y5, 2));
    yf(1,:) = y5(1,:); yf(3,:) = y5(2,:); yf(5,:) = y5(3,:);
    yf(7,:) = y5(4,:); yf(9,:) = y5(5,:);
    ss.y = yf;
end


function bvp_out = seed_bvp_from_candidates(cands, L, par, D, tg, tp, mass_gate)
    % Try each candidate seed through the shared BVP core; accept the FIRST that
    % converges AND passes validate_halfwave AND mass_gate (@(b) true to skip).
    % Returns [] if none pass -- an oscillatory snapshot only becomes a solution
    % when bvp4c converges from it to a validated steady half-peak.
    bvp_out = [];
    for c = 1:numel(cands)
        [bvp, ok] = run_bvp_v7(construct_bvp_guess_v7(cands{c}), L, par, D, tg, tp);
        if ok && validate_halfwave(bvp) && mass_gate(bvp)
            bvp_out = bvp; return;
        end
    end
end


function sol_out = run_pde_polarised(y0, x, L, par, D, T_end, T_chunk, max_chunks)
    % Like run_pde_chunked but does NOT add the symmetry-breaking cos
    % perturbation: the IC is already polarized (it's a half-peak BVP
    % profile) and we want the PDE to relax that exact shape, not flip it.
    sol_out = [];
    if T_end <= 0, return; end
    if T_chunk > T_end, T_chunk = T_end; end
    if T_chunk <= 0, T_chunk = T_end / 2; end
    if max_chunks < 1, max_chunks = 1; end

    if size(y0, 2) == 1
        y0 = repmat(y0, 1, length(x));
    end

    curr_y = y0;
    ws = warning('off', 'all');
    cleanup = onCleanup(@() warning(ws));
    try
        for k = 1:max_chunks
            tspan = [0, T_chunk/2, T_chunk];
            sol = pdepe(0, ...
                @(xq,t,uu,dudx) pde_eq_v7(xq,t,uu,dudx,par,D), ...
                @(xq) pde_ic_local(xq, curr_y, x), ...
                @pde_bc_v7, x, tspan);
            frame_end = squeeze(sol(end,:,:))';
            if max(abs(frame_end(1,:) - curr_y(1,:))) < 1e-6
                curr_y = frame_end;
                break;
            end
            curr_y = frame_end;
        end
    catch
        return;
    end

    sol_out.x = x;
    yf = zeros(14, length(x));
    yf(1,:) = curr_y(1,:); yf(3,:) = curr_y(2,:); yf(5,:) = curr_y(3,:);
    yf(7,:) = curr_y(4,:); yf(9,:) = curr_y(5,:);
    sol_out.y = yf;
end


function [bvp_sol, L_fund] = phase0_interactive_seed(par, D, rho_G, rho_P, p)
    % Interactive seed picker.  Runs PDE at each L in
    % p.L_seed_test_values, plots them in a panel, asks the user to pick the
    % index of the cleanest pattern and how many half-waves are visible.
    L_test = p.L_seed_test_values(:).';   % row vector
    n_test = numel(L_test);
    pde_results = cell(n_test, 1);

    fig = figure('Name', 'Phase 0: Interactive Seed Selection', 'Color', 'w', ...
                 'Position', [50 100 1600 600]);

    fprintf('  Running PDE at %d test L values (this may take a few minutes)...\n', n_test);
    for k = 1:n_test
        L_k = L_test(k);
        nG_k = rho_G * L_k; nP_k = rho_P * L_k;
        try
            [u, q, v, w, p_ss] = get_hss_local(nG_k, nP_k, L_k, par);
        catch
            pde_results{k} = []; continue;
        end
        if any(~isfinite([u q v w p_ss])), pde_results{k} = []; continue; end

        x_k = linspace(0, L_k, 100);
        pde_sol = run_pde_chunked([u; q; v; w; p_ss], x_k, L_k, par, D, ...
                                   p.pde_T_end, p.pde_T_chunk, p.pde_max_chunks);
        pde_results{k} = pde_sol;

        % Live plot
        figure(fig);
        n_cols = ceil(n_test/2);
        subplot(2, n_cols, k);
        if ~isempty(pde_sol)
            u_prof = pde_sol.y(1,:);
            p_prof = pde_sol.y(9,:);
            yyaxis left
            plot(pde_sol.x, u_prof, 'b-', 'LineWidth', 2);
            ylabel('u');
            yyaxis right
            plot(pde_sol.x, p_prof, 'r--', 'LineWidth', 1);
            ylabel('p');
            amp_u = max(u_prof) - min(u_prof);
            title(sprintf('[%d] L=%.2f  A_u=%.3f', k, L_k, amp_u));
            xlabel('x'); grid on;
        else
            text(0.5, 0.5, 'PDE failed', 'HorizontalAlignment','center');
            title(sprintf('[%d] L=%.2f', k, L_k));
            axis off;
        end
        drawnow;
    end

    % Prompt the user
    fprintf('\n  *** Inspect the figure window. ***\n');
    fprintf('  Each subplot shows u(x) (blue, left axis) and p(x) (red dashed, right axis).\n');
    fprintf('  Pick the subplot index whose profile shows a clean polarized pattern.\n');
    fprintf('  (Patterns with amplitude near 0 are HSS; ignore them.)\n');

    L_fund = NaN; bvp_sol = [];
    while true
        idx_str = input(sprintf('  Pick subplot index 1..%d (0 = abort): ', n_test), 's');
        idx_choice = str2double(idx_str);
        if ~isnan(idx_choice) && idx_choice == 0
            close(fig); return;
        end
        if ~isnan(idx_choice) && idx_choice >= 1 && idx_choice <= n_test ...
                && ~isempty(pde_results{idx_choice})
            chosen_idx = round(idx_choice);
            break;
        end
        fprintf('  Invalid selection.\n');
    end

    while true
        nh_str = input('  How many half-waves are in this profile? [1,2,3,...] (Enter = 1): ', 's');
        if isempty(nh_str)
            n_half = 1; break;
        end
        n_half = str2double(nh_str);
        if ~isnan(n_half) && n_half >= 1 && abs(n_half - round(n_half)) < 1e-9
            n_half = round(n_half); break;
        end
        fprintf('  Enter a positive integer.\n');
    end

    L_chosen   = L_test(chosen_idx);
    L_fund     = L_chosen / n_half;
    chosen_sol = pde_results{chosen_idx};
    fprintf('  -> L_chosen = %.3f, %d half-wave(s), L_fund = %.3f.\n', ...
            L_chosen, n_half, L_fund);

    % Extract the fundamental half-wave (slice [0, L_fund] of the chosen profile)
    fund_sol = extract_fundamental_halfwave(chosen_sol, L_fund);

    nG_f = rho_G * L_fund; nP_f = rho_P * L_fund;
    [bvp_sol, ok] = run_bvp_v7(construct_bvp_guess_v7(fund_sol), L_fund, par, D, nG_f, nP_f);
    if ~ok
        fprintf('  -> BVP at L_fund=%.3f failed; cannot use this seed.\n', L_fund);
        bvp_sol = [];
    end
    close(fig);
end


function [bvp_sol, L_fund] = phase0_auto_seed(par, D, L_target, rho_G, rho_P, p)
    % Auto seed: try a small list of L values, accept the first one that gives
    % a polarized profile (amplitude >= 0.5 in u) AND a converged BVP.
    L_candidates = unique([p.L_creep_start, L_target*0.5, L_target*0.7, ...
                           L_target, L_target*1.2, L_target*1.5], 'stable');
    L_candidates = L_candidates(L_candidates > 0);

    bvp_sol = []; L_fund = NaN;
    for L_try = L_candidates
        nG_t = rho_G * L_try; nP_t = rho_P * L_try;
        try
            [u, q, v, w, p_ss] = get_hss_local(nG_t, nP_t, L_try, par);
        catch
            continue;
        end
        if any(~isfinite([u q v w p_ss])), continue; end

        x_t = linspace(0, L_try, 100);
        pde_sol = run_pde_chunked([u; q; v; w; p_ss], x_t, L_try, par, D, ...
                                   p.pde_T_end, p.pde_T_chunk, p.pde_max_chunks);
        if isempty(pde_sol), continue; end

        amp = max(pde_sol.y(1,:)) - min(pde_sol.y(1,:));
        fprintf('  Auto seed: tried L=%.2f, amplitude=%.3f.\n', L_try, amp);
        if amp < 0.5, continue; end

        [bvp_try, ok] = run_bvp_v7(construct_bvp_guess_v7(pde_sol), L_try, par, D, nG_t, nP_t);
        if ok && validate_halfwave(bvp_try)
            bvp_sol = bvp_try; L_fund = L_try;
            return;
        end
    end
end


function fund_sol = extract_fundamental_halfwave(chosen_sol, L_fund)
    % Slice [0, L_fund] of the PDE profile, force edge-peaked orientation.
    x_old = chosen_sol.x;
    if L_fund > x_old(end) + 1e-9
        % Sanity: shouldn't happen if user supplied a sane n_half.
        L_fund = x_old(end);
    end
    x_new = linspace(0, L_fund, 100);

    % Pull the five species; the rest of the 14-row matrix can be left zero —
    % construct_bvp_guess_v7 rebuilds derivatives + integrals from scratch.
    rows = [1, 3, 5, 7, 9];
    species = zeros(numel(rows), numel(x_new));
    for r = 1:numel(rows)
        species(r,:) = interp1(x_old, chosen_sol.y(rows(r),:), x_new, 'pchip');
    end

    % Force edge-peaked: u(0) > u(end). Otherwise mirror.
    if species(1, 1) < species(1, end)
        species = fliplr(species);
    end
    species = max(0, species);   % positivity guard

    y_new = zeros(14, numel(x_new));
    y_new(rows, :) = species;
    fund_sol.x = x_new;
    fund_sol.y = y_new;
end

function sol_out = lcreep_fallback(L_target, par, D, nG, nP, p)
    % Per-point L-creep at the (i,j)-implied density; replaces the brittle
    % fresh-PDE fallback.  Density: rho_eff = nG/L_target, rho_P_eff = nP/L_target
    rho_G_eff = nG / L_target;
    rho_P_eff = nP / L_target;

    L_start = min(p.L_creep_start, L_target * 0.5);
    if L_start <= 0, L_start = L_target * 0.3; end
    nG_s = rho_G_eff * L_start; nP_s = rho_P_eff * L_start;

    try
        [u, q, v, w, p_ss] = get_hss_local(nG_s, nP_s, L_start, par);
    catch
        sol_out = []; return;
    end
    if any(~isfinite([u q v w p_ss])), sol_out = []; return; end

    x = linspace(0, L_start, 80);
    pde_sol = run_pde_chunked([u; q; v; w; p_ss], x, L_start, par, D, ...
                               p.pde_fallback_T, p.pde_T_chunk, ...
                               max(2, ceil(p.pde_max_chunks/2)));
    if isempty(pde_sol), sol_out = []; return; end

    [bvp_sol, ok] = run_bvp_v7(construct_bvp_guess_v7(pde_sol), ...
                                L_start, par, D, nG_s, nP_s);
    if ~ok, sol_out = []; return; end

    % Reject if the small-L BVP is itself near-flat (no PSS branch at this density).
    if ~validate_halfwave(bvp_sol), sol_out = []; return; end

    sol_out = lcreep_to_target(bvp_sol, L_start, L_target, rho_G_eff, rho_P_eff, par, D, p);
    if isempty(sol_out), return; end
    % Final amplitude check: if creep crossed back into HSS, mark as failed.
    if ~validate_halfwave(sol_out), sol_out = []; end
end

function sol_out = lcreep_to_target(sol_in, L_start, L_target, rho_G, rho_P, par, D, p)
    % Drive sol from L_start up to L_target in steps of p.L_creep_step.
    L_step = max(p.L_creep_step, 1e-3);
    if abs(L_target - L_start) < 1e-9
        sol_out = sol_in; return;
    end
    if L_target > L_start
        L_vals = (L_start + L_step) : L_step : L_target;
        if isempty(L_vals) || L_vals(end) < L_target - 1e-9
            L_vals = [L_vals, L_target];
        end
    else
        L_vals = (L_start - L_step) : -L_step : L_target;
        if isempty(L_vals) || L_vals(end) > L_target + 1e-9
            L_vals = [L_vals, L_target];
        end
    end
    curr = sol_in;
    for L_now = L_vals
        nG_now = rho_G * L_now;
        nP_now = rho_P * L_now;
        guess  = rescale_solution_v7(curr, L_now);
        guess  = construct_bvp_guess_v7(guess);
        [bvp_sol, ok] = run_bvp_v7(guess, L_now, par, D, nG_now, nP_now);
        if ~ok
            sol_out = []; return;
        end
        % Crossing-the-bifurcation guard: if the amplitude collapses to HSS
        % during the creep, stop immediately rather than silently propagating
        % a flat solution to the target.
        if ~validate_halfwave(bvp_sol)
            sol_out = []; return;
        end
        curr = bvp_sol;
    end
    sol_out = curr;
end


%% ========================================================================
%  BVP CORE  (per-point solver; mass integrals are not rescaled)
%  ========================================================================





function sol_scaled = rescale_solution_v7(sol_in, L_new)
    % Stretch the mesh to L_new, scale derivatives accordingly. Integrals are
    % rebuilt by construct_bvp_guess_v7 and the s1,s2 channel handles any drift.
    sol_scaled = sol_in;
    L_old = sol_in.x(end);
    if L_old <= 0
        sol_scaled.x = linspace(0, L_new, length(sol_in.x));
        return;
    end
    ratio = L_new / L_old;
    sol_scaled.x = sol_in.x * ratio;
    if size(sol_in.y, 1) >= 10
        sol_scaled.y([2,4,6,8,10], :) = sol_in.y([2,4,6,8,10], :) / ratio;
    end
end


%% ========================================================================
%  PREDICTOR-CORRECTOR
%  ========================================================================

function guess = predictor_step(prev_sol, prev_prev_sol)
    % 2-point linear extrapolation when both previous solutions are available;
    % otherwise use the most recent solution directly. All BVPs in this sweep
    % share L = L_half, so the meshes can be aligned by linear resampling.
    if isempty(prev_prev_sol)
        guess = construct_bvp_guess_v7(prev_sol);
        return;
    end

    x_ref  = prev_sol.x;
    y_prev = prev_sol.y;
    y_pp   = zeros(size(y_prev));
    for k = 1:size(y_prev, 1)
        y_pp(k,:) = interp1(prev_prev_sol.x, prev_prev_sol.y(k,:), x_ref, 'pchip', 'extrap');
    end
    y_pred = y_prev + (y_prev - y_pp);    % one extra step in the same direction
    y_pred([1,3,5,7,9],:) = max(0, y_pred([1,3,5,7,9],:));   % positivity guard
    pred.x = x_ref;
    pred.y = y_pred;
    guess  = construct_bvp_guess_v7(pred);
end


%% ========================================================================
%  CHUNKED PDE
%  ========================================================================






%% ========================================================================
%  HSS ALGEBRAIC SOLVE (used to seed the chunked PDE only)
%  ========================================================================




%% ========================================================================
%  SEED VALIDATION
%  ========================================================================



%% ========================================================================
%  STABILITY (eigenvalues of the BVP Jacobian)
%  ========================================================================

function is_stable = check_stability_v7(sol, par, D)
    if isempty(sol) || ~isfield(sol, 'y'), is_stable = false; return; end
    N = 500;
    x_uni = linspace(min(sol.x), max(sol.x), N);
    dx = x_uni(2) - x_uni(1);

    U = interp1(sol.x, sol.y(1,:), x_uni, 'pchip')';
    Q = interp1(sol.x, sol.y(3,:), x_uni, 'pchip')';
    V = interp1(sol.x, sol.y(5,:), x_uni, 'pchip')';
    W = interp1(sol.x, sol.y(7,:), x_uni, 'pchip')';
    P = interp1(sol.x, sol.y(9,:), x_uni, 'pchip')';

    dRu_du = 2*par.a.*U.*V - par.b - par.c.*P;
    dRu_dv = par.a.*U.^2 + par.f;
    dRu_dp = -par.c.*U;
    dRq_du = par.c.*P;
    dRq_dq = -par.d * ones(N,1);
    dRq_dp = par.c.*U;
    dRv_du = -(2*par.a.*U.*V - par.b);
    dRv_dv = -(par.a.*U.^2 + par.f);
    dRv_dw = par.e * ones(N,1);
    dRw_dq = par.d * ones(N,1);
    dRw_dw = -par.e * ones(N,1);
    dRp_du = -par.c.*P;
    dRp_dq = par.d * ones(N,1);
    dRp_dp = -par.c.*U;

    e  = ones(N,1);
    D2 = spdiags([e -2*e e], -1:1, N, N) / dx^2;
    D2(1,2)     = 2/dx^2;
    D2(end,end-1) = 2/dx^2;
    Z  = sparse(N, N);

    J = [D(1)*D2 + spdiags(dRu_du,0,N,N), Z, spdiags(dRu_dv,0,N,N), Z, spdiags(dRu_dp,0,N,N);
         spdiags(dRq_du,0,N,N), D(2)*D2 + spdiags(dRq_dq,0,N,N), Z, Z, spdiags(dRq_dp,0,N,N);
         spdiags(dRv_du,0,N,N), Z, D(3)*D2 + spdiags(dRv_dv,0,N,N), spdiags(dRv_dw,0,N,N), Z;
         Z, spdiags(dRw_dq,0,N,N), Z, D(4)*D2 + spdiags(dRw_dw,0,N,N), Z;
         spdiags(dRp_du,0,N,N), spdiags(dRp_dq,0,N,N), Z, Z, D(5)*D2 + spdiags(dRp_dp,0,N,N)];

    evs = eig(full(J));
    is_stable = (max(real(evs)) <= 1e-4);
end


%% ========================================================================
%  TROUGH EXTRACTION
%  ========================================================================

function [u_t, q_t, v_t, w_t, p_t] = extract_trough(sol)
    u_t = sol.y(1, end);
    q_t = sol.y(3, end);
    v_t = sol.y(5, end);
    w_t = sol.y(7, end);
    p_t = sol.y(9, end);
end


function [seed1, seed2] = find_neighbour_seeds(sol_grid, fail_mask, i, j)
    % Pick a POLARIZED neighbour seed from (i,j)'s cardinal grid neighbours.
    % Flat (HSS-collapsed) neighbours are skipped: a PDE relax of a flat IC
    % stays flat, so seeding from one is wasted work that always fails the
    % validate_halfwave check at the end.
    %
    % Preference: a collinear pair (closer + further in the same axis) for
    % linear extrapolation toward (i,j). Fall back to the closest single
    % polarized cardinal neighbour if no pair is available.
    [ni, nj] = size(fail_mask);
    seed1 = []; seed2 = [];

    % Collinear pairs first; both must be polarized.
    if i >= 3 && is_polarized(sol_grid, fail_mask, i-1, j) && ...
                 is_polarized(sol_grid, fail_mask, i-2, j)
        seed1 = sol_grid{i-1, j}; seed2 = sol_grid{i-2, j}; return;
    end
    if i <= ni-2 && is_polarized(sol_grid, fail_mask, i+1, j) && ...
                    is_polarized(sol_grid, fail_mask, i+2, j)
        seed1 = sol_grid{i+1, j}; seed2 = sol_grid{i+2, j}; return;
    end
    if j >= 3 && is_polarized(sol_grid, fail_mask, i, j-1) && ...
                 is_polarized(sol_grid, fail_mask, i, j-2)
        seed1 = sol_grid{i, j-1}; seed2 = sol_grid{i, j-2}; return;
    end
    if j <= nj-2 && is_polarized(sol_grid, fail_mask, i, j+1) && ...
                    is_polarized(sol_grid, fail_mask, i, j+2)
        seed1 = sol_grid{i, j+1}; seed2 = sol_grid{i, j+2}; return;
    end

    % Single-neighbour fallback: closest polarized cardinal cell.
    cand = {[i-1,j], [i+1,j], [i,j-1], [i,j+1]};
    for k = 1:numel(cand)
        ii = cand{k}(1); jj = cand{k}(2);
        if ii >= 1 && ii <= ni && jj >= 1 && jj <= nj && ...
                is_polarized(sol_grid, fail_mask, ii, jj)
            seed1 = sol_grid{ii, jj}; return;
        end
    end
end


function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end


function tf = is_polarized(sol_grid, fail_mask, i, j)
    tf = false;
    if i < 1 || j < 1 || i > size(fail_mask,1) || j > size(fail_mask,2), return; end
    if fail_mask(i, j), return; end
    s = sol_grid{i, j};
    if isempty(s) || ~isfield(s, 'y'), return; end
    u = s.y(1, :);
    tf = (max(u) - min(u)) > 1e-3;
end
