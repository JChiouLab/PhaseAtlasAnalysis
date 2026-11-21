function result = run_quarter_hpss_sweeper_for_mpss(p, out_dir)
    % RUN_QUARTER_HPSS_SWEEPER_FOR_MPSS  Continuation-based quarter-domain half-peak (HPSS)
    % sweep, the search engine for the full-peak (MPSS) constructor build_mpss_branch.m.
    %
    % Why this exists: solving every quarter-domain half-peak cell as an independent cold start
    % (no-peak state -> cos-perturbed PDE relaxation -> BVP) fails wherever the relaxation does
    % not land on the polarized branch at that density, so the accepted full-peak region came out
    % fragmented. This function instead uses the continuation strategy of the half-peak sweeper:
    % seed one centre cell, then predictor-corrector continue across the grid so each cell is
    % started from a neighbour already on the branch.
    %
    % What it produces: a grid of quarter-domain half-peak profiles (14-variable bvp4c structs,
    % peak forced to x = 0) on [0, L_quarter], indexed like the (Delta n_PAK, Delta n_GEF) grid,
    % at half-scaled masses. build_mpss_branch.m reflects and stitches each solved quarter profile
    % peak-to-peak into a full-peak profile; the full-peak state is constructed, never solved as a
    % direct interior-peak BVP.
    %
    % Mass-grid convention:
    %   grid cell (i,j):   nG_half = M_BG/2 + GEF_DIFF(i,j)/2
    %                      nP_half = M_P /2 + PAK_DIFF(i,j)/2
    %   quarter target:    nG_q    = nG_half/2 ,  nP_q = nP_half/2
    %   solve domain:      L_quarter = L_total/4
    % Same grid indices as the half-domain grid; half the absolute mass range; twice the
    % mass density. (Equivalent to running the half-peak sweeper with M_BG/2, M_P/2, L/2.)
    %
    % Reuse: calls only the shared BVP core (run_bvp_v7, construct_bvp_guess_v7,
    % run_pde_chunked, get_hss_local, validate_halfwave). The continuation orchestration
    % (predictor-corrector, neighbour retry) is a self-contained copy of the pattern in
    % run_hybrid_sweeper_pdefb.m. Writes only its own files; never touches the half-peak or
    % no-peak library files.
    %
    % Controls:
    %  * checkpoint after Phase 3 and every Phase 4 pass -> checkpoint_quarter_hpss_<id>.mat
    %  * resume via p.quarter_resume_from (a quarter checkpoint); solved cells are skipped
    %  * Phase 4 retry: parfeval + per-cell timeout + cancel (no growing broadcast;
    %    only the 1-2 neighbour seeds each cell needs are pre-sliced)
    %  * failed cells are masked, never filled in
    %
    % Output: quarter_hpss_for_mpss_<id>.mat in out_dir (+ checkpoint). Returns a
    % struct with SOL_grid (quarter half-peak structs
    % atlas grid, and provenance.

    % ---- parameters ----------------------------------------------------------
    if nargin < 1 || isempty(p)
        p.a = 0.59; p.b = 0.13; p.c = 21.0;
        p.d = 1.3;  p.e = 33.0; p.f = 1.9;
        p.Dm = 0.001; p.Dc = 1.0;
        p.L = 4.0; p.M_BG = 16.0; p.M_P = 8.0;
        p.res = 3;
    end
    if nargin < 2 || isempty(out_dir), out_dir = 'data/outputs'; end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    if ~isfield(p, 'res'),            p.res = 3;       end
    if ~isfield(p, 'pde_T_end'),      p.pde_T_end   = 2e6;  end
    if ~isfield(p, 'pde_T_chunk'),    p.pde_T_chunk = 1e6;  end
    if ~isfield(p, 'pde_max_chunks'), p.pde_max_chunks = 12; end
    if ~isfield(p, 'pde_fallback_T'), p.pde_fallback_T = max(1, p.pde_T_end/4); end
    if ~isfield(p, 'mass_tol'),       p.mass_tol = 1e-4; end
    % Oscillatory-regime PDE fallback (default OFF).
    % Same design as run_hybrid_sweeper_pdefb.m: bounded relaxation with early
    % non-relaxation detection + oscillation-phase snapshots as bvp4c SEED
    % candidates. Snapshots are seeds only; validate_halfwave and the quarter
    % mass gate still apply to every accepted solve.
    if ~isfield(p, 'osc_fallback'),   p.osc_fallback   = false; end
    if ~isfield(p, 'osc_conv_tol'),   p.osc_conv_tol   = 1e-6;  end
    if ~isfield(p, 'osc_min_chunks'), p.osc_min_chunks = 3;     end
    if ~isfield(p, 'osc_rel'),        p.osc_rel        = 0.5;   end
    if ~isfield(p, 'osc_n_snap'),     p.osc_n_snap     = 6;     end
    % Accept a quarter solve only if its RECOMPUTED integral mass is on target
    % (the 14-variable BVP's s1,s2 sources can absorb a mass
    % discrepancy, so a shape-valid half-peak may still sit at the wrong mass;
    % such a solve must not be stored OR propagated as a continuation seed).
    % Gate at mass_tol/2 so the reflected (doubled) MPSS mass stays within
    % mass_tol in classify_mpss_profile. Off-mass cells fall through to the
    % cold-start rescue floor instead of polluting the field.
    if ~isfield(p, 'quarter_mass_tol'), p.quarter_mass_tol = p.mass_tol/2; end
    % Continuation controls
    if ~isfield(p, 'use_parallel'),        p.use_parallel = true;  end
    if ~isfield(p, 'cell_timeout'),        p.cell_timeout = 60;    end   % s, rescue only
    if ~isfield(p, 'checkpoint_every'),    p.checkpoint_every = 25; end
    if ~isfield(p, 'quarter_resume_from'), p.quarter_resume_from = ''; end
    if ~isfield(p, 'rescue_max_passes'),   p.rescue_max_passes = 5; end
    if ~isfield(p, 'max_reroot_tries'),    p.max_reroot_tries = p.res*p.res; end
    % Phase 3 (row continuation) bounds each row by a time limit: a plain
    % parfor blocks until every row returns, so one row stuck near a fold
    % (repeated failed predictor-corrector -> PDE fallback across many cells)
    % could stall the whole phase for hours with no progress
    % visible. A timed-out row's unset cells simply stay in failure_mask; Phase 4's
    % retry (also parfeval+timeout-guarded) picks them up via neighbour/cold
    % start, so nothing is lost, only deferred.
    if ~isfield(p, 'phase3_row_timeout'),  p.phase3_row_timeout = 120; end  % s, use_parallel only

    par = p;
    D   = [p.Dm; p.Dm; p.Dc; p.Dc; p.Dc];
    res = p.res;
    L_total   = p.L;
    L_half    = L_total / 2;
    L_quarter = L_total / 4;
    n_GEF_bar = p.M_BG / 2;
    n_PAK_bar = p.M_P  / 2;
    M_BG = p.M_BG; M_P = p.M_P;

    deln_pak = linspace(-0.95, 0.95, res) * p.M_P;
    deln_G   = linspace(-0.95, 0.95, res) * p.M_BG;
    [PAK_DIFF, GEF_DIFF] = meshgrid(deln_pak, deln_G);

    % Quarter-domain integral targets per atlas cell (nG_half/2, nP_half/2).
    NGq = (n_GEF_bar + GEF_DIFF/2) / 2;
    NPq = (n_PAK_bar + PAK_DIFF/2) / 2;

    % ---- fresh vs resume -----------------------------------------------------
    if ~isempty(p.quarter_resume_from)
        if ~exist(p.quarter_resume_from, 'file')
            error('run_quarter_hpss_sweeper_for_mpss: quarter_resume_from not found: %s', ...
                  p.quarter_resume_from);
        end
        ck = load(p.quarter_resume_from);
        if ck.res ~= res || ck.M_BG ~= M_BG || ck.M_P ~= M_P || ck.L_total ~= L_total
            error(['run_quarter_hpss_sweeper_for_mpss: checkpoint grid (res=%d,M_BG=%g,M_P=%g,', ...
                   'L=%g) != requested (res=%d,M_BG=%g,M_P=%g,L=%g). Refusing to resume.'], ...
                   ck.res, ck.M_BG, ck.M_P, ck.L_total, res, M_BG, M_P, L_total);
        end
        computation_id = ck.computation_id;
        SOL_grid       = ck.SOL_grid;
        failure_mask   = ck.failure_mask;
        did_phase0to3  = true;   % resume goes straight to the rescue passes
        fprintf('Resuming quarter sweep from %s (id=%s): %d/%d cells already solved.\n', ...
                p.quarter_resume_from, computation_id, sum(~failure_mask(:)), res*res);
    else
        computation_id = datestr(now, 'yyyymmdd_HHMMSS');
        SOL_grid     = cell(res, res);
        failure_mask = true(res, res);   % true == not solved (mark the cell, never fill it in)
        did_phase0to3 = false;
    end

    t_start = tic;
    fprintf('=== run_quarter_hpss_sweeper_for_mpss ===\n');
    fprintf('L_total=%.3f  L_quarter=%.3f  res=%d  use_parallel=%d\n', ...
            L_total, L_quarter, res, p.use_parallel);

    if ~did_phase0to3
        % =====================================================================
        % PHASE 0/1 - CENTER SEED (with re-rooting if the center is off-branch)
        % =====================================================================
        mid = max(1, round(res/2));
        [oi, oj, seed_sol] = find_root_cell(mid, mid, res, NGq, NPq, L_quarter, par, D, p);
        if isempty(seed_sol)
            warning(['run_quarter_hpss_sweeper_for_mpss: no quarter-domain HPSS could be ', ...
                     'cold-started at ANY grid cell (tried up to %d). The quarter branch ', ...
                     'does not exist on this grid; returning all-failed.'], p.max_reroot_tries);
            result = pack_result(SOL_grid, failure_mask, PAK_DIFF, GEF_DIFF, NGq, NPq, ...
                                 L_total, L_half, L_quarter, M_BG, M_P, res, par, ...
                                 computation_id, '', toc(t_start));
            % still write an (empty) output for provenance
            finalize_and_save(result, out_dir);
            return;
        end
        SOL_grid{oi, oj}    = seed_sol;
        failure_mask(oi, oj) = false;
        fprintf('PHASE 1: root cell (%d,%d) solved. u_max=%.3f (center was (%d,%d)).\n', ...
                oi, oj, max(seed_sol.y(1,:)), mid, mid);

        % =====================================================================
        % PHASE 2 - COLUMN CONTINUATION along the root column oj
        % =====================================================================
        fprintf('PHASE 2: column continuation (j=%d)...\n', oj);
        column_seeds = cell(res, 1);
        column_succ  = false(res, 1);
        column_seeds{oi} = seed_sol; column_succ(oi) = true;
        for dir_step = [1, -1]
            if dir_step > 0, i_range = oi+1:res; else, i_range = oi-1:-1:1; end
            prev_sol = seed_sol; prev_prev = [];
            for i = i_range
                [sol_i, ok] = solve_cell_continued(prev_sol, prev_prev, ...
                    NGq(i,oj), NPq(i,oj), L_quarter, par, D, p);
                if ok
                    SOL_grid{i,oj} = sol_i; failure_mask(i,oj) = false;
                    column_seeds{i} = sol_i; column_succ(i) = true;
                    prev_prev = prev_sol; prev_sol = sol_i;
                else
                    column_seeds{i} = prev_sol; column_succ(i) = false;
                end
            end
        end
        fprintf('  -> column solved %d/%d.\n', sum(column_succ), res);

        % =====================================================================
        % PHASE 3 - ROW CONTINUATION from each column seed
        % =====================================================================
        fprintf('PHASE 3: row continuation...\n');
        row_SOL  = cell(res, res);
        row_fail = true(res, res);
        if p.use_parallel
            % parfeval + per-row wall-clock timeout, NOT a plain parfor: a row
            % that repeatedly falls through predictor-corrector -> PDE fallback
            % near a fold can otherwise block the entire phase for hours with
            % no progress visibility (a plain parfor can't be interrupted or
            % polled per-iteration). A timed-out row's unsolved cells simply
            % stay in failure_mask; Phase 4's rescue (already timeout-guarded)
            % picks them up via neighbour/cold start, so nothing is lost.
            pool = gcp;
            futures = parallel.FevalFuture.empty(0, res);
            for i = 1:res
                futures(i) = parfeval(pool, @sweep_one_row, 2, i, oj, column_seeds{i}, ...
                    column_succ(i), NGq, NPq, L_quarter, par, D, p);
            end
            pending = true(1, res);
            n_row_timeout = 0;
            while any(pending)
                live = find(pending);
                [rel, val1, val2] = fetchNext(futures(live), 2);
                if isempty(rel)
                    for li = live
                        % Snapshot StartDateTime ONCE (queued->running is an
                        % unzoned NaT -> zoned-datetime transition; reading it
                        % twice can mix zoned/unzoned and error; the same
                        % guard is used in Phase 4 and in build_mpss_branch.m).
                        sd = futures(li).StartDateTime;
                        if isnat(sd), continue; end
                        tz = sd.TimeZone;
                        if isempty(tz), now_dt = datetime('now');
                        else, now_dt = datetime('now', 'TimeZone', tz); end
                        if seconds(now_dt - sd) > p.phase3_row_timeout
                            cancel(futures(li)); pending(li) = false;
                            n_row_timeout = n_row_timeout + 1;
                            fprintf(['  row %d: TIMEOUT (> %.0fs); cells left failed ', ...
                                     'for Phase 4 rescue.\n'], li, p.phase3_row_timeout);
                        end
                    end
                    continue;
                end
                ii = live(rel); pending(ii) = false;
                row_SOL(ii,:) = val1; row_fail(ii,:) = val2;
            end
            if n_row_timeout > 0
                fprintf('  Phase 3: %d/%d row(s) timed out (> %.0fs); unset cells deferred to Phase 4.\n', ...
                        n_row_timeout, res, p.phase3_row_timeout);
            end
        else
            for i = 1:res
                [row_SOL(i,:), row_fail(i,:)] = sweep_one_row(i, oj, column_seeds{i}, ...
                    column_succ(i), NGq, NPq, L_quarter, par, D, p);
            end
        end
        % Merge row results into the grid (row sweep only overwrites its own row).
        for i = 1:res
            for j = 1:res
                if ~row_fail(i,j)
                    SOL_grid{i,j} = row_SOL{i,j};
                    failure_mask(i,j) = false;
                end
            end
        end
        fprintf('  -> after row sweep: solved %d/%d.\n', sum(~failure_mask(:)), res*res);

        save_quarter_checkpoint(out_dir, computation_id, SOL_grid, failure_mask, ...
            PAK_DIFF, GEF_DIFF, L_total, L_half, L_quarter, M_BG, M_P, res, par, 'post_phase3');
    end

    % =========================================================================
    % PHASE 4 - MULTI-PASS RETRY (neighbour warm-start + cold-start floor)
    % =========================================================================
    % Each remaining failed cell is attempted with:
    %   tier 1: neighbour-seeded continuation (warm start) -- adds coherence,
    %   tier 2: cold start at that cell's quarter density (the 'independent' strategy of
    %           build_mpss_branch.m), so coverage is at least that of independent solves.
    % Newly solved cells become seeds for the next pass. This both fills gaps in
    % the connected component and discovers disconnected components (the full-peak
    % region can have two), which a single center-rooted continuation misses.
    n_before_rescue = sum(~failure_mask(:));
    total_rescued = 0;
    for pass = 1:p.rescue_max_passes
        fail_idx = find(failure_mask);
        if isempty(fail_idx), break; end
        n_fail = numel(fail_idx);

        % pre-slice: neighbour seeds + targets (bounded; no growing broadcast)
        neigh = cell(n_fail,1); es1 = cell(n_fail,1); es2 = cell(n_fail,1);
        nGc = zeros(n_fail,1);  nPc = zeros(n_fail,1);
        for k = 1:n_fail
            [i,j] = ind2sub([res res], fail_idx(k));
            neigh{k} = accepted_neighbor_sols(SOL_grid, failure_mask, i, j);
            [es1{k}, es2{k}] = neighbour_seeds(SOL_grid, failure_mask, i, j);  % collinear pair for extrapolation
            nGc(k) = NGq(i,j); nPc(k) = NPq(i,j);
        end

        new_sols = cell(n_fail,1);
        n_timeout = 0;
        if p.use_parallel
            pool = gcp;
            futures = parallel.FevalFuture.empty(0, n_fail);
            for k = 1:n_fail
                futures(k) = parfeval(pool, @rescue_one_cell, 1, ...
                    neigh{k}, es1{k}, es2{k}, nGc(k), nPc(k), L_quarter, par, D, p);
            end
            pending = true(1, n_fail);
            while any(pending)
                live = find(pending);
                [rel, val] = fetchNext(futures(live), 2);
                if isempty(rel)
                    for li = live
                        % Snapshot StartDateTime ONCE: a future can transition
                        % queued->running (unzoned NaT -> zoned timestamp) between
                        % reads, which would mix a zoned and an unzoned datetime in
                        % the subtraction below and error. Only time out futures that
                        % have actually started.
                        sd = futures(li).StartDateTime;
                        if isnat(sd), continue; end   % not started yet -> no timeout
                        tz = sd.TimeZone;
                        if isempty(tz), now_dt = datetime('now');
                        else, now_dt = datetime('now','TimeZone',tz); end
                        if seconds(now_dt - sd) > p.cell_timeout
                            cancel(futures(li)); pending(li) = false;
                            n_timeout = n_timeout + 1;
                        end
                    end
                    continue;
                end
                li = live(rel); pending(li) = false;
                new_sols{li} = val;
            end
        else
            for k = 1:n_fail
                new_sols{k} = rescue_one_cell(neigh{k}, es1{k}, es2{k}, nGc(k), nPc(k), ...
                                              L_quarter, par, D, p);
            end
        end

        n_pass = 0;
        for k = 1:n_fail
            if ~isempty(new_sols{k})
                [i,j] = ind2sub([res res], fail_idx(k));
                SOL_grid{i,j} = new_sols{k}; failure_mask(i,j) = false;
                n_pass = n_pass + 1;
            end
        end
        total_rescued = total_rescued + n_pass;
        fprintf('  rescue pass %d: %d cells, rescued %d, timed out %d. Remaining: %d.\n', ...
                pass, n_fail, n_pass, n_timeout, sum(failure_mask(:)));
        save_quarter_checkpoint(out_dir, computation_id, SOL_grid, failure_mask, ...
            PAK_DIFF, GEF_DIFF, L_total, L_half, L_quarter, M_BG, M_P, res, par, ...
            sprintf('post_rescue_pass%d', pass));
        if n_pass == 0, break; end
    end
    fprintf('PHASE 4: rescue added %d cell(s) (from %d to %d solved).\n', ...
            total_rescued, n_before_rescue, sum(~failure_mask(:)));

    wall = toc(t_start);
    result = pack_result(SOL_grid, failure_mask, PAK_DIFF, GEF_DIFF, NGq, NPq, ...
                         L_total, L_half, L_quarter, M_BG, M_P, res, par, ...
                         computation_id, '', wall);
    finalize_and_save(result, out_dir);
    fprintf('=== quarter sweep complete: %d/%d quarter-HPSS cells solved. Wall: %.1fs ===\n', ...
            sum(~failure_mask(:)), res*res, wall);
end


% ==========================================================================
%  ORCHESTRATION HELPERS
% ==========================================================================

function [oi, oj, seed_sol] = find_root_cell(ci, cj, res, NGq, NPq, L_q, par, D, p)
    % Cold-start the center cell; if it fails, spiral outward to the nearest cell
    % that DOES cold-start, and use it as the continuation origin. Returns [] seed
    % if no cell on the grid yields a quarter-HPSS (branch absent everywhere).
    order = cells_by_distance(ci, cj, res);
    n_try = min(size(order,1), p.max_reroot_tries);
    for t = 1:n_try
        i = order(t,1); j = order(t,2);
        [sol, ok] = seed_cold_quarter(NGq(i,j), NPq(i,j), L_q, par, D, p);
        if ok
            oi = i; oj = j; seed_sol = sol;
            if t > 1
                fprintf('  (center off-branch; re-rooted to (%d,%d) after %d cold tries)\n', i, j, t);
            end
            return;
        end
    end
    oi = ci; oj = cj; seed_sol = [];
end

function order = cells_by_distance(ci, cj, res)
    % (i,j) pairs sorted by Euclidean distance from (ci,cj), nearest first.
    [J, I] = meshgrid(1:res, 1:res);
    d = (I-ci).^2 + (J-cj).^2;
    [~, idx] = sort(d(:));
    order = [I(idx), J(idx)];
end

function [row_SOL, row_fail] = sweep_one_row(i, oj, seed, seed_ok, NGq, NPq, L_q, par, D, p)
    % Predictor-corrector across columns of row i, starting from its column seed
    % at (i,oj) and marching outward both directions. Self-contained so it is
    % safe as a parfor body (no shared-state writes).
    res = size(NGq, 2);
    row_SOL  = cell(1, res);
    row_fail = true(1, res);
    if ~seed_ok || isempty(seed), return; end
    row_SOL{oj} = seed; row_fail(oj) = false;
    for dir_step = [1, -1]
        if dir_step > 0, j_range = oj+1:res; else, j_range = oj-1:-1:1; end
        prev_sol = seed; prev_prev = [];
        for j = j_range
            [sol_j, ok] = solve_cell_continued(prev_sol, prev_prev, ...
                NGq(i,j), NPq(i,j), L_q, par, D, p);
            if ok
                row_SOL{j} = sol_j; row_fail(j) = false;
                prev_prev = prev_sol; prev_sol = sol_j;
            end
        end
    end
end

function [sol_out, ok] = solve_cell_continued(prev_sol, prev_prev, nG, nP, L_q, par, D, p)
    % One continuation step: predictor-corrector BVP, then a PDE-extrapolation
    % fallback (same two-tier idea as run_hybrid_sweeper_pdefb Phase 2). Accepts
    % only a validated half-peak.
    sol_out = []; ok = false;
    guess = predictor_guess(prev_sol, prev_prev);
    [bvp, okb] = run_bvp_v7(guess, L_q, par, D, nG, nP);
    if okb && validate_halfwave(bvp) && quarter_mass_ok(bvp, nG, nP, p.quarter_mass_tol)
        sol_out = force_peak_at_zero(bvp); ok = true; return;
    end
    bvp = pde_extrap_fallback(prev_sol, prev_prev, nG, nP, L_q, par, D, p);
    if ~isempty(bvp) && quarter_mass_ok(bvp, nG, nP, p.quarter_mass_tol)
        sol_out = force_peak_at_zero(bvp); ok = true;
    end
end

function sol_out = rescue_one_cell(neigh_sols, seed1, seed2, nG, nP, L_q, par, D, p)
    % Three-tier retry for a still-failed cell:
    %   tier 1: try EACH accepted cardinal-neighbour profile as a full-NMax
    %           run_bvp_v7 warm start (validated + mass-gated). Trying every
    %           accepted neighbour -- not just one, at the full NMax the core uses
    %           -- is what fills the interior holes (it recovered about 85% of the holes
    %           that a single tight-NMax retry left behind).
    %   tier 2: collinear-pair linear extrapolation (gap-bridging). A warm start
    %           sits AT a neighbour's mass; extrapolation predicts one step PAST
    %           the solved frontier, reaching boundary-frontier cells a warm start
    %           cannot (it recovered about 75% of the frontier cells that a
    %           warm-start-only retry left behind).
    %   tier 3: cold start at this density (coverage floor = the independent strategy).
    sol_out = [];
    for s = 1:numel(neigh_sols)
        [bvp, ok] = run_bvp_v7(construct_bvp_guess_v7(neigh_sols{s}), L_q, par, D, nG, nP);
        if ok && validate_halfwave(bvp) && quarter_mass_ok(bvp, nG, nP, p.quarter_mass_tol)
            sol_out = force_peak_at_zero(bvp); return;
        end
    end
    if ~isempty(seed1)
        bvp = bvp_seed_extrap(seed1, seed2, nG, nP, L_q, par, D);
        if ~isempty(bvp) && quarter_mass_ok(bvp, nG, nP, p.quarter_mass_tol)
            sol_out = force_peak_at_zero(bvp); return;
        end
    end
    [bvp, ok] = seed_cold_quarter(nG, nP, L_q, par, D, p);
    if ok, sol_out = bvp; end   % seed_cold_quarter already forces peak-at-zero + mass gate
end

function sols = accepted_neighbor_sols(sol_grid, fail_mask, i, j)
    % All accepted (solved, non-empty) cardinal-neighbour solutions of (i,j).
    % Bounded (<=4) so the parfeval broadcast stays small (no growing
    % accumulator).
    sols = {};
    [ni, nj] = size(fail_mask);
    cand = {[i-1 j],[i+1 j],[i j-1],[i j+1]};
    for k = 1:numel(cand)
        ii = cand{k}(1); jj = cand{k}(2);
        if ii>=1 && ii<=ni && jj>=1 && jj<=nj && ~fail_mask(ii,jj) && ~isempty(sol_grid{ii,jj})
            sols{end+1} = sol_grid{ii,jj}; %#ok<AGROW>
        end
    end
end


% ==========================================================================
%  PER-CELL SOLVERS (call the shared standalone BVP core only)
% ==========================================================================

function [sol_q, ok] = seed_cold_quarter(nG_q, nP_q, L_q, par, D, p)
    % Cold start: HSS -> cos-perturbed PDE relax -> BVP -> validate -> force peak
    % at x=0. Same recipe as quarter_hpss_solve in build_mpss_branch.m.
    %
    % When p.osc_fallback is true the relaxation is oscillation-aware: it bounds
    % the run (early-stops once it detects the cos-perturbed trajectory is not
    % relaxing) and offers the time-mean + phase snapshots as bvp4c SEED
    % candidates, so a quarter regime that oscillates (rather than settling to a
    % steady quarter peak) can still recover a validated steady seed. Snapshots
    % are seeds only; validate_halfwave and the quarter mass gate still apply.
    sol_q = []; ok = false;
    try
        [u, q, v, w, pv] = get_hss_local(nG_q, nP_q, L_q, par);
    catch
        return;
    end
    if any(~isfinite([u q v w pv])), return; end
    x = linspace(0, L_q, 120);

    if p.osc_fallback
        % Same cos-perturbed HSS IC as run_pde_chunked, then the oscillation-
        % aware relax -> candidates -> first candidate that yields a validated,
        % on-mass bvp4c solve (validate_halfwave enforces amplitude >= 0.5 on
        % the SOLUTION, so no separate amplitude pre-check of the relaxed profile is needed).
        y0 = cos_perturbed_ic([u; q; v; w; pv], x, L_q);
        cands = pde_relax_candidates_q(y0, x, par, D, ...
                    p.pde_T_end, p.pde_T_chunk, p.pde_max_chunks, p);
        bvp = seed_bvp_from_candidates_q(cands, L_q, par, D, nG_q, nP_q, ...
                    @(b) quarter_mass_ok(b, nG_q, nP_q, p.quarter_mass_tol));
        if isempty(bvp), return; end
        sol_q = force_peak_at_zero(bvp); ok = true;
        return;
    end

    pde = run_pde_chunked([u; q; v; w; pv], x, L_q, par, D, ...
                          p.pde_T_end, p.pde_T_chunk, p.pde_max_chunks);
    if isempty(pde), return; end
    if (max(pde.y(1,:)) - min(pde.y(1,:))) < 0.5, return; end
    [sol, okb] = run_bvp_v7(construct_bvp_guess_v7(pde), L_q, par, D, nG_q, nP_q);
    if ~okb, return; end
    if ~validate_halfwave(sol), return; end
    if ~quarter_mass_ok(sol, nG_q, nP_q, p.quarter_mass_tol), return; end
    sol_q = force_peak_at_zero(sol); ok = true;
end

function y0 = cos_perturbed_ic(y_hss, x, L)
    % Mass-conserving cos perturbation of a uniform HSS IC -- same as
    % run_pde_chunked's (0.5*u*cos(pi*x/L) added to u, subtracted from p), so
    % the osc-aware cold start seeds from the same symmetry-broken state the
    % non-osc path does.
    N = length(x);
    if size(y_hss, 2) == 1, y0 = repmat(y_hss, 1, N); else, y0 = y_hss; end
    u_ss = y0(1,:); p_ss = y0(5,:);
    perturb = 0.5 * u_ss .* cos(pi * x / L);
    y0(1,:) = u_ss + perturb;
    y0(5,:) = p_ss - perturb;
end

function ok = quarter_mass_ok(sol, nG_q, nP_q, tol)
    ok = quarter_mass_res(sol, nG_q, nP_q) <= tol;
end

function mres = quarter_mass_res(sol, nG_q, nP_q)
    % MESH-ACCURATE quarter mass residual. The 14-var BVP enforces the integral
    % BC via the s1,s2 multipliers (rows 11,12), so recomputing the mass is right
    % (s1,s2 can absorb an off-target mass) -- BUT a plain trapz on bvp4c's
    % NON-UNIFORM mesh conflates mesh-discretization error (empirically ~1e-4)
    % with true mass error. Genuinely on-mass quarter solves
    % (|s1,s2|~1e-17, deval mass ~1e-13) show coarse-trapz residuals up to 3.5e-4,
    % which a gate on the coarse value would wrongly reject. pchip-densifying the species onto a
    % uniform fine mesh before trapz measures ~1e-6 -- an accuracy fix, not a
    % tolerance change. Species rows: u=1,q=3,v=5,w=7,p=9.
    x = sol.x(:).';
    xf = linspace(x(1), x(end), 2000);
    u = interp1(x, sol.y(1,:), xf, 'pchip');
    q = interp1(x, sol.y(3,:), xf, 'pchip');
    v = interp1(x, sol.y(5,:), xf, 'pchip');
    w = interp1(x, sol.y(7,:), xf, 'pchip');
    pp = interp1(x, sol.y(9,:), xf, 'pchip');
    gef = trapz(xf, u + q + v + w);
    pak = trapz(xf, pp + q);
    mres = max(abs(gef - nG_q), abs(pak - nP_q));
end

function sol = force_peak_at_zero(sol)
    % Ensure u(0) >= u(end) so the quarter peak-end maps to the full-peak seam. If we
    % must flip, mirror the 5 species in space and drop the derivative rows
    % (reflect_stitch then finite-differences the seam slope) -- same as the
    % orientation guard of quarter_hpss_solve in build_mpss_branch.m. Usually a no-op because the
    % cos IC and continuation both keep the peak at x=0.
    if sol.y(1,1) < sol.y(1,end)
        yr = zeros(size(sol.y));
        yr([1,3,5,7,9],:) = fliplr(sol.y([1,3,5,7,9],:));
        sol.y = yr;
    end
end

function guess = predictor_guess(prev_sol, prev_prev)
    % 2-point linear extrapolation of the profile (positivity-guarded); falls back
    % to the last solution when only one is available. All cells share L=L_quarter
    % so meshes align by resampling. (Same as predictor_step in run_hybrid_sweeper_pdefb.m.)
    if isempty(prev_prev)
        guess = construct_bvp_guess_v7(prev_sol); return;
    end
    x_ref = prev_sol.x; y_prev = prev_sol.y;
    y_pp = zeros(size(y_prev));
    for k = 1:size(y_prev,1)
        y_pp(k,:) = interp1(prev_prev.x, prev_prev.y(k,:), x_ref, 'pchip', 'extrap');
    end
    y_pred = y_prev + (y_prev - y_pp);
    y_pred([1,3,5,7,9],:) = max(0, y_pred([1,3,5,7,9],:));
    pred.x = x_ref; pred.y = y_pred;
    guess = construct_bvp_guess_v7(pred);
end

function sol_out = pde_extrap_fallback(prev_sol, prev_prev, nG, nP, L_q, par, D, p)
    % Extrapolate the neighbour profile one mass step, relax under the RD PDE
    % (no symmetry-breaking cos -- the IC is already polarized), then BVP. Same as
    % pde_fallback_extrapolated in run_hybrid_sweeper_pdefb.m, on the quarter domain.
    sol_out = [];
    if isempty(prev_sol), return; end
    x_ic = prev_sol.x; y_ic = prev_sol.y;
    if ~isempty(prev_prev)
        y_pp = zeros(size(y_ic));
        for k = 1:size(y_ic,1)
            y_pp(k,:) = interp1(prev_prev.x, prev_prev.y(k,:), x_ic, 'pchip', 'extrap');
        end
        y_ic = y_ic + (y_ic - y_pp);
    end
    y_ic([1,3,5,7,9],:) = max(0, y_ic([1,3,5,7,9],:));
    N = 100; x_pde = linspace(0, L_q, N); rows = [1,3,5,7,9];
    y0 = zeros(5, N);
    for k = 1:5, y0(k,:) = interp1(x_ic, y_ic(rows(k),:), x_pde, 'pchip'); end
    y0 = max(0, y0);
    max_ch = max(2, ceil(p.pde_max_chunks/2));

    if p.osc_fallback
        % Bounded oscillation-aware relaxation -> seed candidates, each tried
        % through the SAME BVP core + validate_halfwave + quarter mass gate.
        cands = pde_relax_candidates_q(y0, x_pde, par, D, ...
                    p.pde_fallback_T, p.pde_T_chunk, max_ch, p);
        sol_out = seed_bvp_from_candidates_q(cands, L_q, par, D, nG, nP, ...
                    @(b) quarter_mass_ok(b, nG, nP, p.quarter_mass_tol));
        return;
    end

    pde_sol = run_pde_polarised_local(y0, x_pde, L_q, par, D, ...
                  p.pde_fallback_T, p.pde_T_chunk, max_ch);
    if isempty(pde_sol), return; end
    [bvp, ok] = run_bvp_v7(construct_bvp_guess_v7(pde_sol), L_q, par, D, nG, nP);
    if ~ok, return; end
    if ~validate_halfwave(bvp), return; end
    sol_out = bvp;
end


function cands = pde_relax_candidates_q(y0, x, par, D, T_end, T_chunk, max_chunks, p)
    % Quarter-domain twin of run_hybrid_sweeper_pdefb.m's pde_relax_candidates.
    % Bounded, oscillation-aware relaxation returning ORDERED seed-candidate
    % structs (.x, 14-row .y). Converged -> the single steady state (the single-shot seed);
    % not relaxing -> time-mean then phase snapshots. Seeds only, no acceptance.
    cands = {};
    if T_end <= 0, return; end
    if T_chunk > T_end, T_chunk = T_end; end
    if T_chunk <= 0, T_chunk = T_end/2; end
    if max_chunks < 1, max_chunks = 1; end
    if size(y0,2) == 1, y0 = repmat(y0, 1, length(x)); end
    nsnap = max(2, p.osc_n_snap);

    curr_y = y0; deltas = []; last_frames = {}; converged = false;
    ws = warning('off','all'); cleanup = onCleanup(@() warning(ws)); %#ok<NASGU>
    try
        for k = 1:max_chunks
            tspan = linspace(0, T_chunk, nsnap + 1);
            sol = pdepe(0, @(xq,t,uu,dudx) pde_eq_v7(xq,t,uu,dudx,par,D), ...
                        @(xq) pde_ic_local(xq, curr_y, x), @pde_bc_v7, x, tspan);
            last_frames = cell(1, nsnap);
            for s = 1:nsnap, last_frames{s} = squeeze(sol(s+1,:,:))'; end
            frame_end = last_frames{end};
            d = max(abs(frame_end(1,:) - curr_y(1,:)));
            deltas(end+1) = d; %#ok<AGROW>
            curr_y = frame_end;
            if d < p.osc_conv_tol, converged = true; break; end
            if k >= p.osc_min_chunks && ...
               min(deltas(max(1,end-p.osc_min_chunks+1):end)) > p.osc_rel * deltas(1)
                break;
            end
        end
    catch
        cands = {}; return;
    end

    if converged
        cands = { state5_to_seed_q(curr_y, x) }; return;
    end
    if isempty(last_frames), return; end
    ymean = zeros(size(last_frames{1}));
    for s = 1:numel(last_frames), ymean = ymean + last_frames{s}; end
    ymean = max(0, ymean / numel(last_frames));
    cands{end+1} = state5_to_seed_q(ymean, x);
    for s = 1:numel(last_frames)
        cands{end+1} = state5_to_seed_q(max(0, last_frames{s}), x); %#ok<AGROW>
    end
end

function ss = state5_to_seed_q(y5, x)
    ss.x = x;
    yf = zeros(14, size(y5, 2));
    yf(1,:) = y5(1,:); yf(3,:) = y5(2,:); yf(5,:) = y5(3,:);
    yf(7,:) = y5(4,:); yf(9,:) = y5(5,:);
    ss.y = yf;
end

function bvp_out = seed_bvp_from_candidates_q(cands, L, par, D, nG, nP, mass_gate)
    % Try each candidate; accept the FIRST bvp4c solve that converges AND passes
    % validate_halfwave AND the mass gate. [] if none. (force_peak_at_zero is
    % applied by the caller, exactly as for the single-shot path.)
    bvp_out = [];
    for c = 1:numel(cands)
        [bvp, ok] = run_bvp_v7(construct_bvp_guess_v7(cands{c}), L, par, D, nG, nP);
        if ok && validate_halfwave(bvp) && mass_gate(bvp)
            bvp_out = bvp; return;
        end
    end
end

function sol_out = bvp_seed_extrap(prev_sol, prev_prev, nG, nP, L_q, par, D)
    % Fast neighbour-seeded BVP with a tight NMax so a stiff input can't run the
    % mesh away (bounds Phase-4 per-cell wall time). Same as
    % bvp_seed_extrapolated in run_hybrid_sweeper_pdefb.m, on the quarter domain.
    sol_out = [];
    if isempty(prev_sol), return; end
    x_ic = prev_sol.x; y_ic = prev_sol.y;
    if ~isempty(prev_prev)
        y_pp = zeros(size(y_ic));
        for k = 1:size(y_ic,1)
            y_pp(k,:) = interp1(prev_prev.x, prev_prev.y(k,:), x_ic, 'pchip', 'extrap');
        end
        y_ic = y_ic + (y_ic - y_pp);
    end
    y_ic([1,3,5,7,9],:) = max(0, y_ic([1,3,5,7,9],:));
    guess_in.x = x_ic; guess_in.y = y_ic;
    guess = construct_bvp_guess_v7(guess_in);
    options = bvpset('RelTol', 1e-4, 'NMax', 1500);
    ws = warning('off','all'); cleanup = onCleanup(@() warning(ws)); %#ok<NASGU>
    try
        bvp = bvp4c(@(x,y) ode_sys_v7(x,y,par,D), ...
                    @(ya,yb) bc_aug_v7(ya,yb,nG,nP), guess, options);
    catch
        return;
    end
    if isempty(bvp) || numel(bvp.x) >= 1500, return; end
    if any(isnan(bvp.y(:))) || any(isinf(bvp.y(:))), return; end
    if min(bvp.y([1,3,5,7,9],:), [], 'all') < -1e-3, return; end
    if ~validate_halfwave(bvp), return; end
    sol_out = bvp;
end

function sol_out = run_pde_polarised_local(y0, x, L, par, D, T_end, T_chunk, max_chunks)
    % Relax an already-polarized IC under the RD PDE with reflecting BCs, NO
    % cos perturbation. Same as run_pde_polarised in run_hybrid_sweeper_pdefb.m (a local function there),
    % reproduced here so this sweeper stays self-contained.
    sol_out = [];
    if T_end <= 0, return; end
    if T_chunk > T_end, T_chunk = T_end; end
    if T_chunk <= 0, T_chunk = T_end/2; end
    if max_chunks < 1, max_chunks = 1; end
    if size(y0,2) == 1, y0 = repmat(y0, 1, length(x)); end
    curr_y = y0;
    ws = warning('off','all'); cleanup = onCleanup(@() warning(ws)); %#ok<NASGU>
    try
        for k = 1:max_chunks
            tspan = [0, T_chunk/2, T_chunk];
            sol = pdepe(0, @(xq,t,uu,dudx) pde_eq_v7(xq,t,uu,dudx,par,D), ...
                        @(xq) pde_ic_local(xq, curr_y, x), @pde_bc_v7, x, tspan);
            frame_end = squeeze(sol(end,:,:))';
            if max(abs(frame_end(1,:) - curr_y(1,:))) < 1e-6
                curr_y = frame_end; break;
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

function [seed1, seed2] = neighbour_seeds(sol_grid, fail_mask, i, j)
    % Prefer a collinear polarized pair (for linear extrapolation toward (i,j));
    % else the closest single polarized cardinal neighbour. Flat (no-peak) neighbours
    % are skipped. Same as find_neighbour_seeds in run_hybrid_sweeper_pdefb.m.
    [ni, nj] = size(fail_mask);
    seed1 = []; seed2 = [];
    if i >= 3 && is_pol(sol_grid,fail_mask,i-1,j) && is_pol(sol_grid,fail_mask,i-2,j)
        seed1 = sol_grid{i-1,j}; seed2 = sol_grid{i-2,j}; return; end
    if i <= ni-2 && is_pol(sol_grid,fail_mask,i+1,j) && is_pol(sol_grid,fail_mask,i+2,j)
        seed1 = sol_grid{i+1,j}; seed2 = sol_grid{i+2,j}; return; end
    if j >= 3 && is_pol(sol_grid,fail_mask,i,j-1) && is_pol(sol_grid,fail_mask,i,j-2)
        seed1 = sol_grid{i,j-1}; seed2 = sol_grid{i,j-2}; return; end
    if j <= nj-2 && is_pol(sol_grid,fail_mask,i,j+1) && is_pol(sol_grid,fail_mask,i,j+2)
        seed1 = sol_grid{i,j+1}; seed2 = sol_grid{i,j+2}; return; end
    cand = {[i-1,j],[i+1,j],[i,j-1],[i,j+1]};
    for k = 1:numel(cand)
        ii = cand{k}(1); jj = cand{k}(2);
        if ii>=1 && ii<=ni && jj>=1 && jj<=nj && is_pol(sol_grid,fail_mask,ii,jj)
            seed1 = sol_grid{ii,jj}; return;
        end
    end
end

function tf = is_pol(sol_grid, fail_mask, i, j)
    tf = false;
    if i<1 || j<1 || i>size(fail_mask,1) || j>size(fail_mask,2), return; end
    if fail_mask(i,j), return; end
    s = sol_grid{i,j};
    if isempty(s) || ~isfield(s,'y'), return; end
    tf = (max(s.y(1,:)) - min(s.y(1,:))) > 1e-3;
end


% ==========================================================================
%  CHECKPOINT / SAVE
% ==========================================================================

function ckpt_path = save_quarter_checkpoint(out_dir, computation_id, SOL_grid, failure_mask, ...
        PAK_DIFF, GEF_DIFF, L_total, L_half, L_quarter, M_BG, M_P, res, par, tag)
    ckpt_path = fullfile(out_dir, sprintf('checkpoint_quarter_hpss_%s.mat', computation_id));
    n_solved = sum(~failure_mask(:));
    try
        save(ckpt_path, 'SOL_grid', 'failure_mask', 'PAK_DIFF', 'GEF_DIFF', ...
             'L_total', 'L_half', 'L_quarter', 'M_BG', 'M_P', 'res', 'par', ...
             'computation_id', 'tag', 'n_solved');
        fprintf('  [quarter checkpoint: %s, tag=%s, solved=%d/%d]\n', ...
                ckpt_path, tag, n_solved, res*res);
    catch ME
        warning('run_quarter_hpss_sweeper_for_mpss: checkpoint save failed (%s): %s', tag, ME.message);
    end
end

function result = pack_result(SOL_grid, failure_mask, PAK_DIFF, GEF_DIFF, NGq, NPq, ...
        L_total, L_half, L_quarter, M_BG, M_P, res, par, computation_id, out_path, wall)
    result = struct();
    result.SOL_grid = SOL_grid;
    result.failure_mask = failure_mask;
    result.PAK_DIFF = PAK_DIFF;
    result.GEF_DIFF = GEF_DIFF;
    result.NGq = NGq; result.NPq = NPq;
    result.L_total = L_total; result.L_half = L_half; result.L_quarter = L_quarter;
    result.M_BG = M_BG; result.M_P = M_P; result.res = res;
    result.par = par; result.computation_id = computation_id;
    result.out_path = out_path; result.wall_time = wall;
    result.n_solved = sum(~failure_mask(:));
    result.n_total = res*res;
end

function finalize_and_save(result, out_dir)
    % Write a timestamped quarter-grid file for provenance (not a library file).
    out_path = fullfile(out_dir, sprintf('quarter_hpss_for_mpss_%s.mat', result.computation_id));
    SOL_grid = result.SOL_grid; failure_mask = result.failure_mask; %#ok<NASGU>
    PAK_DIFF = result.PAK_DIFF; GEF_DIFF = result.GEF_DIFF; %#ok<NASGU>
    NGq = result.NGq; NPq = result.NPq; %#ok<NASGU>
    L_total = result.L_total; L_half = result.L_half; L_quarter = result.L_quarter; %#ok<NASGU>
    M_BG = result.M_BG; M_P = result.M_P; res = result.res; %#ok<NASGU>
    par = result.par; computation_id = result.computation_id; %#ok<NASGU>
    try
        save(out_path, 'SOL_grid', 'failure_mask', 'PAK_DIFF', 'GEF_DIFF', 'NGq', 'NPq', ...
             'L_total', 'L_half', 'L_quarter', 'M_BG', 'M_P', 'res', 'par', 'computation_id');
        fprintf('  wrote %s\n', out_path);
    catch ME
        warning('run_quarter_hpss_sweeper_for_mpss: final save failed: %s', ME.message);
    end
end
