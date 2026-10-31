function report = build_mpss_branch(p, out_dir)
    % BUILD_MPSS_BRANCH  MPSS branch-library constructor (Packet 1/1B/1C).
    %
    % MPSS = quarter-domain HPSS solve + peak-to-peak reflection (NOT a direct
    % interior-peak BVP). See docs/ai_handoff/mpss_constructor_spec.md for the
    % fixed design; the numerical convention there (L_quarter=L_total/4, quarter
    % masses = half/2, peak-to-peak reflection, flux at the outer trough on the
    % half-domain interface) is UNCHANGED by this pass. Failed quarter solves /
    % failed reflection checks are MASKED, never fabricated (Invariant 2).
    % Stability is annotation only and never gates existence (Invariant 3).
    %
    % --- 2026-07-09 Packet 1C: production-scale batching/checkpointing --------
    % Adds, around the UNCHANGED per-cell math (quarter_hpss_solve,
    % reflect_stitch, classify_mpss_profile, and the shared BVP core):
    %  * PARALLEL execution: each pending cell is an independent parfeval task
    %    (cells have no cross-cell dependency, unlike the sweeper's neighbour-
    %    seeded HPSS rescue -- so unlike that rescue's multi-pass design, MPSS
    %    submits one wave over all pending cells and streams results back).
    %  * per-cell TIMEOUT via fetchNext + cancel (same pattern as
    %    run_hybrid_sweeper_pdefb.m Phase 4; Invariant 8). Timed-out cells are
    %    masked with fail_reason="timeout", never fabricated.
    %  * CHECKPOINTING every p.checkpoint_every resolved cells, to
    %    checkpoint_mpss_<computation_id>.mat (Invariant 8).
    %  * RESUME via p.resume_from = a checkpoint path: already-resolved cells
    %    (existence_mask true OR a definitive fail_reason) are never
    %    resubmitted; the SAME computation_id is continued.
    %  * BATCHING via p.max_cells_this_run: bound how many pending cells one
    %    call attempts, so a run can be deliberately split across calls.
    %  * SAFE OUTPUT: branch_MPSS_<id>.mat is always written with whatever was
    %    resolved (honest partial record); branch_MPSS_latest.mat is updated
    %    only if the run is complete (every cell at least attempted) or the
    %    caller explicitly sets p.accept_partial=true.
    % Nothing in quarter_hpss_solve, reflect_stitch, classify_mpss_profile
    % (acceptance thresholds/order), or the shared BVP core was changed.
    %
    % --- SHARED BVP CORE (Packet 1B refactor, 2026-07-09) --------------------
    % The per-point BVP core (run_bvp_v7, ode_sys_v7, bc_aug_v7,
    % construct_bvp_guess_v7, run_pde_chunked, pde_eq_v7, pde_bc_v7, pde_ic_local,
    % get_hss_local, hss_eqs, validate_halfwave) lives as STANDALONE files in
    % src/solvers/, extracted verbatim from run_hybrid_sweeper_pdefb.m. Both this
    % constructor and the sweeper call those single-source-of-truth files.
    %
    % Output: data/outputs/branch_MPSS_<computation_id>.mat (+ branch_MPSS_latest.mat
    % when the run is complete or p.accept_partial=true), plus
    % data/outputs/checkpoint_mpss_<computation_id>.mat. Never writes to or
    % modifies any canonical source .mat.

    % ---- parameters (default: V7 Golden Standard, L_total=4 -> L_quarter=1) --
    if nargin < 1 || isempty(p)
        p.a = 0.59; p.b = 0.13; p.c = 21.0;
        p.d = 1.3;  p.e = 33.0; p.f = 1.9;
        p.Dm = 0.001; p.Dc = 1.0;
        p.L = 4.0; p.M_BG = 16.0; p.M_P = 8.0;
        p.res = 3;
    end
    if nargin < 2 || isempty(out_dir), out_dir = 'data/outputs'; end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    if ~isfield(p, 'res'),            p.res = 3;      end
    if ~isfield(p, 'pde_T_end'),      p.pde_T_end   = 2e6;  end
    if ~isfield(p, 'pde_T_chunk'),    p.pde_T_chunk = 1e6;  end
    if ~isfield(p, 'pde_max_chunks'), p.pde_max_chunks = 12; end
    if ~isfield(p, 'slope_tol'),      p.slope_tol = 1e-2;   end  % C1 seam-slope tol
    if ~isfield(p, 'sym_tol'),        p.sym_tol   = 1e-6;   end  % symmetry tol
    if ~isfield(p, 'mass_tol'),       p.mass_tol  = 1e-4;   end  % mass-doubling tol
    if ~isfield(p, 'max_res_warn'),   p.max_res_warn = 8;   end  % serial-mode guideline
    % --- Packet 1C options ---
    if ~isfield(p, 'use_parallel'),      p.use_parallel = true;  end
    if ~isfield(p, 'checkpoint_every'),  p.checkpoint_every = 5; end  % cells per checkpoint
    if ~isfield(p, 'cell_timeout'),      p.cell_timeout = 60;    end  % seconds, parallel only
    if ~isfield(p, 'resume_from'),       p.resume_from = '';     end  % checkpoint_mpss_*.mat path
    if ~isfield(p, 'accept_partial'),    p.accept_partial = false; end
    if ~isfield(p, 'max_cells_this_run'), p.max_cells_this_run = Inf; end
    % --- Packet 1D: search strategy ---
    %   'continuation' (default) -> run_quarter_hpss_sweeper_for_mpss (center
    %       seed + predictor-corrector column/row continuation + neighbour rescue),
    %       then reflect/stitch each solved quarter-HPSS profile into an MPSS.
    %   'independent'            -> legacy per-cell cold solve (Packet 1C path).
    if ~isfield(p, 'mpss_strategy'),       p.mpss_strategy = 'continuation'; end
    if ~isfield(p, 'quarter_resume_from'), p.quarter_resume_from = '';       end

    if ~p.use_parallel && p.res > p.max_res_warn
        warning(['build_mpss_branch: res=%d exceeds the serial-mode guideline (%d) with ', ...
                 'use_parallel=false. Serial mode has no per-cell timeout and may be slow ', ...
                 'at this resolution -- see project_invariants.md Invariant 8, or set ', ...
                 'p.use_parallel=true (default) for parfeval + timeout + checkpointing.'], ...
                p.res, p.max_res_warn);
    end

    par = p;
    D = [p.Dm; p.Dm; p.Dc; p.Dc; p.Dc];
    res = p.res;
    L_total   = p.L;
    L_half    = L_total / 2;
    L_quarter = L_total / 4;
    n_GEF_bar = p.M_BG / 2;
    n_PAK_bar = p.M_P  / 2;

    deln_pak = linspace(-0.95, 0.95, res) * p.M_P;
    deln_G   = linspace(-0.95, 0.95, res) * p.M_BG;
    [PAK_DIFF, GEF_DIFF] = meshgrid(deln_pak, deln_G);

    % Packet 1D: continuation mode manages its own quarter-sweep checkpoint/resume
    % (p.quarter_resume_from); the per-cell MPSS resume path below is only for the
    % legacy independent strategy. Force a fresh MPSS array init in continuation.
    if strcmpi(p.mpss_strategy, 'continuation') && ~isempty(p.resume_from)
        warning(['build_mpss_branch: p.resume_from is ignored in continuation mode. ', ...
                 'Use p.quarter_resume_from to resume the quarter sweep.']);
        p.resume_from = '';
    end

    % ---- fresh start vs. resume from an existing checkpoint ------------------
    if ~isempty(p.resume_from)
        if ~exist(p.resume_from, 'file')
            error('build_mpss_branch: resume_from file not found: %s', p.resume_from);
        end
        ck = load(p.resume_from);
        if ck.res ~= res || ck.M_BG ~= p.M_BG || ck.M_P ~= p.M_P || ck.L_total ~= p.L
            error(['build_mpss_branch: checkpoint grid (res=%d,M_BG=%g,M_P=%g,L=%g) does not ', ...
                   'match the requested grid (res=%d,M_BG=%g,M_P=%g,L=%g). Refusing to resume ', ...
                   '-- grids must match exactly (Invariant 5/11).'], ...
                   ck.res, ck.M_BG, ck.M_P, ck.L_total, res, p.M_BG, p.M_P, p.L);
        end
        computation_id = ck.computation_id;   % continue under the SAME id
        existence_mask = ck.existence_mask;
        SOL            = ck.SOL;
        peak_x         = ck.peak_x;
        amp            = ck.amp;
        mass_residual  = ck.mass_residual;
        seam_slope     = ck.seam_slope;
        sym_residual   = ck.sym_residual;
        flux_values    = ck.flux_values;
        flux_source_location = ck.flux_source_location;
        fail_reason    = ck.fail_reason;
        n_already = sum(existence_mask(:) | (fail_reason(:) ~= ""));
        fprintf('Resuming from %s (computation_id=%s): %d/%d cells already resolved.\n', ...
                p.resume_from, computation_id, n_already, res*res);
    else
        computation_id = datestr(now, 'yyyymmdd_HHMMSS');
        existence_mask = false(res, res);
        SOL            = cell(res, res);
        peak_x         = nan(res, res);
        amp            = nan(res, res);
        mass_residual  = nan(res, res);
        seam_slope     = nan(res, res);   % diagnostic: |du/dx| at the reflection seam
        sym_residual   = nan(res, res);   % diagnostic: symmetry error about L_half/2
        flux_values    = nan(res, res, 5);
        flux_source_location = strings(res, res);
        fail_reason    = strings(res, res);   % '' where never attempted
    end

    t_start = tic;
    fprintf('=== build_mpss_branch (Packet 1C/1D) ===\n');
    fprintf('L_total=%.3f  L_half=%.3f  L_quarter=%.3f  res=%d  strategy=%s  use_parallel=%d\n', ...
            L_total, L_half, L_quarter, res, p.mpss_strategy, p.use_parallel);

    n_timed_out_this_run = 0;
    n_resolved_since_checkpoint = 0;
    ckpt_args = {out_dir, computation_id};   % shared trailing context for checkpoint calls
    quarter_source = '';

    if strcmpi(p.mpss_strategy, 'continuation')
        % ================================================================
        % Packet 1D: continuation-based quarter-domain HPSS sweep, then
        % reflect/stitch/classify each solved quarter profile into an MPSS.
        % run_quarter_hpss_sweeper_for_mpss replaces the fragmented per-cell
        % cold-solve search. The reflect_stitch + classify_mpss_profile math
        % (acceptance thresholds/order) is UNCHANGED; only the *source* of the
        % quarter profiles changed: independent cold solves -> continuation.
        % ================================================================
        qres = run_quarter_hpss_sweeper_for_mpss(p, out_dir);
        quarter_source = sprintf('quarter_hpss_for_mpss_%s.mat', qres.computation_id);
        for lin = 1:res*res
            [ii, jj] = ind2sub([res res], lin);
            nG_half = n_GEF_bar + GEF_DIFF(ii,jj)/2;
            nP_half = n_PAK_bar + PAK_DIFF(ii,jj)/2;
            if qres.failure_mask(ii,jj)
                fail_reason(ii,jj) = "quarter_sweep_failed";
                continue;
            end
            sol_q = qres.SOL_grid{ii,jj};
            [xm, ym, seam_du] = reflect_stitch(sol_q);
            [valid, why2, this_amp, pk_x, mass_res, sym_res] = ...
                classify_mpss_profile(xm, ym, seam_du, nG_half, nP_half, p);
            if ~valid
                fail_reason(ii,jj) = why2;
                continue;
            end
            existence_mask(ii,jj) = true;
            SOL{ii,jj}.x = xm; SOL{ii,jj}.y = ym;
            peak_x(ii,jj) = pk_x; amp(ii,jj) = this_amp;
            mass_residual(ii,jj) = mass_res; seam_slope(ii,jj) = seam_du;
            sym_residual(ii,jj) = sym_res;
            flux_values(ii,jj,:) = ym(:, end)';   % outer trough at half-domain interface
            flux_source_location(ii,jj) = "outer_trough_at_half_domain_interface";
            log_cell_result(ii, jj, true, this_amp, pk_x, mass_res, seam_du, sym_res, "");
        end

    else
    % Only cells with NO recorded outcome at all are eligible for (re)attempt.
    % A "timeout" IS a recorded outcome (Invariant 8: bounded per-cell work),
    % so timed-out cells are not auto-retried on resume -- rerun explicitly
    % with a fresh computation_id (no resume_from) if a longer timeout is wanted.
    pending_mask = ~existence_mask & (fail_reason == "");
    pend_idx = find(pending_mask);
    if numel(pend_idx) > p.max_cells_this_run
        pend_idx = pend_idx(1:p.max_cells_this_run);   % batch: only the first N this call
    end
    n_this_run = numel(pend_idx);
    fprintf('  [independent] attempting %d/%d pending cell(s) this run\n', ...
            n_this_run, sum(pending_mask(:)));

    if n_this_run > 0
        nG_list = n_GEF_bar + GEF_DIFF(pend_idx)/2;
        nP_list = n_PAK_bar + PAK_DIFF(pend_idx)/2;

        pool = [];
        if p.use_parallel
            try
                pool = gcp;
            catch ME
                warning('build_mpss_branch: could not obtain a parallel pool (%s); falling back to serial.', ...
                        ME.message);
                pool = [];
            end
        end

        if ~isempty(pool)
            % ---- parallel: one parfeval task per pending cell, streamed ----
            futures = parallel.FevalFuture.empty(0, n_this_run);
            for t = 1:n_this_run
                futures(t) = parfeval(pool, @mpss_cell_worker, 1, ...
                    nG_list(t), nP_list(t), L_quarter, par, D, p);
            end
            live = true(1, n_this_run);
            while any(live)
                idxs = find(live);
                [rel_idx, value] = fetchNext(futures(idxs), 2);
                if isempty(rel_idx)
                    for li = idxs
                        % Match datetime('now') to whatever TimeZone (if any)
                        % StartDateTime carries -- MATLAB errors subtracting a
                        % zoned datetime from an unzoned one, and Future
                        % objects' TimeZone can vary by version/platform.
                        % Snapshot StartDateTime ONCE: a future can transition
                        % queued->running (unzoned NaT -> zoned timestamp)
                        % between reads, mixing zoned and unzoned datetimes in
                        % the subtraction below and erroring. Skip futures that
                        % have not actually started.
                        sd = futures(li).StartDateTime;
                        if isnat(sd), continue; end
                        tz = sd.TimeZone;
                        if isempty(tz)
                            now_dt = datetime('now');
                        else
                            now_dt = datetime('now', 'TimeZone', tz);
                        end
                        if seconds(now_dt - sd) > p.cell_timeout
                            cancel(futures(li));
                            live(li) = false;
                            [ii, jj] = ind2sub([res res], pend_idx(li));
                            fail_reason(ii,jj) = "timeout";
                            n_timed_out_this_run = n_timed_out_this_run + 1;
                            n_resolved_since_checkpoint = n_resolved_since_checkpoint + 1;
                            fprintf('  cell(%d,%d): TIMEOUT (> %.0fs)\n', ii, jj, p.cell_timeout);
                            if mod(n_resolved_since_checkpoint, p.checkpoint_every) == 0
                                save_mpss_checkpoint(ckpt_args{:}, PAK_DIFF, GEF_DIFF, ...
                                    existence_mask, SOL, peak_x, amp, mass_residual, seam_slope, ...
                                    sym_residual, flux_values, flux_source_location, fail_reason, ...
                                    par, L_total, L_half, L_quarter, p.M_BG, p.M_P, res, p);
                            end
                        end
                    end
                    continue;
                end
                li = idxs(rel_idx);
                live(li) = false;
                [ii, jj] = ind2sub([res res], pend_idx(li));
                [existence_mask(ii,jj), SOL{ii,jj}, peak_x(ii,jj), amp(ii,jj), ...
                 mass_residual(ii,jj), seam_slope(ii,jj), sym_residual(ii,jj), ...
                 flux_values(ii,jj,:), flux_source_location(ii,jj), fail_reason(ii,jj)] = ...
                    unpack_mpss_result(value);
                log_cell_result(ii, jj, existence_mask(ii,jj), amp(ii,jj), peak_x(ii,jj), ...
                                mass_residual(ii,jj), seam_slope(ii,jj), sym_residual(ii,jj), fail_reason(ii,jj));
                n_resolved_since_checkpoint = n_resolved_since_checkpoint + 1;
                if mod(n_resolved_since_checkpoint, p.checkpoint_every) == 0
                    save_mpss_checkpoint(ckpt_args{:}, PAK_DIFF, GEF_DIFF, ...
                        existence_mask, SOL, peak_x, amp, mass_residual, seam_slope, ...
                        sym_residual, flux_values, flux_source_location, fail_reason, ...
                        par, L_total, L_half, L_quarter, p.M_BG, p.M_P, res, p);
                end
            end
        else
            % ---- serial (use_parallel=false, or no pool available) ----
            % No preemptive per-cell timeout is possible without a worker
            % process to cancel; cell_timeout is not enforced in this path
            % (documented limitation, not a silent gap).
            for t = 1:n_this_run
                [ii, jj] = ind2sub([res res], pend_idx(t));
                value = mpss_cell_worker(nG_list(t), nP_list(t), L_quarter, par, D, p);
                [existence_mask(ii,jj), SOL{ii,jj}, peak_x(ii,jj), amp(ii,jj), ...
                 mass_residual(ii,jj), seam_slope(ii,jj), sym_residual(ii,jj), ...
                 flux_values(ii,jj,:), flux_source_location(ii,jj), fail_reason(ii,jj)] = ...
                    unpack_mpss_result(value);
                log_cell_result(ii, jj, existence_mask(ii,jj), amp(ii,jj), peak_x(ii,jj), ...
                                mass_residual(ii,jj), seam_slope(ii,jj), sym_residual(ii,jj), fail_reason(ii,jj));
                n_resolved_since_checkpoint = n_resolved_since_checkpoint + 1;
                if mod(n_resolved_since_checkpoint, p.checkpoint_every) == 0
                    save_mpss_checkpoint(ckpt_args{:}, PAK_DIFF, GEF_DIFF, ...
                        existence_mask, SOL, peak_x, amp, mass_residual, seam_slope, ...
                        sym_residual, flux_values, flux_source_location, fail_reason, ...
                        par, L_total, L_half, L_quarter, p.M_BG, p.M_P, res, p);
                end
            end
        end
    end
    end   % strategy switch (continuation vs. independent)
    wall_time = toc(t_start);

    % ---- stability: annotation only, never computed for constructed MPSS ----
    is_stable          = false(res, res);   % no MPSS stability check exists
    is_stable_computed = false(res, res);
    stability_source   = 'not_computed';

    % ---- diagnostic masks (mirror Packet 2 style) ----
    mask_quarter_failed = (fail_reason == "hss") | (fail_reason == "pde") | ...
                          (fail_reason == "bvp") | (fail_reason == "not_polarized") | ...
                          (fail_reason == "quarter_sweep_failed");   % Packet 1D continuation
    mask_not_halfpeak   = (fail_reason == "not_halfpeak") | (fail_reason == "orientation");
    mask_kinked_seam    = (fail_reason == "kinked_seam");
    mask_mass_fail      = (fail_reason == "mass_doubling_fail");
    mask_timeout        = (fail_reason == "timeout");
    mask_other_fail     = (fail_reason ~= "") & ~existence_mask & ...
                          ~mask_quarter_failed & ~mask_not_halfpeak & ...
                          ~mask_kinked_seam & ~mask_mass_fail & ~mask_timeout;
    mask_never_attempted = ~existence_mask & (fail_reason == "");

    run_complete = ~any(mask_never_attempted(:));   % every cell at least attempted

    % ---- Level 1/2 checks (validation_protocol.md; mirrors build_branch_library.m) ----
    checks = struct();
    checks.grid_shape_matches = isequal(size(PAK_DIFF), [res res]) && ...
                                isequal(size(GEF_DIFF), [res res]);
    all_masks_partition = existence_mask | mask_quarter_failed | mask_not_halfpeak | ...
                          mask_kinked_seam | mask_mass_fail | mask_timeout | ...
                          mask_other_fail | mask_never_attempted;
    checks.grid_partitioned_by_classes = all(all_masks_partition(:));
    checks.valid_cells_have_SOL = all(~cellfun(@isempty, SOL(existence_mask)));
    checks.invalid_cells_empty_SOL = all(cellfun(@isempty, SOL(~existence_mask)));
    checks.invalid_cells_nan_scalars = all(isnan(amp(~existence_mask))) && ...
                                       all(isnan(mass_residual(~existence_mask))) && ...
                                       all(isnan(peak_x(~existence_mask))) && ...
                                       all(isnan(seam_slope(~existence_mask))) && ...
                                       all(isnan(sym_residual(~existence_mask)));
    fv = reshape(flux_values, res*res, 5);
    checks.flux_finite_only_at_valid = all(isfinite(fv(existence_mask(:),:)), 'all') && ...
                                       all(isnan(fv(~existence_mask(:),:)), 'all');
    checks.flux_labeled_only_at_valid = ...
        all(flux_source_location(existence_mask) == "outer_trough_at_half_domain_interface") && ...
        all(flux_source_location(~existence_mask) == "");
    checks.existence_excludes_amp_below_0p5 = all(amp(existence_mask) >= 0.5);
    ex_peak_frac = peak_x(existence_mask) / L_half;
    checks.valid_mpss_are_interior_peaks = all(ex_peak_frac > 0.4 & ex_peak_frac < 0.6);
    checks.mass_residual_within_tol = all(mass_residual(existence_mask) <= p.mass_tol);
    checks.seam_slope_within_tol = all(seam_slope(existence_mask) <= p.slope_tol);
    checks.symmetry_within_tol = all(sym_residual(existence_mask) <= p.sym_tol);
    n_neg = 0;
    for ii = 1:res
        for jj = 1:res
            if existence_mask(ii,jj) && any(SOL{ii,jj}.y < -1e-9, 'all')
                n_neg = n_neg + 1;
            end
        end
    end
    checks.no_negative_concentrations = (n_neg == 0);
    checks.n_negative_concentration_cells = n_neg;
    checks.is_stable_is_annotation_not_existence = ...
        ~any(is_stable_computed(:)) && strcmp(stability_source, 'not_computed');
    checks.timed_out_cells_are_masked = all(~existence_mask(mask_timeout)) && ...
                                        all(isnan(amp(mask_timeout))) && ...
                                        all(cellfun(@isempty, SOL(mask_timeout)));
    checks.n_valid_MPSS = sum(existence_mask(:));
    checks.n_quarter_failed = sum(mask_quarter_failed(:));
    checks.n_not_halfpeak = sum(mask_not_halfpeak(:));
    checks.n_kinked_seam = sum(mask_kinked_seam(:));
    checks.n_mass_fail = sum(mask_mass_fail(:));
    checks.n_timeout = sum(mask_timeout(:));
    checks.n_other_fail = sum(mask_other_fail(:));
    checks.n_never_attempted = sum(mask_never_attempted(:));
    checks.n_total = res*res;
    checks.info_only_stability_was_computed = false;
    checks.info_only_run_complete = run_complete;

    % Cross-grid check (informational): does this MPSS grid line up with the
    % existing HPSS/HSS branch libraries, if present? Read-only.
    hpss_latest = fullfile(out_dir, 'branch_HPSS_latest.mat');
    hss_latest  = fullfile(out_dir, 'branch_HSS_latest.mat');
    if exist(hpss_latest, 'file')
        hh = load(hpss_latest, 'PAK_DIFF', 'GEF_DIFF');
        checks.info_only_grid_matches_branch_HPSS_latest = ...
            isequal(hh.PAK_DIFF, PAK_DIFF) && isequal(hh.GEF_DIFF, GEF_DIFF);
    end
    if exist(hss_latest, 'file')
        ss = load(hss_latest, 'PAK_DIFF', 'GEF_DIFF');
        checks.info_only_grid_matches_branch_HSS_latest = ...
            isequal(ss.PAK_DIFF, PAK_DIFF) && isequal(ss.GEF_DIFF, GEF_DIFF);
    end

    print_checks('MPSS', checks);

    % ---- provenance ----
    branch_type  = 'MPSS';
    construction = 'reflected_quarter_HPSS';
    M_BG = p.M_BG; M_P = p.M_P;
    if strcmpi(p.mpss_strategy, 'continuation')
        source_path = sprintf('continuation quarter sweep: %s', quarter_source);
    else
        source_path = '(constructed from quarter-domain HPSS solves; no external branch file read)';
    end
    adapter_run_id = computation_id;
    run_status = ternary(run_complete, 'complete', 'partial');

    % ---- final checkpoint (always, so a completed run is also resumable-in-spirit) ----
    checkpoint_path = save_mpss_checkpoint(ckpt_args{:}, PAK_DIFF, GEF_DIFF, ...
        existence_mask, SOL, peak_x, amp, mass_residual, seam_slope, ...
        sym_residual, flux_values, flux_source_location, fail_reason, ...
        par, L_total, L_half, L_quarter, p.M_BG, p.M_P, res, p);

    % ---- save timestamped output (ALWAYS -- honest record of what's resolved) ----
    out_path = fullfile(out_dir, sprintf('branch_MPSS_%s.mat', computation_id));
    save(out_path, 'branch_type', 'construction', 'PAK_DIFF', 'GEF_DIFF', ...
         'existence_mask', 'SOL', 'peak_x', 'amp', 'mass_residual', ...
         'seam_slope', 'sym_residual', 'flux_values', 'flux_source_location', ...
         'is_stable', 'is_stable_computed', 'stability_source', ...
         'mask_quarter_failed', 'mask_not_halfpeak', 'mask_kinked_seam', ...
         'mask_mass_fail', 'mask_timeout', 'mask_other_fail', 'mask_never_attempted', ...
         'fail_reason', 'run_status', ...
         'par', 'L_total', 'L_half', 'L_quarter', 'M_BG', 'M_P', ...
         'computation_id', 'res', 'source_path', 'adapter_run_id');

    % ---- promote to "latest" only if complete, or explicitly accepted partial ----
    promoted = false;
    latest_path = fullfile(out_dir, 'branch_MPSS_latest.mat');
    if run_complete || p.accept_partial
        copyfile(out_path, latest_path);
        promoted = true;
        fprintf('  wrote %s\n  wrote %s (run_status=%s)\n', out_path, latest_path, run_status);
    else
        fprintf(['  wrote %s (run_status=%s; NOT promoted to %s -- %d cell(s) never attempted.\n', ...
                 '  Resume with p.resume_from=''%s'' to continue, or set p.accept_partial=true\n', ...
                 '  to promote this partial result anyway.)\n'], ...
                out_path, run_status, latest_path, sum(mask_never_attempted(:)), checkpoint_path);
    end
    checks.info_only_promoted_to_latest = promoted;

    % ---- report ----
    report = struct();
    report.out_path = out_path;
    report.checkpoint_path = checkpoint_path;
    report.computation_id = computation_id;
    report.res = res;
    report.L_total = L_total; report.L_half = L_half; report.L_quarter = L_quarter;
    report.n_valid_MPSS = sum(existence_mask(:));
    report.n_total = res*res;
    report.n_quarter_failed = sum(mask_quarter_failed(:));
    report.n_not_halfpeak = sum(mask_not_halfpeak(:));
    report.n_kinked_seam = sum(mask_kinked_seam(:));
    report.n_mass_fail = sum(mask_mass_fail(:));
    report.n_timeout = sum(mask_timeout(:));
    report.n_other_fail = sum(mask_other_fail(:));
    report.n_never_attempted = sum(mask_never_attempted(:));
    report.run_complete = run_complete;
    report.run_status = run_status;
    report.promoted_to_latest = promoted;
    report.wall_time = wall_time;
    report.checks = checks;
    % example (first valid cell)
    [vi, vj] = find(existence_mask, 1);
    if ~isempty(vi)
        report.example = struct('i', vi, 'j', vj, ...
            'nG_half', n_GEF_bar + GEF_DIFF(vi,vj)/2, 'nP_half', n_PAK_bar + PAK_DIFF(vi,vj)/2, ...
            'nG_quarter', (n_GEF_bar + GEF_DIFF(vi,vj)/2)/2, 'nP_quarter', (n_PAK_bar + PAK_DIFF(vi,vj)/2)/2, ...
            'amp', amp(vi,vj), 'peak_x', peak_x(vi,vj), 'mass_residual', mass_residual(vi,vj), ...
            'seam_slope', seam_slope(vi,vj), 'sym_residual', sym_residual(vi,vj), ...
            'flux_values_interface', squeeze(flux_values(vi,vj,:))', ...
            'flux_values_outer_x0', SOL{vi,vj}.y(:,1)');
    end
    fprintf('=== build_mpss_branch complete: %d/%d valid MPSS cells (run_status=%s). Wall: %.1fs ===\n', ...
            report.n_valid_MPSS, report.n_total, run_status, wall_time);
end


%% ========================================================================
%  PER-CELL WORKER  (calls the UNCHANGED per-cell math; parfeval-callable)
%  ========================================================================

function result = mpss_cell_worker(nG_half, nP_half, L_quarter, par, D, p)
    % One MPSS cell, end to end: quarter solve -> reflect -> classify. Single
    % struct output so it can be dispatched with parfeval(pool, @mpss_cell_worker,
    % 1, ...). No side effects; safe to run on any worker. The three calls below
    % (quarter_hpss_solve, reflect_stitch, classify_mpss_profile) are the exact
    % same functions used by the original small-grid prototype -- unchanged.
    result = struct('ok', false, 'why', "", 'sol', [], 'peak_x', NaN, 'amp', NaN, ...
                     'mass_residual', NaN, 'seam_slope', NaN, 'sym_residual', NaN, ...
                     'flux_values', nan(1,5));
    nG_q = nG_half / 2;
    nP_q = nP_half / 2;
    [sol_q, ok, why] = quarter_hpss_solve(L_quarter, nG_q, nP_q, par, D, p);
    if ~ok
        result.why = why;
        return;
    end
    [xm, ym, seam_du] = reflect_stitch(sol_q);
    [valid, why2, this_amp, pk_x, mass_res, sym_res] = ...
        classify_mpss_profile(xm, ym, seam_du, nG_half, nP_half, p);
    if ~valid
        result.why = why2;
        return;
    end
    result.ok = true;
    result.sol.x = xm; result.sol.y = ym;
    result.peak_x = pk_x;
    result.amp = this_amp;
    result.mass_residual = mass_res;
    result.seam_slope = seam_du;
    result.sym_residual = sym_res;
    result.flux_values = ym(:, end)';   % outer trough at the half-domain interface
end

function [ex, sol_out, pk_x, this_amp, mass_res, seam_du, sym_res, flux_vec, flux_loc, reason] = ...
        unpack_mpss_result(value)
    % Map a mpss_cell_worker() result struct onto the per-cell schema fields.
    ex = value.ok;
    if value.ok
        sol_out = value.sol;
        pk_x = value.peak_x;
        this_amp = value.amp;
        mass_res = value.mass_residual;
        seam_du = value.seam_slope;
        sym_res = value.sym_residual;
        flux_vec = value.flux_values;
        flux_loc = "outer_trough_at_half_domain_interface";
        reason = "";
    else
        sol_out = [];
        pk_x = NaN; this_amp = NaN; mass_res = NaN; seam_du = NaN; sym_res = NaN;
        flux_vec = nan(1,5);
        flux_loc = "";
        reason = value.why;
    end
end

function log_cell_result(i, j, ex, this_amp, pk_x, mass_res, seam_du, sym_res, reason)
    if ex
        fprintf('  cell(%d,%d): MPSS OK. amp=%.3f peak_x=%.3f mass_res=%.2e seam|du|=%.2e sym=%.2e\n', ...
                i, j, this_amp, pk_x, mass_res, seam_du, sym_res);
    else
        fprintf('  cell(%d,%d): FAILED (%s)\n', i, j, reason);
    end
end


%% ========================================================================
%  CHECKPOINTING  (Invariant 8: bounded, resumable, no growing broadcast)
%  ========================================================================

function ckpt_path = save_mpss_checkpoint(out_dir, computation_id, PAK_DIFF, GEF_DIFF, ...
        existence_mask, SOL, peak_x, amp, mass_residual, seam_slope, sym_residual, ...
        flux_values, flux_source_location, fail_reason, par, L_total, L_half, L_quarter, ...
        M_BG, M_P, res, p)
    % Best-effort checkpoint: overwrites ONE file per computation_id (not
    % versioned per pass), so cost stays bounded regardless of grid size.
    ckpt_path = fullfile(out_dir, sprintf('checkpoint_mpss_%s.mat', computation_id));
    n_resolved = sum(existence_mask(:) | (fail_reason(:) ~= ""));
    try
        save(ckpt_path, 'PAK_DIFF', 'GEF_DIFF', 'existence_mask', 'SOL', 'peak_x', 'amp', ...
             'mass_residual', 'seam_slope', 'sym_residual', 'flux_values', ...
             'flux_source_location', 'fail_reason', 'par', 'L_total', 'L_half', 'L_quarter', ...
             'M_BG', 'M_P', 'res', 'computation_id', 'p', 'n_resolved');
        fprintf('  [checkpoint saved: %s, resolved=%d/%d]\n', ckpt_path, n_resolved, res*res);
    catch ME
        warning('build_mpss_branch: checkpoint save failed: %s', ME.message);
    end
end


%% ========================================================================
%  MPSS-SPECIFIC LOGIC  (UNCHANGED math -- quarter solve, reflection, and
%  acceptance criteria are byte-identical to the validated prototype)
%  ========================================================================

function [sol_q, ok, why] = quarter_hpss_solve(L_q, nG_q, nP_q, par, D, p)
    % Solve a boundary half-peak (HPSS) on the quarter domain [0, L_q] with
    % integral targets (nG_q, nP_q). Seed via HSS -> cos-perturbed PDE relax ->
    % BVP (the sweeper's Phase-0 auto pattern). Peak forced to x=0.
    sol_q = []; ok = false; why = "";
    try
        [u, q, v, w, pv] = get_hss_local(nG_q, nP_q, L_q, par);
    catch
        why = "hss"; return;
    end
    if any(~isfinite([u q v w pv])), why = "hss"; return; end

    x = linspace(0, L_q, 120);
    pde = run_pde_chunked([u; q; v; w; pv], x, L_q, par, D, ...
                          p.pde_T_end, p.pde_T_chunk, p.pde_max_chunks);
    if isempty(pde), why = "pde"; return; end
    if (max(pde.y(1,:)) - min(pde.y(1,:))) < 0.5, why = "not_polarized"; return; end

    [sol, okb] = run_bvp_v7(construct_bvp_guess_v7(pde), L_q, par, D, nG_q, nP_q);
    if ~okb, why = "bvp"; return; end
    if ~validate_halfwave(sol), why = "not_halfpeak"; return; end

    % Force peak-at-x=0 orientation (spec §3). The cos IC boosts u at x=0 so this
    % is normally already true; guard anyway.
    if sol.y(1,1) < sol.y(1,end)
        % spatially reverse the 5 species onto the same increasing mesh
        yr = zeros(size(sol.y));
        yr([1,3,5,7,9],:) = fliplr(sol.y([1,3,5,7,9],:));
        sol.y = yr;   % derivatives dropped; reflect_stitch only uses species rows
    end
    sol_q = sol; ok = true;
end


function [x_mpss, y_mpss, seam_du] = reflect_stitch(sol_q)
    % Pure peak-to-peak reflection of a quarter-domain HPSS (peak at x=0).
    % Returns the stitched MPSS on [0, 2*L_quarter] = [0, L_half] as a 5-row
    % species profile y_mpss = [u;q;v;w;p], and the |du/dx| at the seam (from the
    % quarter's peak-end derivative) for the C1 check.
    x_q  = sol_q.x(:).';                 % 1 x N, [0 .. L_q], peak at column 1
    yq5  = sol_q.y([1,3,5,7,9], :);      % 5 x N
    L_q  = x_q(end);

    % C1 seam slope: du/dx at the quarter peak end (x_q=0). Prefer the BVP's own
    % derivative channel if present (row 2); else finite-difference.
    if size(sol_q.y,1) >= 2 && any(sol_q.y(2,:) ~= 0)
        seam_du = abs(sol_q.y(2,1));
    else
        seam_du = abs((yq5(1,2) - yq5(1,1)) / (x_q(2) - x_q(1)));
    end

    x_right = L_q + x_q;                 % [L_q .. 2 L_q]
    y_right = yq5;
    x_left  = L_q - fliplr(x_q);         % [0 .. L_q]
    y_left  = fliplr(yq5);

    % drop duplicated seam node (x_left(end) == x_right(1) == L_q)
    x_mpss = [x_left(1:end-1), x_right];
    y_mpss = [y_left(:, 1:end-1), y_right];
end


function [ok, why, this_amp, pk_x, mass_res, sym_res] = ...
        classify_mpss_profile(xm, ym, seam_du, nG_half, nP_half, p)
    % Pure validation of a stitched MPSS profile against spec Sec 2 & 4 -- the
    % single source of truth for MPSS acceptance criteria (identical formulas,
    % order, and tolerances as the original inline version; only extracted for
    % reuse/testability, per the 2026-07-09 generalization pass). ym rows =
    % [u;q;v;w;p]. Returns ok=false with a reason string on the first failing
    % check (Invariant 2: mark, don't fabricate); the scalars are still
    % computed so the caller can log/inspect them even on failure.
    ok = false; why = "";
    u = ym(1,:);
    [~, ipk] = max(u);
    pk_x = xm(ipk);
    pk_frac = pk_x / xm(end);
    this_amp = max(u) - min(u);
    % MESH-ACCURATE mass residual (matches run_quarter_hpss_sweeper_for_mpss's
    % quarter_mass_res): a plain trapz on the stitched profile's non-uniform mesh
    % conflates mesh-discretization error (~1e-4) with true mass error (~1e-13).
    % pchip-densify onto a uniform fine mesh before trapz -> ~1e-6. This is an
    % accuracy fix, NOT a tolerance change (Invariant 4); without it, genuinely
    % on-mass MPSS profiles were falsely rejected as mass_doubling_fail.
    xf = linspace(xm(1), xm(end), 4000);
    uf_ = interp1(xm, ym(1,:), xf, 'pchip');   qf_ = interp1(xm, ym(2,:), xf, 'pchip');
    vf_ = interp1(xm, ym(3,:), xf, 'pchip');   wf_ = interp1(xm, ym(4,:), xf, 'pchip');
    pf_ = interp1(xm, ym(5,:), xf, 'pchip');
    gef_mass = trapz(xf, uf_ + qf_ + vf_ + wf_);
    pak_mass = trapz(xf, pf_ + qf_);
    mass_res = max(abs(gef_mass - nG_half), abs(pak_mass - nP_half));
    % symmetry about centre L_half/2
    uf = interp1(xm, u, xm(end) - xm, 'pchip');
    sym_res = max(abs(u - uf));

    neg = any(ym < -1e-9, 'all');
    mono = all(diff(xm) > 0);
    interior = (pk_frac > 0.4 && pk_frac < 0.6);

    if neg
        why = "negative_conc"; return;
    elseif ~mono
        why = "nonmonotone_mesh"; return;
    elseif ~interior
        why = "peak_not_interior"; return;
    elseif this_amp < 0.5
        why = "amp_below_0p5"; return;
    elseif seam_du > p.slope_tol
        why = "kinked_seam"; return;
    elseif mass_res > p.mass_tol
        why = "mass_doubling_fail"; return;
    elseif sym_res > p.sym_tol
        why = "symmetry_fail"; return;
    end
    ok = true;
end


%% ========================================================================
%  VALIDATION / REPORTING HELPERS  (no numerical content -- console reporting
%  only, mirrors build_branch_library.m)
%  ========================================================================

function print_checks(label, c)
    % info_only_* fields are status, not pass/fail -> printed as INFO.
    fprintf('  [%s checks]\n', label);
    fn = fieldnames(c);
    for k = 1:numel(fn)
        val = c.(fn{k});
        is_info = startsWith(fn{k}, 'info_only_');
        if islogical(val) && isscalar(val)
            if is_info, tag = ternary(val, 'INFO (true)', 'INFO (false)');
            else,       tag = ternary(val, 'PASS', 'FAIL');
            end
            fprintf('    %-46s %s\n', fn{k}, tag);
        else
            fprintf('    %-46s %s\n', fn{k}, mat2str(val));
        end
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end


%% ========================================================================
%  BVP CORE: standalone files in src/solvers/ (Packet 1B, 2026-07-09) --
%  run_bvp_v7, ode_sys_v7, bc_aug_v7, construct_bvp_guess_v7, run_pde_chunked,
%  pde_eq_v7, pde_bc_v7, pde_ic_local, get_hss_local, hss_eqs, validate_halfwave.
%  Shared single source of truth with run_hybrid_sweeper_pdefb.m; not duplicated here.
%  ========================================================================
