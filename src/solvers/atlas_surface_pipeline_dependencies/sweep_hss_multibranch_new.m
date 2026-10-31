function sweep_hss_multibranch_new(p_in)
    % SWEEP_HSS_MULTIBRANCH_NEW
    % Multi-branch HSS sweep with continuity-tracked branch labels.
    %
    % This version keeps the generic multistart fsolve architecture, but it
    % separates root discovery from branch labelling:
    %   1. At each (Delta n_PAK, Delta n_GEF) cell, discover all real-positive
    %      HSS roots using multistart fsolve plus neighbour-seeded guesses.
    %   2. Store roots in roots_cell{i,j} without assigning branch indices.
    %   3. Assign global branch IDs by nearest-neighbour continuity in the
    %      full HSS state [u q v w p], not by local u-rank.
    %   4. Save the usual U/Q/V/W/P_hss_branches arrays, now indexed by
    %      tracked branch labels.
    %
    % Output schema (saved to data/outputs/hss_sweep_multibranch_results.mat):
    %   PAK_DIFF, GEF_DIFF
    %   U_hss_branches, Q_hss_branches, V_hss_branches, W_hss_branches,
    %   P_hss_branches                  - res x res x max_branches arrays.
    %                                     Branch index k is a tracked global
    %                                     branch label, ranked by coverage.
    %   n_branches_grid                 - res x res count of STORED branches
    %                                     per cell (non-NaN slices in
    %                                     U_hss_branches after Stage 3
    %                                     pruning to the top max_branches by
    %                                     coverage). Use n_raw_roots_grid for
    %                                     the raw discovered-root count,
    %                                     which can exceed this if a cell's
    %                                     root belongs to a low-coverage
    %                                     branch that got pruned.
    %   n_raw_roots_grid                - res x res raw discovered root
    %                                     count, pre-pruning.
    %   branch_id_grid                  - res x res x max_branches tracked IDs.
    %   branch_summary                  - table-like matrix:
    %                                     [stored_k, branch_id, coverage,
    %                                      median_u, median_p]
    %   roots_cell, branch_id_cell      - diagnostics for full root discovery.
    %   residual_cell                   - accepted residuals per cell.
    %   par, L_total, M_BG, M_P, computation_id
    %
    % Override-able fields on p:
    %   p.guess_u_grid        - vector of initial-guess u values
    %   p.guess_p_grid        - vector of initial-guess p values
    %   p.max_branches        - max tracked branches stored (default 6)
    %   p.dedup_tol           - duplicate tolerance in log-state distance
    %                           during root discovery (default 1e-6)
    %   p.res_tol             - max |F| residual for accepting fsolve answer
    %                           (default 1e-8)
    %   p.track_tol           - log-state distance cutoff for branch matching
    %                           (default 0.35)
    %   p.positive_floor      - lower positivity bound for [u p] (default 1e-10)
    %   p.use_neighbor_guesses - append roots from left/above cells as guesses
    %                            during discovery (default true)

    computation_id = datestr(now, 'yyyymmdd_HHMMSS');

    % --- 1. Parameters ---
    if nargin < 1
        % V7 Gold Standard
        p.a = 0.59; p.b = 0.13; p.c = 21.0;
        p.d = 1.3;  p.e = 33.0; p.f = 1.9;
        p.Dm = 0.001; p.Dc = 1.0;
        p.L = 4.0;
        p.M_BG = 16.0;
        p.M_P  = 8.0;
        p.res = 25;
    else
        p = p_in;
    end

    if ~isfield(p, 'guess_u_grid'),         p.guess_u_grid = [0.005 0.05 0.2 0.5 1.0 2.0 5.0 10.0]; end
    if ~isfield(p, 'guess_p_grid'),         p.guess_p_grid = [0.005 0.05 0.2 0.5 1.0 2.0 5.0 10.0]; end
    if ~isfield(p, 'max_branches'),         p.max_branches = 6; end
    if ~isfield(p, 'dedup_tol'),            p.dedup_tol = 1e-6; end
    if ~isfield(p, 'res_tol'),              p.res_tol = 1e-8; end
    if ~isfield(p, 'track_tol'),            p.track_tol = 0.35; end
    if ~isfield(p, 'positive_floor'),       p.positive_floor = 1e-10; end
    if ~isfield(p, 'use_neighbor_guesses'), p.use_neighbor_guesses = true; end

    par = p;
    L_total = p.L;
    res = p.res;
    L_half = L_total / 2;

    n_GEF_bar = p.M_BG / 2;
    n_PAK_bar = p.M_P / 2;

    deln_pak = linspace(-0.95, 0.95, res) * p.M_P;
    deln_G   = linspace(-0.95, 0.95, res) * p.M_BG;
    [PAK_DIFF, GEF_DIFF] = meshgrid(deln_pak, deln_G);

    M_max = p.max_branches;

    % Initial-guess grid: Cartesian product of u and p guess vectors.
    [GU, GP] = meshgrid(p.guess_u_grid, p.guess_p_grid);
    base_guesses = [GU(:), GP(:)];
    n_base_guesses = size(base_guesses, 1);

    roots_cell = cell(res, res);       % each entry: n_roots x 5, [u q v w p]
    residual_cell = cell(res, res);    % accepted max-norm residuals
    n_branches_grid = zeros(res, res);
    n_raw_roots_grid = zeros(res, res);

    fprintf('Tracked multibranch HSS sweep: %d x %d, %d base guesses/cell, max_branches=%d\n', ...
            res, res, n_base_guesses, M_max);
    fprintf('  dedup_tol=%g, res_tol=%g, track_tol=%g\n', ...
            p.dedup_tol, p.res_tol, p.track_tol);

    opts = optimset('Display', 'off', ...
                    'TolFun', 1e-10, ...
                    'TolX', 1e-10, ...
                    'MaxIter', 1000, ...
                    'MaxFunEvals', 5000);

    % ------------------------------------------------------------------
    % Stage 1. Discover all roots at every cell without branch labelling.
    % ------------------------------------------------------------------
    for i = 1:res
        for j = 1:res
            TG = n_GEF_bar + GEF_DIFF(i, j) / 2;
            TP = n_PAK_bar + PAK_DIFF(i, j) / 2;

            local_guesses = base_guesses;

            % Add already-discovered neighbouring roots as guesses. This
            % improves root recall while keeping the generic fsolve structure.
            if p.use_neighbor_guesses
                if j > 1 && ~isempty(roots_cell{i, j-1})
                    local_guesses = [local_guesses; roots_cell{i, j-1}(:, [1, 5])]; %#ok<AGROW>
                end
                if i > 1 && ~isempty(roots_cell{i-1, j})
                    local_guesses = [local_guesses; roots_cell{i-1, j}(:, [1, 5])]; %#ok<AGROW>
                end
                if i > 1 && j > 1 && ~isempty(roots_cell{i-1, j-1})
                    local_guesses = [local_guesses; roots_cell{i-1, j-1}(:, [1, 5])]; %#ok<AGROW>
                end
                if i > 1 && j < res && ~isempty(roots_cell{i-1, j+1})
                    local_guesses = [local_guesses; roots_cell{i-1, j+1}(:, [1, 5])]; %#ok<AGROW>
                end
            end

            roots_state = zeros(0, 5);
            roots_resid = zeros(0, 1);

            for g = 1:size(local_guesses, 1)
                y0 = local_guesses(g, :);

                [vals, fval, exitflag] = fsolve(@(y) hss_eqs(y, par, L_half, TG, TP), ...
                                                y0, opts);

                if exitflag <= 0,                    continue; end
                if any(~isfinite(vals)),             continue; end
                if any(vals <= p.positive_floor),    continue; end

                resid = max(abs(fval));
                if resid > p.res_tol,                continue; end

                state = make_hss_state(vals, par);   % [u q v w p]
                if any(~isfinite(state)) || any(state <= 0)
                    continue;
                end

                % Deduplicate in full log-state distance, not only [u p].
                is_dup = false;
                for k = 1:size(roots_state, 1)
                    dist = normalized_state_distance(state, roots_state(k, :));
                    if dist < p.dedup_tol
                        is_dup = true;
                        % Keep the numerically better duplicate.
                        if resid < roots_resid(k)
                            roots_state(k, :) = state;
                            roots_resid(k) = resid;
                        end
                        break;
                    end
                end

                if ~is_dup
                    roots_state(end+1, :) = state; %#ok<AGROW>
                    roots_resid(end+1, 1) = resid; %#ok<AGROW>
                end
            end

            % Local sorting is only for deterministic diagnostics. Final
            % branch labels are assigned later by continuity tracking.
            if ~isempty(roots_state)
                [~, sort_idx] = sort(roots_state(:, 1));
                roots_state = roots_state(sort_idx, :);
                roots_resid = roots_resid(sort_idx);
            end

            roots_cell{i, j} = roots_state;
            residual_cell{i, j} = roots_resid;
            n_raw_roots_grid(i, j) = size(roots_state, 1);
            n_branches_grid(i, j) = size(roots_state, 1);
        end

        if mod(i, 5) == 0 || i == res
            fprintf('  Discovery row %d/%d done. Cells with: 1=%d, 2=%d, 3+=%d, max=%d\n', ...
                    i, res, ...
                    sum(n_branches_grid(:) == 1), ...
                    sum(n_branches_grid(:) == 2), ...
                    sum(n_branches_grid(:) >= 3), ...
                    max(n_branches_grid(:)));
        end
    end

    % ------------------------------------------------------------------
    % Stage 2. Assign global branch IDs by continuity.
    % ------------------------------------------------------------------
    % If every cell has exactly one accepted positive root, there is no
    % branch-assignment ambiguity: the connected parameter grid carries one
    % unique HSS representative everywhere.  Applying a finite log-distance
    % cutoff in that case can split a steep-but-unique surface into artificial
    % IDs; subsequent max_branches pruning then creates false nonexistence
    % seams.  Bypass tracking only in this unambiguous all-unique case.
    if all(n_raw_roots_grid(:) == 1)
        branch_id_cell = cell(res, res);
        for i = 1:res
            for j = 1:res
                branch_id_cell{i,j} = 1;
            end
        end
        fprintf(['  Every cell has exactly one raw HSS root: assigned one ', ...
                 'global branch ID without distance-threshold splitting.\n']);
    else
        branch_id_cell = track_hss_branches(roots_cell, p.track_tol);
    end

    all_ids = collect_branch_ids(branch_id_cell);
    n_global_branches = numel(all_ids);

    fprintf('  Branch tracking complete. Global tracked branches found: %d\n', ...
            n_global_branches);

    % Rank tracked branches for storage. Prefer branches covering more cells;
    % use median u as a deterministic tie-breaker.
    branch_stats = compute_branch_stats(roots_cell, branch_id_cell, all_ids);
    % branch_stats columns: [branch_id, coverage, median_u, median_p]

    if isempty(branch_stats)
        kept_ids = [];
        branch_summary = zeros(0, 5);
    else
        [~, order] = sortrows([-branch_stats(:, 2), branch_stats(:, 3)]);
        branch_stats = branch_stats(order, :);

        if size(branch_stats, 1) > M_max
            warning(['Found %d tracked branches, but max_branches=%d. ', ...
                     'Stored the %d highest-coverage branches. Increase ', ...
                     'p.max_branches if lower-coverage branches matter.'], ...
                    size(branch_stats, 1), M_max, M_max);
        end

        n_keep = min(M_max, size(branch_stats, 1));
        kept_stats = branch_stats(1:n_keep, :);
        kept_ids = kept_stats(:, 1);

        branch_summary = zeros(n_keep, 5);
        for k = 1:n_keep
            branch_summary(k, :) = [k, kept_stats(k, :)];
        end
        % branch_summary columns: [stored_k, branch_id, coverage, median_u, median_p]
    end

    % ------------------------------------------------------------------
    % Stage 3. Fill output arrays using tracked branch labels.
    % ------------------------------------------------------------------
    U_hss_branches = nan(res, res, M_max);
    Q_hss_branches = nan(res, res, M_max);
    V_hss_branches = nan(res, res, M_max);
    W_hss_branches = nan(res, res, M_max);
    P_hss_branches = nan(res, res, M_max);
    branch_id_grid = nan(res, res, M_max);

    for i = 1:res
        for j = 1:res
            roots = roots_cell{i, j};
            ids = branch_id_cell{i, j};

            for r = 1:size(roots, 1)
                stored_k = find(kept_ids == ids(r), 1);
                if isempty(stored_k)
                    continue;
                end

                U_hss_branches(i, j, stored_k) = roots(r, 1);
                Q_hss_branches(i, j, stored_k) = roots(r, 2);
                V_hss_branches(i, j, stored_k) = roots(r, 3);
                W_hss_branches(i, j, stored_k) = roots(r, 4);
                P_hss_branches(i, j, stored_k) = roots(r, 5);
                branch_id_grid(i, j, stored_k) = ids(r);
            end
        end
    end

    % n_branches_grid (as saved) must reflect what's ACTUALLY stored, not
    % the raw pre-pruning discovery count set during Stage 1 -- otherwise a
    % cell whose only root belongs to a pruned (low-coverage) branch reports
    % a nonzero branch count while every slice at that cell is NaN. The raw
    % count remains available separately in n_raw_roots_grid.
    n_branches_grid = sum(~isnan(U_hss_branches), 3);

    % ------------------------------------------------------------------
    % Save.
    % ------------------------------------------------------------------
    if isfield(p, 'out_dir') && ~isempty(p.out_dir), out_dir = p.out_dir; else, out_dir = 'data/outputs'; end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    filename = sprintf('HSS_Multibranch_%s.mat', computation_id);
    out_path = fullfile(out_dir, filename);

    M_BG = p.M_BG;
    M_P = p.M_P;

    save(out_path, 'PAK_DIFF', 'GEF_DIFF', ...
         'U_hss_branches', 'Q_hss_branches', 'V_hss_branches', ...
         'W_hss_branches', 'P_hss_branches', ...
         'n_branches_grid', 'n_raw_roots_grid', ...
         'branch_id_grid', 'branch_summary', ...
         'roots_cell', 'branch_id_cell', 'residual_cell', ...
         'L_total', 'M_BG', 'M_P', 'par', 'computation_id');

    copyfile(out_path, fullfile(out_dir, 'hss_sweep_multibranch_results.mat'));

    if exist('log_metadata', 'file') == 2
        log_metadata(computation_id, filename, 'Tracked Multibranch HSS Sweep', ...
                     par, L_total, M_BG, M_P, M_BG/L_total, M_P/L_total, ...
                     sprintf(['max_branches=%d, base_guesses=%d, ', ...
                              'global_branches=%d, max_found=%d, track_tol=%g'], ...
                             M_max, n_base_guesses, n_global_branches, ...
                             max(n_branches_grid(:)), p.track_tol));
    else
        warning('log_metadata.m not found on path; skipped metadata logging.');
    end

    fprintf('Tracked multibranch HSS sweep complete.\n');
    fprintf('  Raw max roots/cell: %d\n', max(n_branches_grid(:)));
    fprintf('  Global tracked branches: %d\n', n_global_branches);
    fprintf('  Stored branches: %d/%d\n', numel(kept_ids), M_max);
    fprintf('  Cells: 1 branch=%d, 2=%d, 3=%d, 4+=%d.\n', ...
            sum(n_branches_grid(:) == 1), ...
            sum(n_branches_grid(:) == 2), ...
            sum(n_branches_grid(:) == 3), ...
            sum(n_branches_grid(:) >= 4));
    fprintf('  Output: %s\n', out_path);
end


%% ------------------------------------------------------------------------

function F = hss_eqs(y, par, L, TG, TP)
    u  = y(1);
    pv = y(2);

    q  = (par.c / par.d) * u * pv;
    w  = (par.d / par.e) * q;
    v  = (par.b * u + par.d * q) / (par.a * u^2 + par.f);

    F = [L * (u + q + v + w) - TG;
         L * (pv + q) - TP];
end


function state = make_hss_state(y, par)
    u  = y(1);
    pv = y(2);

    q  = (par.c / par.d) * u * pv;
    w  = (par.d / par.e) * q;
    v  = (par.b * u + par.d * q) / (par.a * u^2 + par.f);

    state = [u, q, v, w, pv];
end


function d = normalized_state_distance(x, y)
    % Log-state distance for positive concentration vectors.
    % This avoids dominance by the largest species and is more robust when
    % concentrations span orders of magnitude.
    floor_val = 1e-12;
    lx = log(max(x, floor_val));
    ly = log(max(y, floor_val));
    d = norm(lx - ly) / sqrt(numel(x));
end


function branch_id_cell = track_hss_branches(roots_cell, track_tol)
    [res_i, res_j] = size(roots_cell);
    branch_id_cell = cell(res_i, res_j);
    next_branch_id = 1;

    for i = 1:res_i
        % Serpentine scan reduces long jumps at row transitions.
        if mod(i, 2) == 1
            j_list = 1:res_j;
        else
            j_list = res_j:-1:1;
        end

        for jj = 1:numel(j_list)
            j = j_list(jj);
            roots = roots_cell{i, j};
            n_roots = size(roots, 1);
            ids = zeros(n_roots, 1);

            if n_roots == 0
                branch_id_cell{i, j} = ids;
                continue;
            end

            [candidate_states, candidate_ids] = collect_neighbor_candidates( ...
                i, j, roots_cell, branch_id_cell);

            if isempty(candidate_states)
                for r = 1:n_roots
                    ids(r) = next_branch_id;
                    next_branch_id = next_branch_id + 1;
                end
                branch_id_cell{i, j} = ids;
                continue;
            end

            n_candidates = size(candidate_states, 1);
            C = inf(n_roots, n_candidates);
            for r = 1:n_roots
                for c = 1:n_candidates
                    C(r, c) = normalized_state_distance(roots(r, :), candidate_states(c, :));
                end
            end

            % Greedy minimum-distance matching. Number of HSS roots per cell is
            % small, so this is sufficient and avoids a toolbox dependency.
            assigned_roots = false(n_roots, 1);
            blocked_candidates = false(n_candidates, 1);

            while true
                C_work = C;
                C_work(assigned_roots, :) = inf;
                C_work(:, blocked_candidates) = inf;

                [best_val, linear_idx] = min(C_work(:));
                if ~isfinite(best_val) || best_val > track_tol
                    break;
                end

                [r_best, c_best] = ind2sub(size(C_work), linear_idx);
                proposed_id = candidate_ids(c_best);

                if ~any(ids == proposed_id)
                    ids(r_best) = proposed_id;
                    assigned_roots(r_best) = true;
                    % Block every candidate carrying this branch ID to enforce
                    % one root per branch per cell.
                    blocked_candidates = blocked_candidates | (candidate_ids == proposed_id);
                else
                    blocked_candidates(c_best) = true;
                end
            end

            for r = 1:n_roots
                if ids(r) == 0
                    ids(r) = next_branch_id;
                    next_branch_id = next_branch_id + 1;
                end
            end

            branch_id_cell{i, j} = ids;
        end
    end
end


function [candidate_states, candidate_ids] = collect_neighbor_candidates(i, j, roots_cell, branch_id_cell)
    [res_i, res_j] = size(roots_cell);

    candidate_offsets = [
        -1,  0;
        -1, -1;
        -1,  1;
         0, -1;
         0,  1
    ];

    candidate_states = zeros(0, 5);
    candidate_ids = zeros(0, 1);

    for k = 1:size(candidate_offsets, 1)
        ii = i + candidate_offsets(k, 1);
        jj = j + candidate_offsets(k, 2);

        if ii < 1 || ii > res_i || jj < 1 || jj > res_j
            continue;
        end

        % Already visited means branch_id_cell{ii,jj} has been assigned.
        if isempty(branch_id_cell{ii, jj})
            continue;
        end

        nb_roots = roots_cell{ii, jj};
        nb_ids = branch_id_cell{ii, jj};

        for r = 1:size(nb_roots, 1)
            if nb_ids(r) > 0
                candidate_states(end+1, :) = nb_roots(r, :); %#ok<AGROW>
                candidate_ids(end+1, 1) = nb_ids(r); %#ok<AGROW>
            end
        end
    end
end


function all_ids = collect_branch_ids(branch_id_cell)
    all_ids = [];
    for i = 1:size(branch_id_cell, 1)
        for j = 1:size(branch_id_cell, 2)
            ids = branch_id_cell{i, j};
            if ~isempty(ids)
                all_ids = [all_ids; ids(:)]; %#ok<AGROW>
            end
        end
    end
    all_ids = unique(all_ids);
    all_ids = all_ids(all_ids > 0);
end


function branch_stats = compute_branch_stats(roots_cell, branch_id_cell, all_ids)
    branch_stats = zeros(numel(all_ids), 4);

    for b = 1:numel(all_ids)
        bid = all_ids(b);
        u_values = [];
        p_values = [];

        for i = 1:size(roots_cell, 1)
            for j = 1:size(roots_cell, 2)
                roots = roots_cell{i, j};
                ids = branch_id_cell{i, j};

                if isempty(roots) || isempty(ids)
                    continue;
                end

                idx = find(ids == bid, 1);
                if ~isempty(idx)
                    u_values(end+1, 1) = roots(idx, 1); %#ok<AGROW>
                    p_values(end+1, 1) = roots(idx, 5); %#ok<AGROW>
                end
            end
        end

        branch_stats(b, :) = [bid, numel(u_values), median(u_values), median(p_values)];
    end
end
