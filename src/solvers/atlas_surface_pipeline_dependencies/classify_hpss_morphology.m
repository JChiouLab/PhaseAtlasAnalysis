function info = classify_hpss_morphology(sol, opts)
%CLASSIFY_HPSS_MORPHOLOGY Classify a half-domain BVP profile (before it is used as a starting profile).
%
% A valid half-peak (HPSS) state has one peak at a boundary and its minimum at the
% opposite boundary.  The peak may lie slightly inside the boundary,
% but an interior minimum, a second peak, or a dominant
% interior peak is not a half-peak shape.
%
% Classes: HPSS (half-peak), HSS (flat, no-peak), MPSS_LIKE (single interior peak),
% MULTIPEAK (several peaks), UNKNOWN (polarized but not a boundary peak, or ill-formed).

    if nargin < 2, opts = struct(); end
    opts = set_default(opts, 'amplitude_tol', 0.5);
    opts = set_default(opts, 'edge_fraction', 0.10);
    opts = set_default(opts, 'trough_rel_tol', 0.02);
    opts = set_default(opts, 'prominence_rel', 0.05);
    opts = set_default(opts, 'absolute_tol', 1e-3);
    opts = set_default(opts, 'autocorr_min', 0.5);
    opts = set_default(opts, 'negative_tol', -1e-3);

    info = empty_info(opts);
    if isempty(sol) || ~isstruct(sol) || ~isfield(sol, 'x') || ...
            ~isfield(sol, 'y') || isempty(sol.x) || isempty(sol.y) || ...
            size(sol.y, 1) < 1 || numel(sol.x) ~= size(sol.y, 2)
        info.reason = "malformed_solution";
        return;
    end
    if any(~isfinite(sol.x(:))) || any(~isfinite(sol.y(:)))
        info.reason = "nonfinite_solution";
        return;
    end
    concentration_rows = [1 3 5 7 9];
    concentration_rows = concentration_rows(concentration_rows <= size(sol.y, 1));
    if any(sol.y(concentration_rows, :) < opts.negative_tol, 'all')
        info.reason = "negative_concentration";
        return;
    end

    x = sol.x(:).';
    u = sol.y(1, :);
    span = x(end) - x(1);
    if span <= 0
        info.reason = "nonpositive_domain_length";
        return;
    end

    % Inspect either mirror orientation in one canonical orientation:
    % boundary peak at the left and minimum at the far boundary on the right.
    if u(1) < u(end)
        u = fliplr(u);
        info.orientation_flipped = true;
    end

    info.amp = max(u) - min(u);
    [~, ipk] = max(u);
    info.peak_frac = (x(ipk) - x(1)) / span;
    info.peak_x_canonical = x(ipk);

    if info.amp < opts.amplitude_tol
        info.class = "HSS";
        info.reason = "amplitude_below_threshold";
        return;
    end

    prominence = max(opts.prominence_rel * info.amp, opts.absolute_tol);
    [~, loc_max] = findpeaks(u, 'MinPeakProminence', prominence);
    [~, loc_min] = findpeaks(-u, 'MinPeakProminence', prominence);
    info.n_interior_max = numel(loc_max);
    info.n_interior_min = numel(loc_min);

    trough_tol = max(opts.trough_rel_tol * info.amp, opts.absolute_tol);
    info.far_trough_gap = u(end) - min(u);
    info.far_boundary_trough = info.far_trough_gap <= trough_tol;
    info.boundary_associated_peak = info.peak_frac < opts.edge_fraction;

    near_edge_crown = false;
    if isscalar(loc_max)
        crown_frac = (x(loc_max) - x(1)) / span;
        near_edge_crown = crown_frac < opts.edge_fraction;
    end

    uc = u - mean(u);
    if std(uc) < 1e-9
        info.autocorr_lag1 = 1;
    else
        ac = xcorr(uc, 'coeff');
        half = ac(numel(uc):end);
        if numel(half) >= 2, info.autocorr_lag1 = half(2); end
    end
    smooth_enough = isfinite(info.autocorr_lag1) && ...
                    info.autocorr_lag1 > opts.autocorr_min;

    if info.boundary_associated_peak && info.far_boundary_trough && ...
            info.n_interior_min == 0 && ...
            (info.n_interior_max == 0 || near_edge_crown) && smooth_enough
        info.class = "HPSS";
        info.accepted = true;
        if info.n_interior_max == 1
            info.reason = "half_volcano_ok";
        else
            info.reason = "boundary_half_peak_ok";
        end
    elseif info.boundary_associated_peak && ...
            (~info.far_boundary_trough || info.n_interior_min > 0 || ...
             info.n_interior_max >= 2)
        info.class = "MULTIPEAK";
        if ~info.far_boundary_trough
            info.reason = "trough_not_at_far_boundary";
        elseif info.n_interior_min > 0
            info.reason = "significant_interior_trough";
        else
            info.reason = "multiple_significant_crowns";
        end
    elseif info.n_interior_max == 1 && ~info.boundary_associated_peak
        info.class = "MPSS_LIKE";
        info.reason = "dominant_peak_is_interior";
    elseif ~smooth_enough
        info.reason = "insufficient_profile_smoothness";
    else
        info.reason = "polarized_unclassified_morphology";
    end
end


function info = empty_info(opts)
    info = struct( ...
        'accepted', false, ...
        'class', "UNKNOWN", ...
        'reason', "not_evaluated", ...
        'amp', NaN, ...
        'peak_frac', NaN, ...
        'peak_x_canonical', NaN, ...
        'orientation_flipped', false, ...
        'boundary_associated_peak', false, ...
        'far_boundary_trough', false, ...
        'far_trough_gap', NaN, ...
        'n_interior_max', NaN, ...
        'n_interior_min', NaN, ...
        'autocorr_lag1', NaN, ...
        'thresholds', opts);
end


function s = set_default(s, name, value)
    if ~isfield(s, name), s.(name) = value; end
end
