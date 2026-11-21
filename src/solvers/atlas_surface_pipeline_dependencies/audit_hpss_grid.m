function [class_grid, reason_grid, diagnostics_grid, invalid_mask] = ...
    audit_hpss_grid(SOL_grid, failure_mask, class_grid, reason_grid, diagnostics_grid)
%AUDIT_HPSS_GRID Revalidate every half-peak profile that may serve as a starting profile before continuation.
%
% The caller must discard invalid_mask cells: clear SOL_grid, set the
% failure mask, and clear derived trough/stability values. Existing diagnostic
% labels are retained for already-quarantined cells.

    grid_size = size(SOL_grid);
    if nargin < 3 || isempty(class_grid)
        class_grid = repmat("NOT_FOUND", grid_size);
    end
    if nargin < 4 || isempty(reason_grid)
        reason_grid = repmat("no_accepted_solution", grid_size);
    end
    if nargin < 5 || isempty(diagnostics_grid)
        diagnostics_grid = cell(grid_size);
    end

    invalid_mask = false(grid_size);
    solved_idx = find(~failure_mask);
    for k = 1:numel(solved_idx)
        idx = solved_idx(k);
        [ok, info] = validate_halfwave(SOL_grid{idx});
        class_grid(idx) = info.class;
        reason_grid(idx) = info.reason;
        diagnostics_grid{idx} = info;
        invalid_mask(idx) = ~ok;
    end
end
