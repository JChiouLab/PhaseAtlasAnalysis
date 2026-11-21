function y = reconstruct_stitched_sheet(A, branches, sheet_name, i, j, x_full)
    % RECONSTRUCT_STITCHED_SHEET  Build the full-domain [0, L_total] spatial
    % profile of sheet <L>+<R> at grid cell (i,j), on the mesh x_full.
    %
    % Note on the paper: this stitching belongs to the uncoupled approximation (two isolated
    % half-cell states glued at the middle) and is not used for any figure of the paper; the
    % paper solves the coupled quasi-steady states directly.
    %
    % The single definition of sheet reconstruction,
    % shared by run_atlas_pde_trace.m (which uses it to build a PDE initial
    % condition) and classify_pde_to_atlas_sheets.m (which uses it to compare a
    % PDE frame against candidate sheets). Because both call this function,
    % a PDE trace whose initial condition is sheet S at cell C reproduces, at t=0, exactly the
    % classifier's reconstruction of S at C, so the frame-0 residual is ~0.
    %
    % Stitching convention:
    %   * left half occupies [0, L_half]: library L's half-domain profile at cell
    %     (i,j), peak at the outer edge x=0, trough at the interface x=L_half.
    %   * right half occupies [L_half, L_total]: library R's half-domain profile at
    %     the mirror cell (res+1-i, res+1-j), spatially flipped via x_eval =
    %     L_total - x, so its trough sits at the interface and its peak at the
    %     outer edge x=L_total. (Exact discrete mirror, see
    %     build_atlas_sheets.m; no interpolation across grid cells.)
    %
    % No cross-cell interpolation: only the one stored per-cell
    % profile is evaluated at mesh points (via interp1 for half-peak/full-peak states,
    % constant for no-peak states). If the sheet does not
    % exist at (i,j) (either half's existence_mask false), returns [] (missing
    % cells stay unavailable).
    %
    % Args:
    %   A          - loaded atlas_sheets_*.mat struct (needs .sheets, .L_total, .res)
    %   branches   - struct with fields HSS, HPSS, MPSS, each a loaded
    %                branch_<TYPE>_latest.mat struct (needs .SOL, .existence_mask)
    %   sheet_name - one of A.sheet_names, e.g. 'HPSS_MPSS'
    %   i, j       - grid indices
    %   x_full     - 1xN mesh over [0, L_total] to evaluate on
    %
    % Returns y : 5 x numel(x_full) species profile [u;q;v;w;p], or [] if absent.

    y = [];
    sh = A.sheets.(sheet_name);
    if ~sh.sheet_exists(i, j), return; end

    Lt = sh.left_type; Rt = sh.right_type;
    res = A.res;
    mi = res + 1 - i; mj = res + 1 - j;      % exact discrete mirror cell
    L_total = A.L_total; L_half = L_total / 2;

    x_full = x_full(:).';
    idx_L = x_full <= L_half;
    idx_R = ~idx_L;

    yL = extract_half_profile(branches.(Lt), Lt, i,  j,  x_full(idx_L));
    yR = extract_half_profile(branches.(Rt), Rt, mi, mj, L_total - x_full(idx_R));
    if isempty(yL) || isempty(yR)
        y = [];     % should not happen when sheet_exists is true, but stay safe
        return;
    end

    y = zeros(5, numel(x_full));
    y(:, idx_L) = yL;
    y(:, idx_R) = yR;
    y = max(0, y);      % positivity guard
end


function yq = extract_half_profile(branch, type, i, j, xq)
    % Species [u;q;v;w;p] of one half-domain branch profile at cell (i,j),
    % evaluated at query points xq in [0, L_half]. Returns [] if that cell's
    % stored solution is empty.
    s = branch.SOL{i, j};
    if isempty(s), yq = []; return; end
    xq = xq(:).';

    switch type
        case 'HSS'
            % Uniform: SOL{i,j} is a 1x5 (or 5x1) constant [u,q,v,w,p].
            v = s(:);
            yq = repmat(v, 1, numel(xq));

        case 'HPSS'
            % 14-var bvp struct: species at rows 1,3,5,7,9 on mesh s.x=[0,L_half].
            sx = s.x(:).';
            spec = s.y([1 3 5 7 9], :);
            xc = min(max(xq, sx(1)), sx(end));     % clamp to avoid endpoint overshoot NaN
            yq = interp1(sx, spec.', xc, 'pchip').';

        case 'MPSS'
            % Constructed struct: 5-row species on mesh s.x=[0,L_half].
            sx = s.x(:).';
            xc = min(max(xq, sx(1)), sx(end));
            yq = interp1(sx, s.y.', xc, 'pchip').';

        otherwise
            error('reconstruct_stitched_sheet:bad_type', 'Unknown branch type %s', type);
    end
end
