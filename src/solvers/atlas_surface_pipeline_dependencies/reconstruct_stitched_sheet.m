function y = reconstruct_stitched_sheet(A, branches, sheet_name, i, j, x_full)
    % RECONSTRUCT_STITCHED_SHEET  Build the full-domain [0, L_total] spatial
    % profile of atlas sheet <L>+<R> at grid cell (i,j), on the mesh x_full.
    %
    % This is the SINGLE SOURCE OF TRUTH for atlas-sheet spatial reconstruction,
    % shared by run_atlas_pde_trace.m (which uses it to build a PDE initial
    % condition) and classify_pde_to_atlas_sheets.m (which uses it to compare a
    % PDE frame against candidate sheets). Because both call this exact function,
    % a PDE trace whose IC is sheet S at cell C reproduces, at t=0, exactly the
    % classifier's reconstruction of S at C -- so the frame-0 residual is ~0
    % (the Invariant-7 validation anchor).
    %
    % Stitching convention (identical to get_pde_trace.m's IC construction):
    %   * left half occupies [0, L_half]: branch L's half-domain profile at cell
    %     (i,j), peak at the outer edge x=0, trough at the interface x=L_half.
    %   * right half occupies [L_half, L_total]: branch R's half-domain profile at
    %     the MIRROR cell (res+1-i, res+1-j), spatially FLIPPED via x_eval =
    %     L_total - x, so its trough sits at the interface and its peak at the
    %     outer edge x=L_total. (exact discrete mirror -- see
    %     build_atlas_sheets.m; no interpolation across grid cells.)
    %
    % NOT fabrication / not cross-cell interpolation: only the ONE stored per-cell
    % profile is evaluated at mesh points (via interp1 for HPSS/MPSS, constant for
    % HSS) -- exactly what get_pde_trace.m does with deval. If the sheet does not
    % exist at (i,j) (either half's existence_mask false), returns [] (Invariant 2:
    % missing cells stay unavailable).
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
    y = max(0, y);      % positivity guard (matches get_pde_trace.m)
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
