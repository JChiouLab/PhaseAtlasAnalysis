function varargout = build_mode_branch(varargin)
%BUILD_MODE_BRANCH  Mode-n seed library (milestone 8, packet section 3; fission yeast).
%
%   R = BUILD_MODE_BRANCH(p, out_dir, n)
%   R = BUILD_MODE_BRANCH(p, out_dir, n, opts)
%       Mode n (n >= 3) at every cell of the library in out_dir: a half-peak
%       PIECE on [0, l], l = L_half/n, at the cell's half masses / n (the same
%       densities), searched by run_quarter_hpss_sweeper_for_mpss UNCHANGED
%       with the transformed parameters p_n (L = 2L/n, M_BG = 2 M_BG/n,
%       M_P = 2 M_P/n, quarter_mass_tol = mass_tol/n; at n = 2 p_n equals p);
%       then the T-start stitch (n copies, a trough at the pole: drawn n-)
%       and the acceptance of packet section 3.1. p is the struct that
%       atlas_build_libraries passes to build_mpss_branch (the script's p
%       plus out_dir); build_mpss_branch's default fill is applied here
%       verbatim, so the base parameters are the stored MPSS par's field set.
%       Writes branch_MODE<n>_<id>.mat and branch_MODE<n>_latest.mat into
%       out_dir: the grid (PAK_DIFF, GEF_DIFF, res) COPIED from the sibling
%       branch_HSS_latest.mat (bit-equal to the HPSS and MPSS files, checked),
%       never the sweeper's axes (scaled by about 2/n); par = the original
%       parameters; piece_par = p_n. The sweeper's own outputs and
%       checkpoints go to out_dir/mode<n>_work/. An empty mask is written
%       with run_status 'empty'.
%       opts.max_reroot_tries   the sweeper's Phase 0 cap (default: the
%                               sweeper's own default res^2, unchanged)
%       opts.sweeper            injected sweeper (tests; default
%                               @run_quarter_hpss_sweeper_for_mpss)
%   p_n = BUILD_MODE_BRANCH('transform', p, n)
%   p = BUILD_MODE_BRANCH('defaults', p)        build_mpss_branch's default fill
%   C = BUILD_MODE_BRANCH('construct', sol_q, n, nG_half, nP_half, p, branch)
%       One cell's stitch and acceptance from a stored piece (zero solves):
%       branch 'legacy' (n = 2 only: build_mpss_branch's reflect_stitch and
%       classify_mpss_profile verbatim) or 'general' (the n-general code).
%       C: ok, why, x, y, amp, peak_x (the stitch's argmax position, as the
%       MPSS file stores it), xstar_frac (the piece's argmax position / l),
%       mass_res, seam_slope, trough_slope, sym_res, n_dense.
%   [x, y] = BUILD_MODE_BRANCH('stitch', xq, yq5, n, start)   start 'T' | 'P'
%
%   The acceptance, in the legacy order: negativity; mesh monotonicity; the
%   piece's argmax position x*/l < 0.2 (the legacy peak-fraction test at
%   n = 2: pk_frac = (l - x*)/2l in (0.4, 0.6)); amplitude >= 0.5; seam
%   slope at the piece's peak end <= slope_tol (as legacy); n >= 3 only: the
%   trough-end slope <= slope_tol (an interior join); stitched mass within
%   mass_tol, densified to n x 2,000 points (legacy 4,000 at n = 2);
%   symmetry: even n reflection about the half's centre (legacy), odd n the
%   T-start stitch against the reflected P-start stitch of the same piece.
%   The piece mass gate inside the sweeper is mass_tol/n.

if ischar(varargin{1})
    switch varargin{1}
        case 'transform', varargout{1} = transform(varargin{2:end});
        case 'defaults',  varargout{1} = mpss_defaults(varargin{2:end});
        case 'construct', varargout{1} = construct(varargin{2:end});
        case 'stitch',    [varargout{1}, varargout{2}] = stitch(varargin{2:end});
        otherwise, error('BUILD_MODE_BRANCH:ACTION', 'Unknown action %s.', varargin{1});
    end
    return
end
varargout{1} = build(varargin{:});
end

% =========================================================================

function R = build(p, out_dir, n, opts)
if nargin < 4 || isempty(opts), opts = struct(); end
if ~(isscalar(n) && n == fix(n) && n >= 3)
    error('BUILD_MODE_BRANCH:N', 'Mode n must be an integer >= 3 (mode 2 is the MPSS library).');
end
sweeper = @run_quarter_hpss_sweeper_for_mpss;
if isfield(opts, 'sweeper') && ~isempty(opts.sweeper), sweeper = opts.sweeper; end
t0 = tic;
p = mpss_defaults(p);
p_n = transform(p, n);
if isfield(opts, 'max_reroot_tries') && ~isempty(opts.max_reroot_tries)
    p_n.max_reroot_tries = opts.max_reroot_tries;
end
% ---- the library grid: the sibling files, bit-equal --------------------
sib = {'HSS', 'HPSS', 'MPSS'};
G = cell(1, 3);
for k = 1:3
    f = fullfile(out_dir, sprintf('branch_%s_latest.mat', sib{k}));
    if exist(f, 'file') ~= 2
        error('BUILD_MODE_BRANCH:SIBLING', 'Mode files need the three library files; missing %s.', f);
    end
    G{k} = load(f, 'PAK_DIFF', 'GEF_DIFF', 'res');
end
if ~isequal(G{1}, G{2}) || ~isequal(G{1}, G{3})
    error('BUILD_MODE_BRANCH:GRID', 'The three library files in %s do not share one grid.', out_dir);
end
PAK_DIFF = G{1}.PAK_DIFF; GEF_DIFF = G{1}.GEF_DIFF; res = G{1}.res;
if res ~= p.res
    error('BUILD_MODE_BRANCH:GRID', 'The library grid has res %d, the parameters say %d.', res, p.res);
end
% ---- the piece search (the sweeper, unchanged) ---------------------------
work = fullfile(out_dir, sprintf('mode%d_work', n));
if ~exist(work, 'dir'), mkdir(work); end
qres = sweeper(p_n, work);
src_file = fullfile(work, sprintf('quarter_hpss_for_mpss_%s.mat', qres.computation_id));
% ---- stitch and acceptance per cell ---------------------------------------
L_half = p.L / 2;
ell = L_half / n;
n_GEF_bar = p.M_BG / 2;
n_PAK_bar = p.M_P / 2;
existence_mask = false(res, res);
SOL = cell(res, res);
SOL_piece = cell(res, res);
fail_reason = strings(res, res);
acc = struct('amp', nan(res), 'peak_x', nan(res), 'xstar_frac', nan(res), 'mass_res', nan(res), ...
    'seam_slope', nan(res), 'trough_slope', nan(res), 'sym_res', nan(res), 'n_dense', n * 2000);
for lin = 1:res*res
    [ii, jj] = ind2sub([res res], lin);
    if qres.failure_mask(ii, jj)
        fail_reason(ii, jj) = "piece_sweep_failed";
        continue
    end
    sol_q = qres.SOL_grid{ii, jj};
    nG_half = n_GEF_bar + GEF_DIFF(ii, jj) / 2;
    nP_half = n_PAK_bar + PAK_DIFF(ii, jj) / 2;
    C = construct(sol_q, n, nG_half, nP_half, p, 'general');
    for f = {'amp', 'peak_x', 'xstar_frac', 'mass_res', 'seam_slope', 'trough_slope', 'sym_res'}
        acc.(f{1})(ii, jj) = C.(f{1});
    end
    if ~C.ok
        fail_reason(ii, jj) = C.why;
        continue
    end
    existence_mask(ii, jj) = true;
    SOL{ii, jj} = struct('x', C.x, 'y', C.y);
    SOL_piece{ii, jj} = struct('x', sol_q.x(:).', 'y', sol_q.y([1, 3, 5, 7, 9], :));
end
% ---- write ------------------------------------------------------------------
branch_type = sprintf('MODE%d', n);
construction = 'reflected_piece_HPSS_n';
mode_n = n;
piece_length = ell;
par = p;
piece_par = p_n;
acceptance = acc;
if any(existence_mask(:)), run_status = 'complete'; else, run_status = 'empty'; end
computation_id = qres.computation_id;
source_sweep_file = src_file;
source_sweep_sha256 = '';
if exist(src_file, 'file') == 2, source_sweep_sha256 = sha256_file(src_file); end
builder_sha256 = sha256_file([mfilename('fullpath') '.m']);
wall_time = toc(t0);
out_path = fullfile(out_dir, sprintf('branch_MODE%d_%s.mat', n, computation_id));
save(out_path, 'branch_type', 'construction', 'mode_n', 'piece_length', 'PAK_DIFF', 'GEF_DIFF', ...
    'res', 'par', 'piece_par', 'existence_mask', 'fail_reason', 'run_status', 'SOL_piece', 'SOL', ...
    'acceptance', 'source_sweep_file', 'source_sweep_sha256', 'builder_sha256', ...
    'computation_id', 'wall_time');
latest = fullfile(out_dir, sprintf('branch_MODE%d_latest.mat', n));
copyfile(out_path, latest);
reasons = unique(fail_reason(fail_reason ~= ""));
counts = arrayfun(@(r) sum(fail_reason(:) == r), reasons);
fprintf('[MODE%d] %d/%d cells (piece l = %.6g, run_status %s), wall %.1f s; wrote %s\n', n, ...
    nnz(existence_mask), res * res, ell, run_status, wall_time, latest);
for k = 1:numel(reasons), fprintf('   fail %-22s %d\n', reasons(k), counts(k)); end
R = struct('path', latest, 'id_path', out_path, 'n', n, 'n_valid', nnz(existence_mask), ...
    'run_status', run_status, 'wall_time', wall_time, 'ell', ell, 'fail_reasons', {cellstr(reasons)}, ...
    'fail_counts', counts, 'n_pieces', nnz(~qres.failure_mask));
end

function p = mpss_defaults(p)
% build_mpss_branch.m:60-79 verbatim (the default fill before its sweep call)
if ~isfield(p, 'res'),            p.res = 3;      end
if ~isfield(p, 'pde_T_end'),      p.pde_T_end   = 2e6;  end
if ~isfield(p, 'pde_T_chunk'),    p.pde_T_chunk = 1e6;  end
if ~isfield(p, 'pde_max_chunks'), p.pde_max_chunks = 12; end
if ~isfield(p, 'slope_tol'),      p.slope_tol = 1e-2;   end  % C1 seam-slope tol
if ~isfield(p, 'sym_tol'),        p.sym_tol   = 1e-6;   end  % symmetry tol
if ~isfield(p, 'mass_tol'),       p.mass_tol  = 1e-4;   end  % mass-doubling tol
if ~isfield(p, 'max_res_warn'),   p.max_res_warn = 8;   end  % serial-mode guideline
if ~isfield(p, 'use_parallel'),      p.use_parallel = true;  end
if ~isfield(p, 'checkpoint_every'),  p.checkpoint_every = 5; end  % cells per checkpoint
if ~isfield(p, 'cell_timeout'),      p.cell_timeout = 60;    end  % seconds, parallel only
if ~isfield(p, 'resume_from'),       p.resume_from = '';     end  % checkpoint_mpss_*.mat path
if ~isfield(p, 'accept_partial'),    p.accept_partial = false; end
if ~isfield(p, 'max_cells_this_run'), p.max_cells_this_run = Inf; end
if ~isfield(p, 'mpss_strategy'),       p.mpss_strategy = 'continuation'; end
if ~isfield(p, 'quarter_resume_from'), p.quarter_resume_from = '';       end
end

function p_n = transform(p, n)
% The spec's factor-of-2 route with 2 replaced by n: the sweeper solves on
% p_n.L/4 = L/(2n) at (p_n.M/2 + diff/2)/2 = (half masses)/n. At n = 2 the
% arithmetic is exact (x2 then /2) and the sweeper's own default
% quarter_mass_tol = mass_tol/2 is the same double, so p_n is p itself.
p_n = p;
p_n.L = 2 * p.L / n;
p_n.M_BG = 2 * p.M_BG / n;
p_n.M_P = 2 * p.M_P / n;
if n ~= 2
    p_n.quarter_mass_tol = p.mass_tol / n;
end
end

% =========================================================================
% stitch and acceptance

function C = construct(sol_q, n, nG_half, nP_half, p, branch)
C = struct('ok', false, 'why', "", 'x', [], 'y', [], 'amp', NaN, 'peak_x', NaN, ...
    'xstar_frac', NaN, 'mass_res', NaN, 'seam_slope', NaN, 'trough_slope', NaN, 'sym_res', NaN, ...
    'n_dense', n * 2000);
switch branch
    case 'legacy'
        if n ~= 2, error('BUILD_MODE_BRANCH:LEGACY', 'The legacy expressions exist for n = 2 only.'); end
        [xm, ym, seam_du] = reflect_stitch(sol_q);
        [ok, why, this_amp, pk_x, mass_res, sym_res] = classify_mpss_profile(xm, ym, seam_du, ...
            nG_half, nP_half, p);
        C.ok = ok; C.why = why; C.x = xm; C.y = ym; C.amp = this_amp; C.peak_x = pk_x;
        C.mass_res = mass_res; C.seam_slope = seam_du; C.sym_res = sym_res;
        x_q = sol_q.x(:).';
        [~, k] = max(sol_q.y(1, :));
        C.xstar_frac = x_q(k) / x_q(end);
    case 'general'
        C = classify_general(sol_q, n, nG_half, nP_half, p, C);
    otherwise
        error('BUILD_MODE_BRANCH:BRANCH', 'Unknown branch %s.', branch);
end
end

function C = classify_general(sol_q, n, nG_half, nP_half, p, C)
x_q = sol_q.x(:).';
yq5 = sol_q.y([1, 3, 5, 7, 9], :);
ell = x_q(end);
[xm, ym] = stitch(x_q, yq5, n, 'T');
C.x = xm; C.y = ym;
u = ym(1, :);
[~, ipk] = max(u);
C.peak_x = xm(ipk);
C.amp = max(u) - min(u);
[~, k] = max(yq5(1, :));
C.xstar_frac = x_q(k) / ell;
% seam slope at the piece's peak end (legacy rule), trough end (n >= 3)
has_d = size(sol_q.y, 1) >= 2 && any(sol_q.y(2, :) ~= 0);
if has_d
    C.seam_slope = abs(sol_q.y(2, 1));
    C.trough_slope = abs(sol_q.y(2, end));
else
    C.seam_slope = abs((yq5(1, 2) - yq5(1, 1)) / (x_q(2) - x_q(1)));
    C.trough_slope = abs((yq5(1, end) - yq5(1, end-1)) / (x_q(end) - x_q(end-1)));
end
% stitched mass, densified to n x 2,000 points (legacy 4,000 at n = 2)
xf = linspace(xm(1), xm(end), C.n_dense);
uf_ = interp1(xm, ym(1,:), xf, 'pchip');   qf_ = interp1(xm, ym(2,:), xf, 'pchip');
vf_ = interp1(xm, ym(3,:), xf, 'pchip');   wf_ = interp1(xm, ym(4,:), xf, 'pchip');
pf_ = interp1(xm, ym(5,:), xf, 'pchip');
gef_mass = trapz(xf, uf_ + qf_ + vf_ + wf_);
pak_mass = trapz(xf, pf_ + qf_);
C.mass_res = max(abs(gef_mass - nG_half), abs(pak_mass - nP_half));
% symmetry
if mod(n, 2) == 0
    uf = interp1(xm, u, xm(end) - xm, 'pchip');            % about the half's centre (legacy)
else
    [xP, yP] = stitch(x_q, yq5, n, 'P');                   % odd n: T-start vs reflected P-start
    uf = interp1(xP, yP(1, :), xm(end) - xm, 'pchip');
end
C.sym_res = max(abs(u - uf));
% acceptance, in the legacy order (the trough-end gate after the seam gate)
if any(ym < -1e-9, 'all')
    C.why = "negative_conc";
elseif ~all(diff(xm) > 0)
    C.why = "nonmonotone_mesh";
elseif ~(C.xstar_frac < 0.2)
    C.why = "peak_not_interior";
elseif C.amp < 0.5
    C.why = "amp_below_0p5";
elseif C.seam_slope > p.slope_tol
    C.why = "kinked_seam";
elseif n >= 3 && C.trough_slope > p.slope_tol
    C.why = "kinked_trough_join";
elseif C.mass_res > p.mass_tol
    C.why = "mass_doubling_fail";
elseif C.sym_res > p.sym_tol
    C.why = "symmetry_fail";
else
    C.ok = true;
end
end

function [x, y] = stitch(xq, yq, n, start)
% n copies of the piece (peak at its left end) on [0, n l], l = xq(end):
% T-start (n-, a trough at the pole) or P-start (n+). A copy whose right end
% is a peak is the mirrored piece, x = (k+1) l - fliplr(xq); otherwise
% x = k l + xq; each copy's last node (the duplicate join) is dropped, the
% final copy's kept. At n = 2, 'T' is reflect_stitch verbatim.
xq = xq(:).';
l = xq(end);
x = zeros(1, 0); y = zeros(size(yq, 1), 0);
for k = 0:n-1
    if strcmp(start, 'T'), mirrored = mod(k, 2) == 0; else, mirrored = mod(k, 2) == 1; end
    if mirrored
        xk = (k + 1) * l - fliplr(xq); yk = fliplr(yq);
    else
        xk = k * l + xq; yk = yq;
    end
    if k < n - 1
        xk = xk(1:end-1); yk = yk(:, 1:end-1);
    end
    x = [x, xk]; %#ok<AGROW>
    y = [y, yk]; %#ok<AGROW>
end
end

% ---- build_mpss_branch.m:625-711, verbatim (the n = 2 legacy branch) -------

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
    ok = false; why = "";
    u = ym(1,:);
    [~, ipk] = max(u);
    pk_x = xm(ipk);
    pk_frac = pk_x / xm(end);
    this_amp = max(u) - min(u);
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

function h = sha256_file(f)
fid = fopen(f, 'r'); b = fread(fid, Inf, '*uint8'); fclose(fid);
md = java.security.MessageDigest.getInstance('SHA-256');
if ~isempty(b), md.update(typecast(b(:), 'int8')); end
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
end
