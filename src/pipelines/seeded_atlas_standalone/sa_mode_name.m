function varargout = sa_mode_name(action, varargin)
%SA_MODE_NAME  Mode names: the shape label of a coupled root (a quasi-steady
%   state), per half.
%
%   Labels: 0 = no-peak half; 1+ / 1- = half-peak at a pole (left / right
%   end); 2- = full-peak (peak in the middle of the half).
%
%   A root is named '<left>/<right>' from its main variable (sa_main_variable)
%   on the two halves. Each half is read in its DRAWN left-to-right
%   orientation on the full cell (left half: outer pole -> cell middle, the
%   'seam'; right half: seam -> outer pole). A half is
%     '0'   flat: (max - min) / ref < model_options.shape_flat_tol, where ref
%           is the main variable's value in the uniform no-peak (HSS) half-cell
%           state at the anchor ('reference');
%     'n+'  / 'n-' otherwise: n = the number of monotone runs, where a turn
%           counts only when the profile moves back by at least
%           shape_turn_frac x (max - min of the half) (hysteresis); '+' when
%           the half's drawn left end is high (the first run descends), '-'
%           when it is low (the first run rises).
%   So a peak at the half's left end is 1+, at its right end 1-, in its
%   middle 2-, at both ends 2+, 1.5 periods 3+ / 3-. Examples: 1+/0 a pole
%   peak in the left half, 0/1- one in the right half, 1-/1+ a single peak
%   at the cell middle, 2-/0 a middle peak in the left half.
%
%   The seam (the cell middle, x = L/2). Let mx be the largest value within
%   shape_edge_frac of the seam on either side (the seam window) and th
%   = shape_turn_frac x (max - min of the whole cell). The seam is a PEAK
%   when the cell is not flat and, on each side, scanning along the whole
%   cell outward from mx's own position, the profile falls to mx - th
%   before it ever exceeds mx + th/4 (or does neither): a local maximum of
%   the whole cell inside the seam window, wiggles smaller than th ignored,
%   a top up to th/4 higher just beyond the window tolerated
%   (is_seam_peak). Both halves are then named as they
%   are (x-/y+ at the seam). Otherwise a seam rise is the neighbour's tail,
%   not a peak: in a half whose last run (toward the seam) rises, that run
%   is a TAIL when the other half, scanned outward from the seam, rises
%   above the run's top by more than the seam mismatch |m_L - m_R| at the
%   seam before falling th_h below it (th_h the half's own turn threshold);
%   that is, the neighbour keeps rising across the seam. A tail is replaced
%   by its lowest value (the run's start) before the flatness and turn
%   tests, which then see only the rest of the half. Replacing (not
%   truncating) keeps the half's length, so the windows keep their meaning;
%   max, min and the runs of the rest are the same either way.
%
%   Everything is computed on the halves as stored (xi = 0 the outer pole,
%   the last point the seam, both halves); the right half's drawn name takes
%   its sign from its last stored run. Reflecting the cell swaps the stored
%   halves exactly, so the mirror image of a root has the label with the
%   halves swapped and the signs of odd modes flipped (1+/0 <-> 0/1-,
%   2-/1- <-> 1+/2-).
%
%   [lab, info] = SA_MODE_NAME('label', mL, mR, ref, mo)
%       mL, mR: the main variable per half as stored (rows); ref > 0;
%       mo: model_options (shape_flat_tol, shape_turn_frac,
%       shape_edge_frac). info: halves {left right}, seam_peak, tail
%       [left right], flat [left right], amplitude [left right] (max - min
%       after a tail is removed, divided by ref), ref.
%   [ref, src] = SA_MODE_NAME('reference', fs, index)
%       The uniform-state value of the main variable at the anchor: the sum
%       of fs.anchor_hss_values(index) (field_scale's library HSS values
%       interpolated bilinearly to the anchor itself), in the order given.
%       A field scale written by an older version carries only hss_values
%       (the HSS of the native cell nearest the anchor, raw, not floored);
%       they are used then, and src says so ('anchor' | 'nearest_cell').
%       Refused unless finite and > 0.
%   fam = SA_MODE_NAME('library_family', lab, types, map)
%       The library family 'A_B' (A, B in types: the model's library_types;
%       a family names a sheet) a root of name lab corresponds to, from the model's mode_families
%       map (map.left / map.right: n x 2 cells {half mode name, library
%       type}; both models: the flat half, the pole peak of mode 1 (left
%       1+, right 1-) and the middle peak 2- have a library type). A name
%       whose half is not in the map (a seam peak, 2+, 3+/-, ...) -> ''
%       (no library family), as is every mode name when map is absent. A
%       label that already is a library family name 'A_B' is returned as it
%       is (older folders, test stubs).
%   [n, s] = SA_MODE_NAME('parse', lab)
%       Mode numbers [nL nR] and signs {sL sR} ('' for 0); n = [] when lab
%       is not a mode name.
%   [v, why] = SA_MODE_NAME('bilinear', axis_row, axis_col, x, value)
%       Bilinear interpolation at x = [x_col x_row] on a native grid
%       (increasing axes): value(i, j) returns the column of values at row
%       i, column j, or [] where there is none. Cells of zero weight are
%       not read (a node needs only itself). v = [] when x lies outside the
%       grid or a needed cell has no value (why says which).
%   t = SA_MODE_NAME('file_tag', lab)
%       A file-name-safe form of a family name: '+' -> p, '-' -> m,
%       '/' -> _ (1+/0 -> 1p_0, 2-/1- -> 2m_1m); a name made of letters,
%       digits and '_' only (a library family name such as HPSS_HSS) is unchanged.

switch action
    case 'label',          [varargout{1}, varargout{2}] = label(varargin{:});
    case 'reference',      [varargout{1}, varargout{2}] = reference(varargin{:});
    case 'library_family', varargout{1} = library_family(varargin{:});
    case 'parse',          [varargout{1}, varargout{2}] = parse(varargin{:});
    case 'file_tag',       varargout{1} = file_tag(varargin{:});
    case 'bilinear',       [varargout{1}, varargout{2}] = bilinear(varargin{:});
    otherwise, error('SA:MODE_NAME', 'Unknown sa_mode_name action %s.', action);
end
end

% =========================================================================

function [lab, info] = label(mL, mR, ref, mo)
if ~(isnumeric(ref) && isscalar(ref) && isfinite(ref) && ref > 0)
    error('SA:MODE_NAME:E_REFERENCE', 'The uniform-state reference must be a finite positive scalar.');
end
mL = double(mL(:).'); mR = double(mR(:).');
w = [mL, fliplr(mR)];
th = mo.shape_turn_frac * (max(w) - min(w));
seam_peak = (max(w) - min(w)) / ref >= mo.shape_flat_tol && is_seam_peak(mL, mR, th, mo.shape_edge_frac);
tail = [false false];
if ~seam_peak
    gap = abs(mL(end) - mR(end));
    [tail(1), cL] = seam_tail(mL, mR, gap, mo.shape_turn_frac);
    [tail(2), cR] = seam_tail(mR, mL, gap, mo.shape_turn_frac);
    if tail(1), mL = cL; end
    if tail(2), mR = cR; end
end
[nL, fL, aL] = half_name(mL, ref, mo, 'left');
[nR, fR, aR] = half_name(mR, ref, mo, 'right');
lab = [nL '/' nR];
info = struct('halves', {{nL, nR}}, 'seam_peak', seam_peak, 'tail', tail, 'flat', [fL fR], ...
    'amplitude', [aL aR], 'ref', ref);
end

function tf = is_seam_peak(mL, mR, th, frac)
% The largest value of the seam window is a local maximum of the cell:
% scanned along the whole cell from that value's own position, neither side
% rises more than SEAM_MARGIN x th above it before falling th below it.
% Scanning from the maximum's position (not from the seam) keeps a steep
% flank that rises by more than th across the window from reading as a
% fall. The margin lets a broad middle peak whose top lies a few per cent of
% a half beyond the window, or a flat top tilted by round-off, still read as
% a peak at the seam; it stays small enough that a narrow peak just outside
% the window names its own half.
SEAM_MARGIN = 0.25;
eL = max(1, round(frac * numel(mL))); eR = max(1, round(frac * numel(mR)));
w = [mL, fliplr(mR)];
lo = numel(mL) - eL + 1; hi = numel(mL) + eR;
[mx, k] = max(w(lo:hi));
j = lo + k - 1;
tf = ~strcmp(scan(w(j - 1:-1:1), mx, SEAM_MARGIN * th, th), 'up') && ...
    ~strcmp(scan(w(j + 1:end), mx, SEAM_MARGIN * th, th), 'up');
end

function [tf, m] = seam_tail(m, other, gap, frac)
% the last run of m (toward the seam) rises and the other half keeps rising
% across the seam: that run is the other half's tail; m with it replaced by
% its lowest value
tf = false;
th = frac * (max(m) - min(m));
p = pivots(m, th);
if numel(p) < 2 || ~(m(p(end)) > m(p(end - 1))), return; end
if strcmp(scan(fliplr(other), m(p(end)), gap, th), 'up')
    tf = true;
    m(p(end - 1):end) = m(p(end - 1));
end
end

function s = scan(v, level, over, th)
% v ordered outward from the seam: 'up' if it exceeds level + over before it
% falls to level - th, 'down' if the fall comes first, 'none' if neither
k = find(v > level + over | v <= level - th, 1);
if isempty(k), s = 'none';
elseif v(k) > level + over, s = 'up';
else, s = 'down';
end
end

function [name, flat, a] = half_name(m, ref, mo, side)
a = (max(m) - min(m)) / ref;
flat = a < mo.shape_flat_tol;
if flat, name = '0'; return; end
p = pivots(m, mo.shape_turn_frac * (max(m) - min(m)));
n = numel(p) - 1;
if strcmp(side, 'left')
    high = m(p(2)) < m(p(1));            % the first run descends: left end high
else
    high = m(p(end)) > m(p(end - 1));    % stored last run rises = drawn first run descends
end
if high, s = '+'; else, s = '-'; end
name = sprintf('%d%s', n, s);
end

function p = pivots(m, th)
% Turning points with hysteresis th: p(1) the start of the first run, then
% every extremum after which the profile moves back by at least th, and the
% last extremum. numel(p) - 1 monotone runs; none when max - min < th.
p = zeros(1, 0);
if ~(th > 0), return; end
lo = 1; hi = 1; dir = 0; ext = 1;
for k = 2:numel(m)
    if dir == 0
        if m(k) > m(hi), hi = k; end
        if m(k) < m(lo), lo = k; end
        if m(hi) - m(lo) >= th
            if lo < hi, p = lo; dir = 1; ext = hi;
            else, p = hi; dir = -1; ext = lo;
            end
        end
    elseif dir == 1
        if m(k) > m(ext), ext = k;
        elseif m(ext) - m(k) >= th, p(end + 1) = ext; dir = -1; ext = k; %#ok<AGROW>
        end
    else
        if m(k) < m(ext), ext = k;
        elseif m(k) - m(ext) >= th, p(end + 1) = ext; dir = 1; ext = k; %#ok<AGROW>
        end
    end
end
if dir ~= 0, p(end + 1) = ext; end
end

% =========================================================================

function [ref, src] = reference(fs, index)
if isstruct(fs) && isfield(fs, 'anchor_hss_values')
    h = fs.anchor_hss_values(:).'; src = 'anchor';
elseif isstruct(fs) && isfield(fs, 'hss_values')
    h = fs.hss_values(:).'; src = 'nearest_cell';
else
    error('SA:MODE_NAME:E_REFERENCE', ['The field scale carries no uniform-state values: the ' ...
        'reference of the mode names is unavailable.']);
end
ref = h(index(1));
for k = 2:numel(index), ref = ref + h(index(k)); end
if ~(isfinite(ref) && ref > 0)
    error('SA:MODE_NAME:E_REFERENCE', ['The uniform-state value of the main variable is %g: the ' ...
        'flatness test needs a positive reference.'], ref);
end
end

function fam = library_family(lab, types, map)
fam = '';
t = regexp(lab, '^([A-Za-z]\w*)_([A-Za-z]\w*)$', 'tokens', 'once');
if ~isempty(t) && all(ismember(t, types)), fam = lab; return; end
if nargin < 3 || isempty(map) || isempty(parse(lab)), return; end
h = strsplit(lab, '/');
a = map.left(strcmp(map.left(:, 1), h{1}), 2); b = map.right(strcmp(map.right(:, 1), h{2}), 2);
if isempty(a) || isempty(b) || ~all(ismember([a b], types)), return; end
fam = [a{1} '_' b{1}];
end

function [n, s] = parse(lab)
n = []; s = {};
if ~(ischar(lab) && isrow(lab)), return; end
t = regexp(lab, '^(0|[1-9]\d*[+-])/(0|[1-9]\d*[+-])$', 'tokens', 'once');
if isempty(t), return; end
n = zeros(1, 2); s = {'', ''};
for k = 1:2
    if ~strcmp(t{k}, '0')
        n(k) = str2double(t{k}(1:end - 1)); s{k} = t{k}(end);
    end
end
end

function t = file_tag(lab)
t = regexprep(regexprep(regexprep(lab, '\+', 'p'), '-', 'm'), '/', '_');
t = regexprep(t, '[^A-Za-z0-9_]', '_');
end

function [v, why] = bilinear(ar, ac, x, value)
v = []; why = '';
[jc, tc] = bracket(ac, x(1)); [ir, tr] = bracket(ar, x(2));
if isempty(jc) || isempty(ir), why = 'outside the native grid'; return; end
cells = [ir jc (1 - tr) * (1 - tc); ir + 1 jc tr * (1 - tc); ir jc + 1 (1 - tr) * tc; ir + 1 jc + 1 tr * tc];
acc = 0;
for k = 1:4
    if cells(k, 3) == 0, continue; end
    y = value(cells(k, 1), cells(k, 2));
    if isempty(y), why = sprintf('no value at cell (%d,%d)', cells(k, 1), cells(k, 2)); return; end
    acc = acc + cells(k, 3) * y(:);
end
v = acc;
end

function [k, t] = bracket(a, x)
% the lower index k of the interval [a(k), a(k+1)] holding x and the weight t of a(k+1)
k = []; t = [];
if ~(x >= a(1) && x <= a(end)) || numel(a) < 2, return; end
k = min(find(a <= x, 1, 'last'), numel(a) - 1);
t = (x - a(k)) / (a(k + 1) - a(k));
end
