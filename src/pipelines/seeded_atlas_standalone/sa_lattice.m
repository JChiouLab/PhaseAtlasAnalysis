function varargout = sa_lattice(action, varargin)
%SA_LATTICE  The lattice (grid of the phase plane) in coordinates that are
%   the differences between the two halves' amounts of each conserved protein.
%
%   L = SA_LATTICE('build', M, theta, S)
%       Derive the lattice from the settings and the model's half-mass targets
%       (each half's target amount of a conserved protein, (M_c +/- d_c)/2 for
%       total amount M_c; the units of d are those of the model's pool_targets).
%       One coordinate per conserved protein. A coordinate is
%           d = anchor + step*k,          k integer, NEVER accumulated,
%       where the anchor is the lattice origin (default d = 0), and a node
%       (grid point) is admissible when every half-mass target the model
%       derives at that coordinate is at least `lattice.margin`. The extents are the
%       admissibility boundary stepped outward from the anchor on each axis,
%       not a library property: whether starting profiles exist is a separate,
%       later question.
%
%   tf = SA_LATTICE('admissible', M, theta, S, d)
%       The admissibility predicate on one coordinate row.
%
%   [S, SRC] = SA_LATTICE('derive_step', M, theta, S, SRC)
%       lattice.points_per_axis -> the per-axis step (a no-op when it is
%       empty). Per axis c, with n_c the requested odd count:
%         R_c   the admissible half-width from the anchor, by bisection (200
%               halvings) of 'admissible' in each direction; the two
%               directions must agree to 1e-12 relative, else SA:LATTICE:E_STEP
%         cap   lattice.max_abs_coordinate(c) < R_c: the CLOSED interval
%               [-cap, cap], step_c = 2 cap / (n_c - 1) (end nodes on the cap)
%         else  the OPEN interval (-R_c, R_c), step_c = 2 R_c / n_c (the
%               outermost node half a step inside the admissible boundary)
%       S.lattice.step becomes the 1 x n_pools row (SRC 'derived'); the
%       walk and the cap must then give exactly n_c nodes on every axis,
%       else SA:LATTICE:E_STEP. Called by seeded_atlas_run before the run
%       folder is prepared, and by sa_resume('status'), so a refused
%       derivation writes nothing.
%   h = SA_LATTICE('steps', L)
%       The per-axis step row of a lattice (L.steps; a lattice stored by an older
%       version carries the scalar L.step, expanded here).
%
%   Returned fields
%     step, anchor, margin        the settings actually used (step as set:
%                                 scalar or one per axis)
%     steps                       the step of every axis (1 x n_coordinates)
%     n_coordinates               one per conserved protein of the model
%     index_min, index_max        integer extents per axis
%     index                       n_nodes x n_coordinates integer indices
%     coordinate                  n_nodes x n_coordinates, anchor + step*index
%     half_targets                n_nodes x n_half_targets model half-mass targets
%     n_nodes                     numel of the rectangle
%     key                         n_nodes x 1 cellstr target keys
%     rule                        the text of the rule applied

switch action
    case 'build',      varargout{1} = build(varargin{:});
    case 'admissible', varargout{1} = admissible(varargin{:});
    case 'target_key', varargout{1} = target_key(varargin{:});
    case 'derive_step', [varargout{1}, varargout{2}] = derive_step(varargin{:});
    case 'steps',      varargout{1} = steps_of(varargin{:});
    otherwise
        error('SA:LATTICE:E_ACTION', 'Unknown sa_lattice action %s.', action);
end
end

% =========================================================================

function L = build(M, theta, S)
nc     = numel(M.pools);
anchor = S.lattice.anchor;
step   = S.lattice.step;
steps  = expand(step, nc);
cap    = 1e6;                       % walk guard: a decoupled axis cannot reach it

if ~admissible(M, theta, S, anchor)
    error('SA:LATTICE:E_ANCHOR', ...
        ['The anchor coordinate is not admissible at margin %g: every half-pool ' ...
         'target must be at least the margin.'], S.lattice.margin);
end

kmin = zeros(1, nc); kmax = zeros(1, nc);
for c = 1:nc
    kmax(c) = walk(M, theta, S, anchor, c, +1, cap, steps);
    kmin(c) = -walk(M, theta, S, anchor, c, -1, cap, steps);
end

% the full rectangle, ascending on every axis (axis 1 slowest)
sizes = kmax - kmin + 1;
n     = prod(sizes);
idx   = zeros(n, nc);
subs  = cell(1, nc);
[subs{1:nc}] = ndgrid_desc(sizes);
for c = 1:nc
    idx(:, c) = subs{c}(:) + kmin(c) - 1;
end

coord = zeros(n, nc);
for c = 1:nc
    coord(:, c) = anchor(c) + steps(c) * idx(:, c);     % never accumulated
end

% every node of the rectangle must be admissible; fail closed if not
probe = M.pool_targets(theta, coord(1, :));
nt    = numel(probe.values);
ht    = zeros(n, nt);
for r = 1:n
    t = M.pool_targets(theta, coord(r, :));
    ht(r, :) = t.values;
end
bad = find(min(ht, [], 2) < S.lattice.margin);
if ~isempty(bad)
    error('SA:LATTICE:E_RECTANGLE', ...
        ['%d node(s) of the admissible rectangle are not admissible (first at ' ...
         'index [%s]). The per-axis walk assumes the model''s half-pool targets ' ...
         'depend on one coordinate each; this model does not satisfy that, so ' ...
         'the extents rule must be redesigned rather than relaxed.'], ...
        numel(bad), num2str(idx(bad(1), :)));
end

L = struct();
L.step          = step;
L.steps         = steps;
L.anchor        = anchor;
L.margin        = S.lattice.margin;
L.n_coordinates = nc;
L.index_min     = kmin;
L.index_max     = kmax;
L.index         = idx;
L.coordinate    = coord;
L.half_targets  = ht;
L.half_target_labels = probe.labels;
L.n_nodes       = n;
L.key           = cell(n, 1);
for r = 1:n
    L.key{r} = target_key(idx(r, :));
end
L.rule = ['d = anchor + step*k (never accumulated); the extents are the ' ...
    'admissibility boundary stepped outward from the anchor on each axis; ' ...
    'admissible means every model half-pool target >= margin'];
if ~isscalar(step)
    L.rule = [L.rule sprintf('; per-axis step %s (d_c = anchor_c + step_c*k_c)', mat2str(step, 17))];
end
end

function k = walk(M, theta, S, anchor, c, sgn, cap, steps)
% steps: the per-axis step row (as S is passed only for the margin)
k = 0;
while true
    d = anchor;
    d(c) = anchor(c) + steps(c) * sgn * (k + 1);
    if ~admissible(M, theta, S, d), return; end
    k = k + 1;
    if k > cap
        error('SA:LATTICE:E_WALK', ...
            ['The admissibility walk on axis %d passed %d steps without ' ...
             'leaving the admissible set.'], c, cap);
    end
end
end

function tf = admissible(M, theta, S, d)
t  = M.pool_targets(theta, d);
v  = t.values;
tf = all(isfinite(v)) && min(v) >= S.lattice.margin;
end

function h = expand(step, nc)
if isscalar(step), h = repmat(step, 1, nc); else, h = step; end
end

function h = steps_of(L)
if isfield(L, 'steps') && ~isempty(L.steps)
    h = L.steps;
else
    h = expand(L.step, L.n_coordinates);
end
end

function [S, SRC] = derive_step(M, theta, S, SRC)
if nargin < 4, SRC = struct(); end
if ~isfield(S.lattice, 'points_per_axis') || isempty(S.lattice.points_per_axis), return; end
id = 'SA:LATTICE:E_STEP';
nc = numel(M.pools);
n = expand(S.lattice.points_per_axis, nc);
anchor = S.lattice.anchor;
if ~admissible(M, theta, S, anchor)
    error('SA:LATTICE:E_ANCHOR', ...
        ['The anchor coordinate is not admissible at margin %g: every half-pool ' ...
         'target must be at least the margin.'], S.lattice.margin);
end
cap = expand(S.lattice.max_abs_coordinate, nc);
h = zeros(1, nc);
for c = 1:nc
    Rp = half_width(M, theta, S, anchor, c, +1);
    Rm = half_width(M, theta, S, anchor, c, -1);
    if abs(Rp - Rm) > 1e-12 * max(Rp, Rm)
        error(id, ['lattice.points_per_axis: axis %d is not symmetric about the anchor (admissible ' ...
            'half-widths %.17g and %.17g); the per-axis step needs symmetric pools.'], c, Rp, Rm);
    end
    R = min(Rp, Rm);
    if cap(c) < R
        h(c) = 2 * cap(c) / (n(c) - 1);
    else
        h(c) = 2 * R / n(c);
    end
end
S.lattice.step = h;
if isstruct(SRC) && isfield(SRC, 'lattice') && isstruct(SRC.lattice)
    SRC.lattice.step = 'derived';
end
% the walk and the cap must give exactly n_c nodes on every axis
for c = 1:nc
    kp = walk(M, theta, S, anchor, c, +1, 1e6, h);
    km = walk(M, theta, S, anchor, c, -1, 1e6, h);
    k = -km:kp;
    keep = abs(anchor(c) + h(c) * k) <= cap(c) + 1e-9 * h(c);
    if sum(keep) ~= n(c)
        error(id, ['lattice.points_per_axis: axis %d gives %d lattice nodes, not the %d requested ' ...
            '(step %.17g); internal error, report it.'], c, sum(keep), n(c), h(c));
    end
end
end

function R = half_width(M, theta, S, anchor, c, sgn)
% the admissible half-width from the anchor along axis c: bisection between
% an admissible distance lo and an inadmissible distance hi
lo = 0; hi = 1;
d = anchor;
while true
    d(c) = anchor(c) + sgn * hi;
    if ~admissible(M, theta, S, d), break; end
    lo = hi; hi = 2 * hi;
    if hi > 1e6
        error('SA:LATTICE:E_STEP', ['lattice.points_per_axis: axis %d stays admissible beyond %g; ' ...
            'the pool has no finite range.'], c, 1e6);
    end
end
for k = 1:200
    mid = (lo + hi) / 2;
    d(c) = anchor(c) + sgn * mid;
    if admissible(M, theta, S, d), lo = mid; else, hi = mid; end
end
R = lo;
end

function key = target_key(index)
parts = arrayfun(@(v) sprintf('%d', v), index, 'UniformOutput', false);
key = ['k:' strjoin(parts, ',')];
end

function varargout = ndgrid_desc(sizes)
% ndgrid over 1:sizes(c) with axis 1 varying SLOWEST, so the node order is
% ascending lexicographic in the integer indices.
nc = numel(sizes);
args = arrayfun(@(n) 1:n, sizes, 'UniformOutput', false);
out = cell(1, nc);
[out{nc:-1:1}] = ndgrid(args{nc:-1:1});
varargout = out;
end
