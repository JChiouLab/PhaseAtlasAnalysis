function varargout = sa_modes(action, varargin)
%SA_MODES  User modes, layered-settings helpers and target geometry
%   (model-free).
%
%   Lay = SA_MODES('layer', mode)       the mode layer (model-free leaves)
%   P = SA_MODES('layer_leaves', Lay)   leaf paths a layer struct sets
%   P = SA_MODES('mode_leaves')         union of every mode layer's leaves
%   names = SA_MODES('modes')           {'smoke','pilot','production','focus'}
%
%   Mode layer:
%                     smoke            pilot                 production
%     targets.policy  central_nine     pilot_sample          per_seed (one cold
%                                                            call per library
%                                                            seed; set
%                                                            all_seedable for
%                                                            every seedable node)
%     routes          cold, continuation, rescue on in every mode
%     continuation    region box,      region cross, box []  region box, box []
%                     box [1 1]
%     rescue.rungs    {secant}         {secant}              {secant}
%     rescue.ring     false            false                 false
%     budget          100 / 3,600 s    500 / 7,200 s         REQUIRED from the user
%     budget.wall_s   Inf              Inf                   Inf
%     pilot           -                sample_per_axis 7     -
%
%   focus (selected roots of an earlier run as starting profiles, restricted to
%   selected families (sheets)): the
%   production layer with routes.cold false (import only). The resolver
%   (sa_settings) adds the focus rules that depend on the other settings:
%   import.source_dir required (refused without), budgets required from
%   the user (as production), and, unless set by the user,
%   import.source_edges 'auto' and continuation.families 'selected' (both
%   source 'mode'). The focus layer sets the same leaves as production, so
%   mode_leaves (and from_run) are unchanged by it.
%
%   L = SA_MODES('cap_lattice', L, S)
%       lattice.max_abs_coordinate: a per-axis cap on
%       |coordinate| in lattice coordinate units (scalar = every axis; Inf =
%       none). The lattice stays the rectangle of sa_lattice, cut to
%       |anchor_c + step_c*k_c| <= cap_c (+ 1e-9 step_c); fails closed when the
%       cap removes the anchor.
%   q = SA_MODES('central_node', L)
%       The node nearest the symmetric point d = 0 (equal halves) by the
%       Euclidean coordinate distance; ties (within 1e-12 relative) by the
%       smaller Euclidean index distance to the anchor, then lattice order.
%   [nodes, missing] = SA_MODES('central_nine', L)
%       The central node and its eight lattice neighbours (|k - k_c| <= 1 on
%       every axis), lattice order; missing = neighbour indices that are not
%       lattice nodes (reported, never backfilled).
%   nodes = SA_MODES('pilot_sample', L, C, S)
%       pilot.sample_per_axis values evenly spaced over the native (library)
%       range of each axis (the extent of C.in_range), snapped to the
%       nearest lattice index, every combination, de-duplicated, central node
%       included; lattice order.
%   tf = SA_MODES('cross', L)
%       The continuation region 'cross': the lattice lines through the
%       central node (nodes differing from it on at most one axis), out to
%       the lattice extent, ring included.
%   h = SA_MODES('fingerprint', guess, coordinate, solver)
%       Report-only guess fingerprint: SHA-256 (hex)
%       of the IEEE bytes of guess.x, guess.y, the target coordinate and the
%       solver options (reltol, abstol, nmax, mesh_reject); '' when the guess
%       carries no x/y mesh. Two calls with equal fingerprints are repeated
%       calls.

switch action
    case 'layer',        varargout{1} = layer(varargin{:});
    case 'layer_leaves', varargout{1} = leaves_of(varargin{1}, '');
    case 'mode_leaves',  varargout{1} = mode_leaves();
    case 'modes',        varargout{1} = {'smoke', 'pilot', 'production', 'focus'};
    case 'cap_lattice',  varargout{1} = cap_lattice(varargin{:});
    case 'central_node', varargout{1} = central_node(varargin{:});
    case 'central_nine', [varargout{1}, varargout{2}] = central_nine(varargin{:});
    case 'pilot_sample', varargout{1} = pilot_sample(varargin{:});
    case 'cross',        varargout{1} = cross(varargin{:});
    case 'fingerprint',  varargout{1} = fingerprint(varargin{:});
    otherwise, error('SA:MODES', 'Unknown sa_modes action %s.', action);
end
end

% =========================================================================
% modes

function Lay = layer(mode)
routes = struct('cold', true, 'continuation', true, 'rescue', true);
% Rescue trimmed to secant, no ring rescue, in
% every mode (package defaults keep the full ladder and ring true)
rescue = struct('rungs', {{'secant'}}, 'ring', false);
switch mode
    case 'smoke'
        Lay = struct('targets', struct('policy', 'central_nine'), 'routes', routes, ...
            'continuation', struct('region', 'box', 'box', [1 1]), ...
            'rescue', rescue, ...
            'budget', struct('max_calls', 100, 'max_wall_s', 3600, 'wall_s', Inf));
    case 'pilot'
        Lay = struct('targets', struct('policy', 'pilot_sample'), 'routes', routes, ...
            'continuation', struct('region', 'cross', 'box', []), ...
            'rescue', rescue, ...
            'budget', struct('max_calls', 500, 'max_wall_s', 7200, 'wall_s', Inf), ...
            'pilot', struct('sample_per_axis', 7));
    case 'production'
        Lay = struct('targets', struct('policy', 'per_seed'), 'routes', routes, ...
            'continuation', struct('region', 'box', 'box', []), ...
            'rescue', rescue, ...
            'budget', struct('wall_s', Inf));
    case 'focus'
        Lay = struct('targets', struct('policy', 'per_seed'), ...
            'routes', struct('cold', false, 'continuation', true, 'rescue', true), ...
            'continuation', struct('region', 'box', 'box', []), ...
            'rescue', rescue, ...
            'budget', struct('wall_s', Inf));
    otherwise
        error('SA:SETTINGS:E_MODE', 'mode must be one of smoke, pilot, production, focus (got %s).', ...
            char(string(mode)));
end
end

function P = mode_leaves()
P = {};
for m = {'smoke', 'pilot', 'production', 'focus'}
    P = union(P, leaves_of(layer(m{1}), ''));
end
end

function P = leaves_of(s, path)
P = {};
f = fieldnames(s);
for k = 1:numel(f)
    if isempty(path), p = f{k}; else, p = [path '.' f{k}]; end
    v = s.(f{k});
    if isstruct(v) && isscalar(v) && ~isempty(fieldnames(v))
        P = [P, leaves_of(v, p)]; %#ok<AGROW>
    else
        P{end+1} = p; %#ok<AGROW>
    end
end
end

% =========================================================================
% lattice extent and target geometry

function L = cap_lattice(L, S)
cap = S.lattice.max_abs_coordinate;
nc = L.n_coordinates;
if isscalar(cap), cap = repmat(cap, 1, nc); end
if all(isinf(cap)), return; end
h = sa_lattice('steps', L);
keep = true(L.n_nodes, 1);
for c = 1:nc
    keep = keep & abs(L.coordinate(:, c)) <= cap(c) + 1e-9 * h(c);
end
a = find(all(L.index == 0, 2), 1);
if isempty(a) || ~keep(a)
    error('SA:LATTICE:E_CAP', ['lattice.max_abs_coordinate %s removes the anchor %s; the cap ' ...
        'must keep the anchor node.'], mat2str(cap), mat2str(L.anchor));
end
L.index = L.index(keep, :);
L.coordinate = L.coordinate(keep, :);
L.half_targets = L.half_targets(keep, :);
L.key = L.key(keep);
L.n_nodes = sum(keep);
L.index_min = min(L.index, [], 1);
L.index_max = max(L.index, [], 1);
if prod(L.index_max - L.index_min + 1) ~= L.n_nodes
    error('SA:LATTICE:E_CAP', 'the capped lattice is not a rectangle (internal error).');
end
L.cap = cap;
L.rule = [L.rule sprintf('; per-axis cap |d_c| <= %s (lattice.max_abs_coordinate)', mat2str(cap))];
end

function q = central_node(L)
d = sqrt(sum(L.coordinate.^2, 2));
dmin = min(d);
tol = 1e-12 * max(1, max(abs(L.coordinate(:))));
cand = find(d <= dmin + tol);
di = sqrt(sum(L.index(cand, :).^2, 2));
[~, o] = sortrows([di cand]);
q = cand(o(1));
end

function [nodes, missing] = central_nine(L)
q = central_node(L);
kc = L.index(q, :);
nc = size(L.index, 2);
nodes = find(all(abs(L.index - kc) <= 1, 2));
want = 3^nc;
missing = zeros(0, nc);
if numel(nodes) < want
    g = cell(1, nc);
    [g{:}] = ndgrid(-1:1);
    off = cell2mat(cellfun(@(x) x(:), g, 'UniformOutput', false));
    for r = 1:size(off, 1)
        k = kc + off(r, :);
        if ~any(all(L.index == k, 2)), missing(end+1, :) = k; end %#ok<AGROW>
    end
end
end

function nodes = pilot_sample(L, C, S)
n = S.pilot.sample_per_axis;
nc = size(L.index, 2);
h = sa_lattice('steps', L);
in = C.in_range(:);
if ~any(in)
    nodes = central_node(L); return
end
vals = cell(1, nc);
for c = 1:nc
    lo = min(L.coordinate(in, c)); hi = max(L.coordinate(in, c));
    v = linspace(lo, hi, n);
    k = round((v - L.anchor(c)) / h(c));
    k = min(max(k, L.index_min(c)), L.index_max(c));
    vals{c} = unique(k);
end
g = cell(1, nc);
[g{:}] = ndgrid(vals{:});
K = cell2mat(cellfun(@(x) x(:), g, 'UniformOutput', false));
nodes = zeros(0, 1);
for r = 1:size(K, 1)
    hit = find(all(L.index == K(r, :), 2), 1);
    if ~isempty(hit), nodes(end+1, 1) = hit; end %#ok<AGROW>
end
nodes = unique([nodes; central_node(L)]);
end

function tf = cross(L)
q = central_node(L);
kc = L.index(q, :);
tf = sum(L.index ~= kc, 2) <= 1;
end

function h = fingerprint(guess, coordinate, solver)
h = '';
if ~(isstruct(guess) && isfield(guess, 'x') && isfield(guess, 'y') && isnumeric(guess.x) && ...
        isnumeric(guess.y))
    return
end
md = java.security.MessageDigest.getInstance('SHA-256');
put = @(v) md.update(typecast(typecast(double(v(:)), 'uint8'), 'int8'));
put(size(guess.x)); put(guess.x);
put(size(guess.y)); put(guess.y);
put(coordinate);
put([solver.reltol solver.abstol solver.nmax solver.mesh_reject]);
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
end
