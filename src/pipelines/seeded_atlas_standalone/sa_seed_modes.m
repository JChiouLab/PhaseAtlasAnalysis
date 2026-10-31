function varargout = sa_seed_modes(action, varargin)
%SA_SEED_MODES  The seed_modes setting: active seed types and the pre-run count
%   (milestone 8, packet sections 4.3 and 5).
%
%   settings.seed_modes is a cell row of LEFT-drawn half shapes ('0', 'n+',
%   'n-', n >= 1) or 'all'. Each listed shape is seeded on both halves;
%   families are the ordered pairs of the listed shapes' type ids (filtered
%   by the census and by settings.families). The model module owns the
%   mapping shape -> type id, base file and transform, the type order and
%   the mirror rule (optional member seed_shape); this file names none of
%   them.
%
%   [M, S, A] = SA_SEED_MODES('activate', M, S, manifest)
%       After settings are resolved (zero solves). Resolves 'all' from the
%       library folder (every shape it offers; each skipped 'empty' mode
%       file is named in a printed line) and stores the resolved list in
%       S.seed_modes; maps the shapes to type ids in the module's type
%       order; when the types differ from the module's library_types,
%       M.library_types becomes the active list and M.mode_families the map
%       built from it. A = struct(shapes, types, skipped, default). At the
%       default seed_modes M and S are returned unchanged.
%   lib = SA_SEED_MODES('load', M, manifest, theta, parameters, A)
%       The module's library 'load' with the active types (the 0.11 call
%       when A.default).
%   L = SA_SEED_MODES('table', M, C, T, S)
%       The pre-run lines (zero solves): one per seedable family (and one
%       per pair with 0 seeds, '0 (no cell with both halves)'), with the
%       expected drawn label left / mirror(right), the library seeds, the
%       cold calls this run declares and rough hours at the module's
%       measured cold rate; families with a new shape say 'no measured
%       rate; legacy rate used'. Totals last.
%   SA_SEED_MODES('print', M, C, T, S)   prints the table

switch action
    case 'activate', [varargout{1}, varargout{2}, varargout{3}] = activate(varargin{:});
    case 'load',     varargout{1} = do_load(varargin{:});
    case 'table',    varargout{1} = table_lines(varargin{:});
    case 'print'
        L = table_lines(varargin{:});
        fprintf('%s\n', strjoin(L, newline));
    case 'expected', varargout{1} = expected(varargin{:});
    otherwise, error('SA:SEED:E_ACTION', 'Unknown sa_seed_modes action %s.', action);
end
end

% =========================================================================

function [M, S, A] = activate(M, S, manifest)
A = struct('shapes', {S.seed_modes}, 'types', {M.library_types}, ...
    'skipped', struct('file', {}, 'status', {}), 'default', true);
if ~isfield(M, 'seed_shape')
    if ~isequal(S.seed_modes, {'0', '1+', '2-'})
        error('SA:SEED:E_MODEL', 'Model %s offers no seed shapes beyond the default.', M.id);
    end
    return
end
shapes = S.seed_modes;
if ischar(shapes) || isequal(shapes, {'all'})
    [shapes, skipped] = M.seed_shape('all', manifest.directory);
    for k = 1:numel(skipped)
        fprintf('seed_modes all: skipped %s (run_status ''%s'')\n', skipped(k).file, skipped(k).status);
    end
    A.skipped = skipped;
    S.seed_modes = shapes;
end
types = cell(1, numel(shapes));
for k = 1:numel(shapes)
    i = M.seed_shape('resolve', shapes{k});
    types{k} = i.type;
end
types = M.seed_shape('order', types);
shapes = cellfun(@(t) M.seed_shape('shape_of_type', t), types, 'UniformOutput', false);
A.shapes = shapes;
A.types = types;
if isequal(types, M.library_types) && isempty(A.skipped) && ~ischar(S.seed_modes) && ...
        isequal(S.seed_modes, shapes)
    return                                   % the default: 0.11 exactly
end
if ischar(S.seed_modes) || ~isequal(sort(S.seed_modes), sort(shapes))
    S.seed_modes = shapes;
end
A.default = false;
M.library_types = types;
M.mode_families = M.seed_shape('mode_families', types);
% per-piece convenience warning: (guess_nodes - 1) / n < 50
for k = 1:numel(types)
    i = M.seed_shape('resolve', shapes{k});
    if i.n >= 3 && (S.model_options.guess_nodes - 1) / i.n < 50
        warning('SA:SEED:W_PIECE_NODES', ['seed shape %s: %g guess points per piece ' ...
            '((model_options.guess_nodes - 1) / %d < 50).'], shapes{k}, ...
            (S.model_options.guess_nodes - 1) / i.n, i.n);
    end
end
end

function lib = do_load(M, manifest, theta, parameters, A)
if A.default
    lib = M.library('load', manifest, theta, parameters);
else
    lib = M.library('load', manifest, theta, parameters, A.types);
end
end

function e = expected(M, tL, tR)
% the drawn label of a pair: left shape / mirror(right shape)
if isfield(M, 'seed_shape')
    sL = M.seed_shape('shape_of_type', tL);
    sR = M.seed_shape('mirror', M.seed_shape('shape_of_type', tR));
    e = sprintf('%s/%s', sL, sR);
else
    e = '';
end
end

function L = table_lines(M, C, T, S)
rate = NaN;
if isfield(M, 'seed_shape'), rate = M.seed_shape('cold_rate'); end
seeds = sa_targets('seed_counts', C);
legacy = {};
if isfield(M, 'seed_shape')
    legacy = cellfun(@(s) M.seed_shape('resolve', s), M.seed_shape('defaults'), 'UniformOutput', false);
    legacy = cellfun(@(i) i.type, legacy, 'UniformOutput', false);
end
L = {sprintf('pre-run seed count (zero solves; targets.policy %s; seed_modes %s)', ...
    S.targets.policy, strjoin(S.seed_modes, ' ')), ...
    sprintf('%-18s %-9s %7s %11s %8s', 'family', 'expected', 'seeds', 'cold calls', '~h')};
tot = [0 0 0];
note = false;
fam_sel = T.families;
for f = 1:numel(C.names)
    tL = C.family_types{f, 1}; tR = C.family_types{f, 2};
    ncold = sum(strcmp({T.list.family}, C.names{f}));
    e = expected(M, tL, tR);
    if C.counts(f) == 0
        L{end+1} = sprintf('%-18s %-9s %s', C.names{f}, e, '0 (no cell with both halves)'); %#ok<AGROW>
        continue
    end
    if ~any(strcmp(C.names{f}, fam_sel)) && ncold == 0
        continue                              % seedable, not selected (settings.families)
    end
    h = ncold * rate / 3600;
    new = ~(any(strcmp(tL, legacy)) && any(strcmp(tR, legacy)));
    tag = '';
    if new, tag = '  (no measured rate; legacy rate used)'; note = true; end
    L{end+1} = sprintf('%-18s %-9s %7d %11d %8.2f%s', C.names{f}, e, seeds(f), ncold, h, tag); %#ok<AGROW>
    tot = tot + [seeds(f) ncold h];
end
L{end+1} = sprintf('%-18s %-9s %7d %11d %8.2f', 'total', '', tot(1), tot(2), tot(3));
L{end+1} = sprintf(['hours are rough: %g s per cold call (the model''s measured cold rate)%s; ' ...
    'continuation, refinement and the per-root pass are not counted.'], rate, ...
    tern(note, '; families with a new shape have no measured rate', ''));
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
