function varargout = seeded_atlas_compare(varargin)
%SEEDED_ATLAS_COMPARE  Post-run comparison of a run folder with a declared
%   reference. Zero solves. A reference
%   never seeds or steers a run: this tool only reads a finished folder.
%
%   C = SEEDED_ATLAS_COMPARE(run_dir, reference)
%   C = SEEDED_ATLAS_COMPARE(run_dir, reference, opts)
%       reference: a reference struct or the .mat file declaring one
%       (variable `reference`, schema 'seeded_atlas/reference-1'):
%         .schema, .name, .source (text)
%         .nodes  struct array: coordinate (row, lattice coordinate units),
%                 usable (logical), x (row, the reference's own mesh),
%                 profile (one row per identity field of the model, in the
%                 module's identity_fields order, numel(x) columns), label
%       Per lattice node of the run: our roots (seeded_atlas_roots, verified
%       by default; opts.verify) against the usable reference entries at the
%       same coordinate (|difference| <= 1e-9 x step on every axis): the
%       distance of a pair is max over identity fields of
%       max_x |ours(x) - ref(x)| / max_x |ref(x)| on the reference mesh
%       (ours by deval; the myxo comparison convention); a node's class from the
%       smallest pair distance: agree (<= opts.agree, default 1e-6) | band |
%       different (>= opts.different, default 1e-2); ours_only |
%       reference_only. Report classes only: no tolerance is commissioned.
%       Memory: the roots are loaded without their iterates; every
%       representative attempt's sol is first read once, in root order, and
%       released (the file check a load of every sol makes, in its order),
%       then each root at a node with a usable reference entry reads its sol
%       again for its distances; at most one sol is held at a time.
%       run_dir: a run folder of either layout (sa_layout); the state folder
%       is read. Writes <run_dir>/compare/<name>.csv and <name>_map.txt (the reach
%       map with the reference overlaid: A agree, b band, D different,
%       o ours-only, r reference-only, . neither).
%   ref = SEEDED_ATLAS_COMPARE('from_myxo_sheet', file, step, name)
%       Converter (ported from tests/seeded_atlas_standalone/
%       sag_myxo_compare.m): a myxo coupled sheet (variables ok_grid, prof
%       (25 x 25 x 8 x nx), xg); cell (i, j) <-> lattice index [j - 13,
%       i - 13] at the given step (anchor 0); usable = ok_grid and every
%       value finite and > 0 (negative cells solve the clamped problem).
%   ref = SEEDED_ATLAS_COMPARE('from_stageb_registry', registry_csv, raw_dir, M, family, name)
%       Converter for a Stage-B accepted-object registry (rows family, i, j,
%       raw_file; '#' comment lines): each row of `family` becomes a usable
%       entry at [raw_dP raw_dG] with the RAW file's raw_sol on its own mesh
%       (identity rows of the module M). Read only.

if ischar(varargin{1}) && any(strcmp(varargin{1}, {'from_myxo_sheet', 'from_stageb_registry'}))
    switch varargin{1}
        case 'from_myxo_sheet',      varargout{1} = from_myxo(varargin{2:end});
        case 'from_stageb_registry', varargout{1} = from_stageb(varargin{2:end});
    end
    return
end
varargout{1} = compare(varargin{:});
end

% =========================================================================

function C = compare(run_dir, ref, opts)
if nargin < 3, opts = struct(); end
if ischar(ref)
    Z = load(ref, 'reference'); ref = Z.reference;
end
if ~isfield(ref, 'schema') || ~strcmp(ref.schema, 'seeded_atlas/reference-1')
    error('SA:COMPARE:E_REFERENCE', 'reference must declare schema seeded_atlas/reference-1.');
end
ag = getf(opts, 'agree', 1e-6); df = getf(opts, 'different', 1e-2);
E = load(fullfile(sa_layout('state', run_dir), 'settings_effective.mat'), 'settings', 'lattice');
E.settings = sa_settings('legacy_fill', E.settings);   % older folders load with neutral leaves
L = E.lattice;
M = sa_model(E.settings.model);
R = seeded_atlas_roots(run_dir, struct('verify', getf(opts, 'verify', true)));
D = sa_layout('state', run_dir);
for k = 1:numel(R)                               % every sol file, in root order, one at a time
    rep_sol(D, R(k));
end
tol = 1e-9 * sa_lattice('steps', L);
n = numel(L.key);
cls = repmat({''}, n, 1); dmin = nan(n, 1); nref = zeros(n, 1); nours = zeros(n, 1);
refnode = zeros(numel(ref.nodes), 1);
for k = 1:numel(ref.nodes)
    c = ref.nodes(k).coordinate;
    hit = find(all(abs(L.coordinate - c) <= tol, 2), 1);
    if ~isempty(hit), refnode(k) = hit; end
end
rk = {R.key};
for q = 1:n
    mine = R(strcmp(rk, L.key{q}));
    refs = ref.nodes(refnode == q & [ref.nodes.usable].');
    nours(q) = numel(mine); nref(q) = numel(refs);
    if isempty(mine) && isempty(refs), continue; end
    if isempty(refs), cls{q} = 'ours_only'; continue; end
    if isempty(mine), cls{q} = 'reference_only'; continue; end
    d = Inf;
    for a = 1:numel(mine)
        sol = rep_sol(D, mine(a));
        for b = 1:numel(refs)
            d = min(d, dist(sol, refs(b), [M.identity_fields.row]));
        end
        sol = []; %#ok<NASGU>
    end
    dmin(q) = d;
    if d <= ag, cls{q} = 'agree'; elseif d >= df, cls{q} = 'different'; else, cls{q} = 'band'; end
end
C = struct('name', ref.name, 'key', {L.key}, 'class', {cls}, 'distance', dmin, ...
    'n_ours', nours, 'n_reference', nref, 'unmatched_reference', sum(refnode == 0));
for c = {'agree', 'band', 'different', 'ours_only', 'reference_only'}
    C.counts.(c{1}) = sum(strcmp(cls, c{1}));
end
od = fullfile(run_dir, 'compare');
if exist(od, 'dir') ~= 7, mkdir(od); end
lines = {'key,index,coordinate,class,min_distance,n_ours,n_reference'};
for q = 1:n
    if isempty(cls{q}), continue; end
    lines{end+1} = sprintf('"%s","%s","%s",%s,%.6g,%d,%d', L.key{q}, mat2str(L.index(q, :)), ...
        mat2str(L.coordinate(q, :), 10), cls{q}, dmin(q), nours(q), nref(q)); %#ok<AGROW>
end
sa_write('text', fullfile(od, [ref.name '.csv']), [strjoin(lines, newline) newline]);
M2 = {sprintf('reach map vs %s (A agree, b band, D different, o ours-only, r reference-only, . neither)', ref.name)};
if size(L.index, 2) == 2
    k1 = unique(L.index(:, 1)); k2 = unique(L.index(:, 2));
    for a = k1(end:-1:1).'
        s = blanks(numel(k2));
        for b = 1:numel(k2)
            s(b) = sym_of(cls{find(L.index(:, 1) == a & L.index(:, 2) == k2(b), 1)});
        end
        M2{end+1} = sprintf('%4d %s', a, s); %#ok<AGROW>
    end
end
M2{end+1} = sprintf('counts: agree %d, band %d, different %d, ours-only %d, reference-only %d; reference entries off the lattice %d', ...
    C.counts.agree, C.counts.band, C.counts.different, C.counts.ours_only, C.counts.reference_only, ...
    C.unmatched_reference);
sa_write('text', fullfile(od, [ref.name '_map.txt']), [strjoin(M2, newline) newline]);
end

function sol = rep_sol(D, r)
% the root's representative iterate, as seeded_atlas_roots' sol option reads it
Z = load(fullfile(D, 'attempts', [r.representative_attempt_id '.mat']), 'sol');
sol = Z.sol;
end

function d = dist(sol, r, rows)
if r.x(1) < sol.x(1) - 1e-12 || r.x(end) > sol.x(end) + 1e-12, d = NaN; return; end
A = deval(sol, r.x, rows);
B = r.profile;
d = 0;
for k = 1:size(B, 1)
    s = max(abs(B(k, :)));
    if ~(s > 0) || any(~isfinite(A(k, :))), d = NaN; return; end
    d = max(d, max(abs(A(k, :) - B(k, :))) / s);
end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

% =========================================================================
% converters

function ref = from_myxo(file, step, name)
Z = load(file, 'ok_grid', 'prof', 'xg');
nodes = struct('coordinate', {}, 'usable', {}, 'x', {}, 'profile', {}, 'label', {});
[ni, nj] = size(Z.ok_grid);
c0 = (ni + 1) / 2;
for i = 1:ni
    for j = 1:nj
        P = reshape(Z.prof(i, j, :, :), size(Z.prof, 3), size(Z.prof, 4));
        ok = logical(Z.ok_grid(i, j)) && all(isfinite(P(:))) && min(P(:)) > 0;
        nodes(end+1) = struct('coordinate', step * [j - c0, i - c0], 'usable', ok, ...
            'x', reshape(Z.xg, 1, []), 'profile', P, 'label', sprintf('myxo(%d,%d)', i, j)); %#ok<AGROW>
    end
end
ref = struct('schema', 'seeded_atlas/reference-1', 'name', name, 'source', file, 'nodes', nodes);
end

function ref = from_stageb(csvfile, raw_dir, M, family, name)
txt = strsplit(fileread(csvfile), newline);
txt = txt(~startsWith(txt, '#') & ~cellfun(@isempty, strtrim(txt)));
h = split_csv(txt{1});
fi = find(strcmp(h, 'family')); rf = find(strcmp(h, 'raw_file'));
nodes = struct('coordinate', {}, 'usable', {}, 'x', {}, 'profile', {}, 'label', {});
rows = [M.identity_fields.row];
for k = 2:numel(txt)
    r = split_csv(txt{k});
    if ~strcmp(r{fi}, family), continue; end
    B = load(fullfile(raw_dir, r{rf}), 'raw_sol', 'raw_dP', 'raw_dG');
    nodes(end+1) = struct('coordinate', [B.raw_dP B.raw_dG], 'usable', true, ...
        'x', B.raw_sol.x, 'profile', B.raw_sol.y(rows, :), 'label', r{rf}); %#ok<AGROW>
end
ref = struct('schema', 'seeded_atlas/reference-1', 'name', name, 'source', csvfile, 'nodes', nodes);
end

function c = split_csv(s)
c = strrep(regexp(strtrim(s), ',(?=(?:[^"]*"[^"]*")*[^"]*$)', 'split'), '"', '');
end

function c = sym_of(cl)
switch cl
    case 'agree', c = 'A';
    case 'band', c = 'b';
    case 'different', c = 'D';
    case 'ours_only', c = 'o';
    case 'reference_only', c = 'r';
    otherwise, c = '.';
end
end
