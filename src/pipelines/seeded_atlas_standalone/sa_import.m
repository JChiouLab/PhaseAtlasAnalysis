function varargout = sa_import(action, varargin)
%SA_IMPORT  Roots of an earlier run as seeds ("import"; milestone 10).
%
%   A new run reads a finished run folder (the SOURCE, settings
%   import.source_dir). Each source root goes to the nearest node of the new
%   lattice and is solved there ONCE from its stored iterate (route, kind
%   and phase 'import'; S.solver, no cold cap, no halving, no retry); an
%   accepted import is a starting root like a cold root. Optionally the
%   mirror image of each source root is imported too (module member
%   reflection; sa_model). Zero solves here: the calls are made by
%   sa_search's import stage through the one executor.
%
%   P = SA_IMPORT('blank')
%       The inactive plan (no import): active false, no rows.
%   P = SA_IMPORT('plan', S, M, theta, parameters, L, fs)
%       The plan, built once per new run before its folder is prepared and
%       stored as effective.import (never rebuilt on resume):
%       compatibility (SA:IMPORT:E_SOURCE, E_LIVE, E_MODEL, E_PARAMETERS,
%       E_STATES, E_OPTION, E_EMPTY), the source list
%       (seeded_atlas_roots(source, struct('verify', false))), the nearest
%       node of each root (|| (c - c_node) ./ step ||_2, ties to the lowest
%       lattice row; no_node_near beyond import.max_map_steps; outside_box
%       outside continuation.box / region), mirror rows with their twin
%       dedup, the source stamp, the SHA-256 of every planned row's source
%       attempt file, and (import.source_edges not 'none') the source's
%       edges. fs: the new run's field scale values (identity rows).
%       Row fields: i, kind (direct | mirror), source_root_id,
%       source_attempt_id, source_attempt_sha256, source_coordinate,
%       source_label, family, node, key, index, coordinate, map_distance,
%       exact, twin_root_id, twin_distance, status (planned | no_node_near |
%       outside_box | mirror_twin_planned).
%       Selection (import.root_ids, a cell row or a CSV file with a root_id
%       column; import.families, observed family labels): only the source
%       roots that match either list get rows (direct and mirror); both
%       empty = every root, and the plan is then the one built without
%       these settings, field for field. An id or a label the source does
%       not hold is refused (SA:IMPORT:E_SELECTION, listing them). Fields
%       added only when used: selection (the lists, selected / skipped
%       counts, the selected labels), fence (continuation.families
%       'selected': the selected labels plus their mirror labels when the
%       mirror is on; refused SA:IMPORT:E_FENCE when that list is empty,
%       so the fence is never switched off silently), source_edges_auto (import.source_edges 'auto':
%       'same_or_stronger' when the edge conditions hold, else 'none' with
%       a printed note), source_stopped (a killed source accepted with
%       import.source_stopped: the call it was killed at, the calls
%       without an end line that are ignored, the checkpoint used; its
%       identity is the SHA-256 of each planned row's attempt file, not the
%       folder stamp, so the source may be resumed later; its source-edge
%       knowledge comes from its last checkpoint queue.mat only; a journal
%       written in the last 10 minutes is refused).
%       Source-edge reuse (same_or_stronger, all) needs the lattice, the
%       solver and validity groups, continuation guess_mesh / guess_nodes
%       / step_rule and the model version (corrector and solver code) of
%       the source.
%   tf = SA_IMPORT('active', ctx_or_plan)
%   [st, parts] = SA_IMPORT('stamp', source_dir)
%       SHA-256 over the source's internal/ settings_effective.mat,
%       attempts.mat, queue_final.mat, journal.txt, compaction.mat (or ''),
%       and the attempt-file count; plots/, traces/ and the human files are
%       not part of it.
%   st = SA_IMPORT('current_stamp', S [, P])   the stamp of
%       S.import.source_dir ('' when no import is set or the folder holds no
%       readable run); with the stored plan P of a killed source
%       (P.source_stopped) the SHA-256 over its planned attempt files
%       instead ('' when one is missing)
%   st = SA_IMPORT('name_stamp', imp)   the stamp a default run-folder name
%       shows (sa_run_name): the folder stamp, or with
%       imp.source_stopped true the SHA-256 of the source's
%       settings_effective.mat (constant while the source is resumed)
%   f = SA_IMPORT('fence', S, P)   the family fence in force: {} (none),
%       S.continuation.families, or for 'selected' the plan's fence
%   ids = SA_IMPORT('read_root_ids', csv_file)   the root_id column
%   [u, tg] = SA_IMPORT('unit', ctx, row)   the work unit and target row
%   g = SA_IMPORT('guess', ctx, row)
%       The guess of one import call, live and on replay: the source attempt
%       file, its SHA-256 checked against the plan (SA:INTEGRITY), then
%       'source': struct('x', sol.x, 'y', sol.y) (mirror rows: y(p, :)) or
%       'rule': sa_guess_mesh('from', sol, S) of the (mirrored) iterate.
%   K = SA_IMPORT('knowledge', ctx)
%       Source-edge knowledge (import.source_edges 'same_or_stronger' |
%       'all'), after the import stage: source roots map to the new roots
%       their accepted DIRECT imports landed on (merged roots pool their
%       knowledge; refused imports carry none). Sets (keys '<new root>|<row>'):
%       pass_a (a source edge r_s -> q with calls > 0: refused with h_s >=
%       the pass's h under 'same_or_stronger', any h under 'all'; completed
%       only when its child mapped to an accepted import), pass_b (child|p
%       of a completed source edge whose child mapped), lin_ss (refused
%       r_s -> q with h_s >= continuation.max_halvings), lin_weak (refused
%       with h_s below it; used under 'all' only). Edges that ended on a
%       cap or unfinished carry nothing, nor do secant rescue edges (their
%       guess is not the parent's iterate); rescue uses none.
%   L = SA_IMPORT('report_lines', ctx)   the summary.txt Import section
%   f = SA_IMPORT('family_of', label, seed_family, types, map)
%       The shared family rule (sa_root_pass's lookup): the library family
%       of the mode name, else the creating discovery's seed family.
%   m = SA_IMPORT('mirror_label', label)  halves swapped, odd modes' signs
%       flipped (1+/0 <-> 0/1-, 2-/1- <-> 1+/2-); 'A_B' -> 'B_A'
%   f = SA_IMPORT('mirror_family', family)  'A_B' -> 'B_A'
%   n = SA_IMPORT('counter' [, 'reset'])   plan builds in this MATLAB
%       session (T106: a resume does not list the source again)
%   d = SA_IMPORT('distance', Pa, Pb, fs)  scaled_rms_v1 of two comparison
%       profiles (identity rows on one comparison mesh)

switch action
    case 'blank',          varargout{1} = blank_plan();
    case 'plan',           varargout{1} = plan(varargin{:});
    case 'active',         varargout{1} = is_active(varargin{:});
    case 'stamp',          [varargout{1}, varargout{2}] = stamp(varargin{:});
    case 'current_stamp',  varargout{1} = current_stamp(varargin{:});
    case 'name_stamp',     varargout{1} = name_stamp(varargin{:});
    case 'fence',          varargout{1} = fence_of(varargin{:});
    case 'read_root_ids',  varargout{1} = read_root_ids(varargin{:});
    case 'unit',           [varargout{1}, varargout{2}] = unit_of(varargin{:});
    case 'guess',          varargout{1} = guess_of(varargin{:});
    case 'knowledge',      varargout{1} = knowledge(varargin{:});
    case 'report_lines',   varargout{1} = report_lines(varargin{:});
    case 'family_of',      varargout{1} = family_of(varargin{:});
    case 'mirror_label',   varargout{1} = mirror_label(varargin{:});
    case 'mirror_family',  varargout{1} = mirror_family(varargin{:});
    case 'counter',        varargout{1} = counter(varargin{:});
    case 'distance',       varargout{1} = prof_distance(varargin{:});
    case 'box_mask',       varargout{1} = box_mask(varargin{:});
    otherwise, error('SA:IMPORT', 'Unknown sa_import action %s.', action);
end
end

% =========================================================================
% plan

function P = blank_plan()
P = struct('active', false, 'source_dir', '', 'stamp', '', 'stamp_parts', struct(), ...
    'source_software', '', 'source_model_version', '', 'differing_leaves', {{}}, ...
    'guess_mesh', '', 'mirror', false, 'source_edges_mode', 'none', ...
    'rows', blank_rows(), 'edges', blank_edges(), ...
    'source_n_roots', 0, 'source_node_counts', zeros(1, 0), 'source_max_roots_per_node', NaN, ...
    'source_nodes_at_cap', 0, 'plan_s', 0, 'hash_s', 0, 'hash_bytes', 0);
end

function R = blank_rows()
R = struct('i', {}, 'kind', {}, 'source_root_id', {}, 'source_attempt_id', {}, ...
    'source_attempt_sha256', {}, 'source_coordinate', {}, 'source_label', {}, 'family', {}, ...
    'node', {}, 'key', {}, 'index', {}, 'coordinate', {}, 'map_distance', {}, 'exact', {}, ...
    'twin_root_id', {}, 'twin_distance', {}, 'status', {});
end

function E = blank_edges()
E = struct('parent_root', {}, 'p_key', {}, 'q_key', {}, 'status', {}, 'calls', {}, ...
    'child_root', {}, 'kind', {}, 'rule', {}, 'h', {});
end

function tf = is_active(x)
tf = false;
if isstruct(x) && isfield(x, 'import'), x = x.import; end
if isstruct(x) && isfield(x, 'active'), tf = logical(x.active); end
end

function n = counter(varargin)
persistent c
if isempty(c), c = 0; end
if nargin >= 1 && strcmp(varargin{1}, 'reset'), c = 0; end
if nargin >= 1 && strcmp(varargin{1}, 'add'), c = c + 1; end
n = c;
end

function P = plan(S, M, theta, parameters, L, fs) %#ok<INUSL>
t0 = tic;
counter('add');
P = blank_plan();
src = S.import.source_dir;
D = source_state(src, S);
E0 = load(fullfile(D, 'settings_effective.mat'));
KS = check_live(D, leaf(S.import, 'source_stopped', false));   % killed and accepted: KS.killed
% ---- compatibility ---------------------------------------------------------
if ~isfield(E0, 'settings') || ~strcmp(E0.settings.model, M.id)
    error('SA:IMPORT:E_MODEL', 'import refused: the source run is model %s, this run is %s.', ...
        E0.settings.model, M.id);
end
if ~isfield(E0, 'parameters') || ~isequaln(E0.parameters, parameters)
    error('SA:IMPORT:E_PARAMETERS', 'import refused: the parameters differ from the source run''s (%s).', ...
        strjoin(param_diff(E0.parameters, parameters), ', '));
end
if E0.model.n_states ~= M.n_states || ~isequal({E0.model.fields.name}, {M.fields.name})
    error('SA:IMPORT:E_STATES', 'import refused: the state count or field names differ from the source run''s.');
end
Es = sa_settings('legacy_fill', E0.settings);
mode = S.import.source_edges;
AUTO = [];
if ~strcmp(mode, 'none')
    why = edge_conditions(Es, E0.lattice, S, L, E0.model.version, M.version);
    if strcmp(mode, 'auto')
        % the focus default: reuse when the conditions hold, else none
        if isempty(why), mode = 'same_or_stronger'; else, mode = 'none'; end
        AUTO = struct('resolved', mode, 'why', {why});
        if ~isempty(why)
            fprintf('import.source_edges auto -> none (source-edge knowledge not reused: %s)\n', strjoin(why, '; '));
        end
    elseif ~isempty(why)
        error('SA:IMPORT:E_OPTION', ['import.source_edges = ''%s'' needs the source''s lattice, solver and ' ...
            'validity groups, continuation guess mesh / guess nodes / step rule and model version ' ...
            'unchanged: %s.'], mode, strjoin(why, '; '));
    end
end
mirror = S.import.mirror;
if mirror && ~isfield(M, 'reflection')
    error('SA:IMPORT:E_OPTION', 'import.mirror = true but model %s declares no reflection.', M.id);
end
P.active = true;
P.source_dir = src;
if ~KS.killed
    [P.stamp, P.stamp_parts] = stamp(src);
end
P.source_software = Es.software_version;
P.source_model_version = E0.model.version;
P.differing_leaves = sa_resume('scientific_diff', Es, S);
P.guess_mesh = S.import.guess_mesh;
P.mirror = mirror;
P.source_edges_mode = mode;
% ---- the source roots, attempts ---------------------------------------------
R = seeded_atlas_roots(src, struct('verify', false));
Z = load(fullfile(D, 'attempts.mat'), 'attempts');
A = Z.attempts;
ids = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:numel(A), ids(A(k).attempt_id) = k; end
if KS.killed
    % a killed source: roots whose representative call has no end line are
    % ignored; calls after the last checkpoint are read from their own
    % attempt files (append-only)
    keep = isKey(KS.ended, {R.representative_attempt_id});
    KS.ignored_roots = {R(~keep).root_id};
    R = R(keep);
    for r = 1:numel(R)
        id = R(r).representative_attempt_id;
        if ~isKey(ids, id)
            Za = load(fullfile(D, 'attempts', [id '.mat']), 'attempt');
            A(end+1) = Za.attempt; %#ok<AGROW>
            ids(id) = numel(A);
        end
    end
end
types = M.library_types;
map = [];
if isfield(M, 'mode_families'), map = M.mode_families; end
nR = numel(R);
P.source_n_roots = nR;
lab_all = cell(1, nR);
for r = 1:nR, lab_all{r} = A(ids(R(r).representative_attempt_id)).observed_family; end
[sel, SEL] = selection(S, R, lab_all, KS);
[uk, ~, j] = unique({R.key});
cnt = accumarray(j(:), 1).';
P.source_node_counts = cnt;
P.source_max_roots_per_node = Es.multiseed.max_roots_per_node;
P.source_nodes_at_cap = sum(cnt >= Es.multiseed.max_roots_per_node);
% ---- map to the new lattice ---------------------------------------------------
h = sa_lattice('steps', L);
box = box_mask(L, S);
fsv = fs(:);
rows = blank_rows();
for r = find(sel)
    a = A(ids(R(r).representative_attempt_id));
    lab = a.observed_family;
    row = blank_row();
    row.kind = 'direct';
    row.source_root_id = R(r).root_id;
    row.source_attempt_id = R(r).representative_attempt_id;
    row.source_coordinate = R(r).coordinate;
    row.source_label = lab;
    row.family = family_of(lab, a.family, types, map);
    row = map_row(row, R(r).coordinate, L, h, box, S.import.max_map_steps);
    rows(end+1) = row; %#ok<AGROW>
end
direct = rows;
% ---- mirror rows (twin dedup at plan time, zero solves) ------------------------
mrows = blank_rows();
if mirror
    pi_ = M.reflection.identity_perm;
    hs = sa_lattice('steps', E0.lattice);
    byk = containers.Map(uk, num2cell(1:numel(uk)));
    members = accumarray(j(:), (1:nR).', [], @(v) {sort(v).'});
    planned_src = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    for r = 1:numel(direct)
        planned_src(direct(r).source_root_id) = strcmp(direct(r).status, 'planned');
    end
    for r = find(sel)
        a = A(ids(R(r).representative_attempt_id));
        mc = -R(r).coordinate;
        row = blank_row();
        row.kind = 'mirror';
        row.source_root_id = R(r).root_id;
        row.source_attempt_id = R(r).representative_attempt_id;
        row.source_coordinate = R(r).coordinate;
        row.source_label = mirror_label(a.observed_family);
        row.family = family_of(row.source_label, mirror_family(a.family), types, map);
        row = map_row(row, mc, L, h, box, S.import.max_map_steps);
        % the twin: a root at the SOURCE node nearest -c whose comparison
        % profile is within similarity_tolerance of the permuted profile
        dsrc = sqrt(sum(((mc - E0.lattice.coordinate) ./ hs) .^ 2, 2));
        [~, qs] = min(dsrc);
        tk = E0.lattice.key{qs};
        if isKey(byk, tk)
            Pm = R(r).representative_profile(pi_, :);
            best = Inf; bid = '';
            for t = members{byk(tk)}
                d = prof_distance(R(t).representative_profile, Pm, fsv);
                if d <= S.identity.similarity_tolerance && d < best
                    best = d; bid = R(t).root_id;
                end
            end
            if ~isempty(bid)
                row.twin_root_id = bid; row.twin_distance = best;
                if strcmp(row.status, 'planned') && isKey(planned_src, bid) && planned_src(bid)
                    row.status = 'mirror_twin_planned';
                end
            end
        end
        mrows(end+1) = row; %#ok<AGROW>
    end
end
% ---- order: direct by new lattice row then source order; then mirrors -----------
rows = [order_rows(direct), order_rows(mrows)];
for k = 1:numel(rows), rows(k).i = k; end
if ~S.routes.cold && ~any(strcmp({rows.status}, 'planned'))
    error('SA:IMPORT:E_EMPTY', ['import only (routes.cold false) with nothing planned: every source root ' ...
        'is skipped (no node near, outside the box, or a planned twin).']);
end
% ---- per-file SHA-256 of the planned rows' source attempt files -----------------
t1 = tic;
cache = containers.Map('KeyType', 'char', 'ValueType', 'char');
nb = 0;
for k = 1:numel(rows)
    if ~strcmp(rows(k).status, 'planned'), continue; end
    id = rows(k).source_attempt_id;
    if ~isKey(cache, id)
        f = fullfile(D, 'attempts', [id '.mat']);
        if exist(f, 'file') ~= 2
            error('SA:IMPORT:E_SOURCE', 'import refused: the source attempt file %s is missing.', f);
        end
        cache(id) = sa_sha256(f);
        info = dir(f); nb = nb + info.bytes;
    end
    rows(k).source_attempt_sha256 = cache(id);
end
P.hash_s = toc(t1);
P.hash_bytes = nb;
P.rows = rows;
% ---- source edges ----------------------------------------------------------------
if ~strcmp(mode, 'none')
    if KS.killed
        P.edges = source_edges(D, E0.lattice, 'queue.mat');      % the last checkpoint only
    else
        P.edges = source_edges(D, E0.lattice, 'queue_final.mat');
    end
end
% ---- fields present only when their setting is used --------------------------------
if ~isempty(SEL), P.selection = SEL; end
if ischar(S.continuation.families)
    % 'selected': the observed labels of the selected roots, plus their
    % mirror labels when mirror images are imported
    fl = lab_all(sel);
    if mirror, fl = [fl, cellfun(@mirror_label, fl, 'UniformOutput', false)]; end
    fl = unique(fl(~cellfun(@isempty, fl)));
    if isempty(fl)
        % never silently unfenced: an empty fence would switch the fence off
        error('SA:IMPORT:E_FENCE', ['import refused: continuation.families ''selected'' resolves to no label (the ' ...
            '%d selected source root(s) carry no observed family label): set continuation.families to a list ' ...
            'of labels, or to {} to run without a fence.'], sum(sel));
    end
    P.fence = reshape(fl, 1, []);
end
if ~isempty(AUTO), P.source_edges_auto = AUTO; end
if KS.killed
    [P.stamp, n] = rows_stamp(rows);
    P.stamp_parts = struct('kind', 'planned attempt files (killed source)', 'attempt_files', n);
    P.source_stopped = rmfield(KS, 'ended');
end
P.plan_s = toc(t0);
end

function v = leaf(s, name, d)
v = d;
if isstruct(s) && isfield(s, name), v = s.(name); end
end

function [sel, SEL] = selection(S, R, labels, KS)
% The source roots to import: those whose id is in import.root_ids or
% whose observed label is in import.families; both empty = all (SEL = []).
nR = numel(R);
sel = true(1, nR);
SEL = [];
ids_in = leaf(S.import, 'root_ids', {});
fams = leaf(S.import, 'families', {});
csv = '';
if ischar(ids_in)
    csv = ids_in;
    ids_in = read_root_ids(csv);
    if isempty(ids_in)
        error('SA:IMPORT:E_SELECTION', 'import refused: the root id file %s lists no root id.', csv);
    end
end
if isempty(ids_in) && isempty(fams), return; end
rid = {R.root_id};
eid = 'SA:IMPORT:E_SELECTION';
unk = setdiff(ids_in, rid, 'stable');
if ~isempty(unk)
    extra = '';
    if KS.killed && any(ismember(unk, KS.ignored_roots))
        extra = sprintf(' (%d of them belong to calls without an end line, which are ignored)', ...
            sum(ismember(unk, KS.ignored_roots)));
    end
    error(eid, 'import refused: import.root_ids names %d root id(s) the source does not hold%s: %s.', ...
        numel(unk), extra, list_text(unk));
end
unf = setdiff(fams, labels, 'stable');
if ~isempty(unf)
    error(eid, ['import refused: import.families names %d label(s) no source root carries: %s (the ' ...
        'source''s labels: %s).'], numel(unf), list_text(unf), list_text(unique(labels)));
end
sel = ismember(rid, ids_in) | ismember(labels, fams);
SEL = struct('root_ids', {reshape(ids_in, 1, [])}, 'root_ids_file', csv, 'families', {fams}, ...
    'n_source', nR, 'n_selected', sum(sel), 'n_skipped', nR - sum(sel), ...
    'labels', {reshape(unique(labels(sel)), 1, [])});
end

function t = list_text(c)
c = reshape(c, 1, []);
if numel(c) > 20, t = [strjoin(c(1:20), ', ') sprintf(', ... (%d more)', numel(c) - 20)];
else, t = strjoin(c, ', ');
end
end

function ids = read_root_ids(f)
% the root_id column of a CSV file (a header line; quoted or bare fields;
% a quoted field may hold commas), unique, in file order
if exist(f, 'file') ~= 2
    error('SA:IMPORT:E_SELECTION', 'import refused: the root id file %s does not exist.', f);
end
L = regexp(fileread(f), '\r?\n', 'split');
L = L(~cellfun(@(x) isempty(strtrim(x)), L));
if isempty(L)
    error('SA:IMPORT:E_SELECTION', 'import refused: the root id file %s is empty.', f);
end
h = csv_fields(L{1});
c = find(strcmp(strtrim(h), 'root_id'), 1);
if isempty(c)
    error('SA:IMPORT:E_SELECTION', 'import refused: the root id file %s has no root_id column (header: %s).', ...
        f, L{1});
end
ids = {};
for k = 2:numel(L)
    v = csv_fields(L{k});
    if numel(v) >= c && ~isempty(strtrim(v{c})), ids{end+1} = strtrim(v{c}); end %#ok<AGROW>
end
ids = unique(ids, 'stable');
end

function v = csv_fields(line)
% split one CSV line at commas outside double quotes; quotes removed
v = {};
cur = '';
q = false;
k = 1;
while k <= numel(line)
    ch = line(k);
    if ch == '"'
        if q && k < numel(line) && line(k+1) == '"', cur(end+1) = '"'; k = k + 1; %#ok<AGROW>
        else, q = ~q;
        end
    elseif ch == ',' && ~q
        v{end+1} = cur; cur = ''; %#ok<AGROW>
    else
        cur(end+1) = ch; %#ok<AGROW>
    end
    k = k + 1;
end
v{end+1} = cur;
end

function [h, n] = rows_stamp(rows)
% SHA-256 over the planned rows' (source attempt id, attempt-file SHA-256),
% each file once, in id order: a killed source's identity (its attempt files
% are append-only, so resuming the source later leaves it unchanged)
h = ''; n = 0;
if isempty(rows), pl = rows; else, pl = rows(strcmp({rows.status}, 'planned')); end
if isempty(pl)
    ids = {}; shas = {};
else
    [ids, k] = unique({pl.source_attempt_id});
    shas = {pl(k).source_attempt_sha256};
end
if any(cellfun(@isempty, shas)), return; end
n = numel(ids);
txt = ['planned_attempt_files;' strjoin(cellfun(@(a, b) [a '=' b ';'], ids, shas, 'UniformOutput', false), '')];
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(txt), 'int8'));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
end

function row = blank_row()
row = struct('i', 0, 'kind', '', 'source_root_id', '', 'source_attempt_id', '', ...
    'source_attempt_sha256', '', 'source_coordinate', [], 'source_label', '', 'family', '', ...
    'node', 0, 'key', '', 'index', [], 'coordinate', [], 'map_distance', NaN, 'exact', false, ...
    'twin_root_id', '', 'twin_distance', NaN, 'status', '');
end

function row = map_row(row, c, L, h, box, max_steps)
d = sqrt(sum(((c - L.coordinate) ./ h) .^ 2, 2));
[dm, q] = min(d);                       % the first minimum: ties to the lowest lattice row
row.map_distance = dm;
if dm > max_steps
    row.status = 'no_node_near';
    return
end
row.node = q; row.key = L.key{q}; row.index = L.index(q, :); row.coordinate = L.coordinate(q, :);
row.exact = isequal(L.coordinate(q, :), c);
if ~box(q)
    row.status = 'outside_box';
else
    row.status = 'planned';
end
end

function rows = order_rows(rows)
if isempty(rows), return; end
n = [rows.node]; n(n == 0) = Inf;
[~, o] = sortrows([n(:), (1:numel(rows)).']);   % new lattice row, then source order
rows = rows(o);
end

function m = box_mask(L, S)
% the continuation region (sa_continuation init's Q.box)
b = S.continuation.box;
if strcmp(S.continuation.region, 'cross')
    m = sa_modes('cross', L);
elseif isempty(b)
    m = true(L.n_nodes, 1);
else
    m = all(abs(L.index) <= b, 2);
end
m = logical(m(:));
end

function d = prof_distance(Pa, Pb, fs)
dd = (Pa - Pb) ./ fs(:);
d = sqrt(sum(dd(:) .^ 2) / numel(dd));
end

function D = source_state(src, S)
if ~(ischar(src) && ~isempty(src)) || exist(src, 'dir') ~= 7
    error('SA:IMPORT:E_SOURCE', 'import refused: import.source_dir %s is not a folder.', src);
end
% the top or its internal/ folder
top = src;
[up, nm] = fileparts(regexprep(src, '[\\/]+$', ''));
if strcmp(nm, 'internal') && exist(fullfile(src, 'settings_effective.mat'), 'file') == 2, top = up; end
if strcmp(sa_layout('layout', top), 'flat')
    error('SA:IMPORT:E_SOURCE', ['import refused: %s is a flat folder written by software 0.5 or earlier; ' ...
        'import reads the internal/ layout only.'], src);
end
D = sa_layout('state', top);
if exist(fullfile(D, 'settings_effective.mat'), 'file') ~= 2 || exist(fullfile(D, 'attempts.mat'), 'file') ~= 2
    error('SA:IMPORT:E_SOURCE', 'import refused: %s holds no finished run (settings_effective.mat, attempts.mat).', src);
end
if isfield(S, 'run_dir') && ~isempty(S.run_dir) && strcmp(canon(top), canon(S.run_dir))
    error('SA:IMPORT:E_SOURCE', 'import refused: the source is the run folder itself.');
end
end

function c = canon(d)
c = regexprep(char(java.io.File(d).getAbsolutePath()), '[\\/]+$', '');
try
    c = char(java.io.File(c).getCanonicalPath());
catch
end
end

function K = check_live(D, stopped)
% A running or killed source is refused; a budget-stopped one is not. With
% import.source_stopped (stopped) a killed source is accepted (K.killed),
% unless its journal was written in the last 10 minutes.
J = sa_journal('read', D);
st = {J.starts.id}; en = {J.ends.id};
open_ = setdiff(st, en);
f = dir(fullfile(D, 'attempts', 'A*.mat'));
n = sum(~cellfun(@isempty, regexp({f.name}, '^A\d{6}\.mat$', 'once')));
why = '';
qf = fullfile(D, 'queue_final.mat');
if ~isempty(open_)
    why = sprintf(['import refused: the source journal has a call started without an end (%s): ' ...
        'the source is running or was killed; resume it to a stop first.'], open_{1});
elseif exist(qf, 'file') ~= 2
    why = 'import refused: the source has no queue_final.mat (running or killed).';
else
    Q = load(qf, 'attempt_count');
    nq = n;
    if exist(fullfile(D, 'compaction.mat'), 'file') == 2
        Zc = load(fullfile(D, 'compaction.mat'), 'attempt_count_before');
        nq = Zc.attempt_count_before;
    end
    if Q.attempt_count ~= nq
        why = sprintf(['import refused: the source''s queue_final.mat counts %d calls but %d attempt ' ...
            'files exist (running or killed).'], Q.attempt_count, nq);
    end
end
K = struct('killed', false);
if isempty(why), return; end
if ~stopped
    error('SA:IMPORT:E_LIVE', ['%s If the source was killed and will not write again, set ' ...
        'import.source_stopped = true: its calls with an end line are imported (calls without one are ' ...
        'ignored); a source that may still be running must not be imported.'], why);
end
age = journal_age_s(D);
if age < 600
    error('SA:IMPORT:E_LIVE', ['import refused: import.source_stopped is set but the source journal was ' ...
        'written %.0f s ago (less than 10 minutes): the source may still be running; wait, or check that its ' ...
        'process has ended.'], age);
end
ended = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for k = 1:numel(en), ended(en{k}) = true; end
last = '';
if ~isempty(st), last = st{end}; end
qc = NaN; ac = NaN;
if exist(fullfile(D, 'queue.mat'), 'file') == 2
    Zq = load(fullfile(D, 'queue.mat'), 'attempt_count'); qc = Zq.attempt_count;
end
if exist(fullfile(D, 'attempts.mat'), 'file') == 2
    Za = load(fullfile(D, 'attempts.mat'), 'attempt_count'); ac = Za.attempt_count;
end
K = struct('killed', true, 'why', why, 'last_started', last, 'open_calls', {reshape(open_, 1, [])}, ...
    'attempt_files', n, 'checkpoint_queue', qc, 'checkpoint_attempts', ac, 'journal_age_s', age, ...
    'ignored_roots', {{}}, 'ended', ended);
end

function a = journal_age_s(D)
% seconds since the source journal was last written (file modification time)
f = fullfile(D, 'journal.txt');
a = Inf;
if exist(f, 'file') ~= 2, return; end
a = (double(java.lang.System.currentTimeMillis()) - double(java.io.File(f).lastModified())) / 1000;
end

function why = edge_conditions(Es, Lat0, S, L, v0, v1)
why = {};
if ~(isequal(Lat0.key, L.key) && isequal(Lat0.coordinate, L.coordinate))
    why{end+1} = 'the lattice differs (node keys or coordinates)';
end
if ~isequaln(Es.solver, S.solver), why{end+1} = 'the solver group differs'; end
if ~isequaln(Es.validity, S.validity), why{end+1} = 'the validity group differs'; end
for f = {'guess_mesh', 'guess_nodes', 'step_rule'}
    if ~isequaln(Es.continuation.(f{1}), S.continuation.(f{1}))
        why{end+1} = sprintf('continuation.%s differs', f{1}); %#ok<AGROW>
    end
end
if ~strcmp(v0, v1)
    % a refusal depends on the corrector and solver code (owner ruling)
    why{end+1} = 'the model version differs (corrector or solver code)';
end
end

function E = source_edges(D, Lat0, file)
E = blank_edges();
if exist(fullfile(D, file), 'file') ~= 2, return; end          % a killed source before its first checkpoint
Z = load(fullfile(D, file), 'Q');
if ~isstruct(Z.Q) || ~isfield(Z.Q, 'edges'), return; end
Q = Z.Q.edges;
for k = 1:numel(Q)
    e = Q(k);
    E(end+1) = struct('parent_root', e.parent_root, 'p_key', Lat0.key{e.p}, 'q_key', Lat0.key{e.q}, ...
        'status', e.status, 'calls', e.calls, 'child_root', e.child_root, 'kind', e.kind, ...
        'rule', e.rule, 'h', e.h); %#ok<AGROW>
end
end

function L = param_diff(a, b)
L = {};
if ~isstruct(a) || ~isstruct(b), L = {'(not comparable)'}; return; end
for n = union(fieldnames(a), fieldnames(b)).'
    if ~isfield(a, n{1}) || ~isfield(b, n{1}) || ~isequaln(a.(n{1}), b.(n{1}))
        L{end+1} = n{1}; %#ok<AGROW>
    end
end
if isempty(L), L = {'(field order)'}; end
end

% =========================================================================
% stamp

function [h, parts] = stamp(src)
top = src;
[up, nm] = fileparts(regexprep(src, '[\\/]+$', ''));
if strcmp(nm, 'internal') && exist(fullfile(src, 'settings_effective.mat'), 'file') == 2, top = up; end
D = sa_layout('state', top);
names = {'settings_effective.mat', 'attempts.mat', 'queue_final.mat', 'journal.txt', 'compaction.mat'};
parts = struct();
txt = '';
for k = 1:numel(names)
    f = fullfile(D, names{k});
    if exist(f, 'file') == 2, v = sa_sha256(f); else, v = ''; end
    parts.(matlab.lang.makeValidName(names{k})) = v;
    txt = [txt names{k} '=' v ';']; %#ok<AGROW>
end
f = dir(fullfile(D, 'attempts', 'A*.mat'));
n = sum(~cellfun(@isempty, regexp({f.name}, '^A\d{6}\.mat$', 'once')));
parts.attempt_files = n;
txt = [txt sprintf('attempt_files=%d;', n)];
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(txt), 'int8'));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
end

function h = current_stamp(S, P)
h = '';
if ~(isfield(S, 'import') && isstruct(S.import) && ~isempty(S.import.source_dir)), return; end
src = S.import.source_dir;
if exist(src, 'dir') ~= 7, return; end
try
    top = src;
    [up, nm] = fileparts(regexprep(src, '[\\/]+$', ''));
    if strcmp(nm, 'internal'), top = up; end
    D = sa_layout('state', top);
    if exist(fullfile(D, 'settings_effective.mat'), 'file') ~= 2, return; end
    if nargin >= 2 && isstruct(P) && isfield(P, 'source_stopped')
        % a killed source: its planned attempt files, re-hashed
        rows = P.rows;
        for k = 1:numel(rows)
            if ~strcmp(rows(k).status, 'planned'), continue; end
            f = fullfile(D, 'attempts', [rows(k).source_attempt_id '.mat']);
            if exist(f, 'file') ~= 2, return; end
            rows(k).source_attempt_sha256 = file_sha(f);
        end
        h = rows_stamp(rows);
        return
    end
    h = stamp(src);
catch
    h = '';
end
end

function h = file_sha(f)
% SHA-256 of a file, cached in this MATLAB process by (path, bytes,
% modification time in ms): a status / prepare sequence hashes each planned
% attempt file of a killed source once; a rewritten or replaced file has a
% new modification time and is hashed again (each import call checks its
% file's SHA-256 again in any case)
persistent C
if isempty(C), C = containers.Map('KeyType', 'char', 'ValueType', 'char'); end
jf = java.io.File(f);
key = sprintf('%s|%d|%d', f, jf.length(), jf.lastModified());
if isKey(C, key), h = C(key); return; end
h = sa_sha256(f);
if C.Count > 200000, C = containers.Map('KeyType', 'char', 'ValueType', 'char'); end
C(key) = h;
end

function h = name_stamp(imp)
% the stamp in a default run-folder name: the folder stamp, or for a source
% marked stopped the SHA-256 of its settings_effective.mat, which stays the
% same while the source is resumed
if isstruct(imp) && isfield(imp, 'source_stopped') && isequal(imp.source_stopped, true)
    src = imp.source_dir;
    top = src;
    [up, nm] = fileparts(regexprep(src, '[\\/]+$', ''));
    if strcmp(nm, 'internal'), top = up; end
    h = sa_sha256(fullfile(sa_layout('state', top), 'settings_effective.mat'));
else
    h = stamp(imp.source_dir);
end
end

function f = fence_of(S, P)
% the family fence in force ({} = none)
f = {};
if ~(isfield(S, 'continuation') && isfield(S.continuation, 'families')), return; end
v = S.continuation.families;
if ischar(v)
    if nargin >= 2 && isstruct(P) && isfield(P, 'fence'), f = P.fence; end
else
    f = reshape(v, 1, []);
end
end

% =========================================================================
% the call

function [u, tg] = unit_of(ctx, row)
P = ctx.import;
u = sa_search('unit');
u.route = 'import'; u.kind = 'import'; u.phase = 'import';
u.target_number = row.i; u.family = row.family; u.key = row.key; u.index = row.index;
u.fraction = 1; u.coordinate = row.coordinate;
u.guess_source = sprintf('import %s %s', P.stamp(1:8), row.source_attempt_id);
if strcmp(row.kind, 'mirror'), u.guess_source = [u.guess_source ' mirror']; end
if strcmp(P.guess_mesh, 'rule'), u.guess_source = [u.guess_source sa_guess_mesh('tag', ctx.S)]; end
u.seed_left = row.source_root_id;
u.seed_right = row.source_attempt_id;
tg = struct('family', row.family, 'left_type', '', 'right_type', '', 'node', row.node, ...
    'index', row.index, 'coordinate', row.coordinate, 'key', row.key, 'cell_left', [NaN NaN], ...
    'cell_right', [NaN NaN]);
end

function g = guess_of(ctx, row)
src = ctx.S.import.source_dir;
top = src;
[up, nm] = fileparts(regexprep(src, '[\\/]+$', ''));
if strcmp(nm, 'internal'), top = up; end
f = fullfile(sa_layout('state', top), 'attempts', [row.source_attempt_id '.mat']);
if exist(f, 'file') ~= 2
    integrity(ctx, sprintf('import row %d: the source attempt file %s is missing', row.i, f));
end
h = sa_sha256(f);
if ~strcmp(h, row.source_attempt_sha256)
    integrity(ctx, sprintf('import row %d: the source attempt file %s differs from the plan (SHA-256 %s, planned %s)', ...
        row.i, row.source_attempt_id, h(1:12), row.source_attempt_sha256(1:min(12, end))));
end
Z = load(f, 'sol');
sol = Z.sol;
if strcmp(row.kind, 'mirror')
    p = ctx.M.reflection.state_perm;
    sol.y = sol.y(p, :);
    if isfield(sol, 'yp') && ~isempty(sol.yp), sol.yp = sol.yp(p, :); end
end
if strcmp(ctx.import.guess_mesh, 'rule')
    g = sa_guess_mesh('from', sol, ctx.S);
else
    g = struct('x', sol.x, 'y', sol.y);
end
end

function integrity(ctx, why)
try
    sa_journal('note', ctx.state_dir, ['integrity ' why]);
catch
end
error('SA:INTEGRITY', 'Integrity stop at attempt count %d: %s', numel(ctx.attempts), why);
end

% =========================================================================
% source-edge knowledge

function K = knowledge(ctx)
K = struct('active', false, 'mode', 'none', 'src2new', [], 'new2src', [], 'pass_a', [], 'pass_b', [], ...
    'lin_ss', [], 'lin_weak', [], 'h_pass', NaN, 'h_lin', NaN);
if ~is_active(ctx) || strcmp(ctx.import.source_edges_mode, 'none'), return; end
P = ctx.import;
S = ctx.S;
K.active = true; K.mode = P.source_edges_mode;
K.h_lin = S.continuation.max_halvings;
K.h_pass = S.root_pass.max_halvings;
mk = @() containers.Map('KeyType', 'char', 'ValueType', 'logical');
K.src2new = containers.Map('KeyType', 'char', 'ValueType', 'char');
K.new2src = containers.Map('KeyType', 'char', 'ValueType', 'any');
K.pass_a = mk(); K.pass_b = mk(); K.lin_ss = mk(); K.lin_weak = mk();
A = ctx.attempts;
direct = containers.Map('KeyType', 'double', 'ValueType', 'logical');
for k = 1:numel(P.rows), direct(P.rows(k).i) = strcmp(P.rows(k).kind, 'direct'); end
for k = 1:numel(A)
    if ~strcmp(A(k).kind, 'import') || ~strcmp(A(k).outcome, 'accepted'), continue; end
    a = A(k);
    i = row_number(a);
    if isnan(i) || ~direct(i), continue; end
    s = P.rows(i).source_root_id;
    K.src2new(s) = a.root_id;
    if isKey(K.new2src, a.root_id), K.new2src(a.root_id) = [K.new2src(a.root_id), {s}];
    else, K.new2src(a.root_id) = {s}; end
end
E = P.edges;
for k = 1:numel(E)
    e = E(k);
    if e.calls <= 0 || isempty(e.parent_root) || ~isKey(K.src2new, e.parent_root), continue; end
    % a secant rescue edge's guess is extrapolated from two lineage roots, not stepped from the parent's
    % iterate: its refusal says nothing about an edge from that parent (G22 finding at V1; L = 30 has 161)
    if strcmp(e.kind, 'rescue_secant'), continue; end
    if ~isKey(ctx.node_of, e.q_key) || ~isKey(ctx.node_of, e.p_key), continue; end
    np = K.src2new(e.parent_root);
    q = ctx.node_of(e.q_key); p = ctx.node_of(e.p_key);
    switch e.status
        case 'completed'
            if ~isempty(e.child_root) && isKey(K.src2new, e.child_root)
                K.pass_a(sprintf('%s|%d', np, q)) = true;
                K.pass_b(sprintf('%s|%d', K.src2new(e.child_root), p)) = true;
            end
        case 'refused'
            if strcmp(K.mode, 'all') || e.h >= K.h_pass
                K.pass_a(sprintf('%s|%d', np, q)) = true;
            end
            if e.h >= K.h_lin
                K.lin_ss(sprintf('%s|%d', np, q)) = true;
            else
                K.lin_weak(sprintf('%s|%d', np, q)) = true;
            end
    end
end
end

function i = row_number(a)
% the plan row of an import attempt: its unit key's target number
t = regexp(a.unit_key, '^import\|import\|import\|p\d+\|t([0-9.eE+-]+)\|', 'tokens', 'once');
if isempty(t), i = NaN; else, i = str2double(t{1}); end
end

% =========================================================================
% names and families

function f = family_of(label, seed_family, types, map)
f = '';
if ~isempty(label), f = sa_mode_name('library_family', label, types, map); end
if isempty(f), f = seed_family; end
end

function m = mirror_label(lab)
m = lab;
t = regexp(lab, '^(0|[1-9][0-9]*[+-])/(0|[1-9][0-9]*[+-])$', 'tokens', 'once');
if ~isempty(t)
    m = [flip_half(t{2}) '/' flip_half(t{1})];
    return
end
m = mirror_family(lab);
end

function h = flip_half(h)
if strcmp(h, '0'), return; end
n = str2double(h(1:end-1));
if mod(n, 2) == 1
    if h(end) == '+', h(end) = '-'; else, h(end) = '+'; end
end
end

function f = mirror_family(f)
parts = strsplit(f, '_');
if numel(parts) == 2 && ~isempty(parts{1}) && ~isempty(parts{2})
    f = [parts{2} '_' parts{1}];
end
end

% =========================================================================
% report

function L = report_lines(ctx)
L = {};
if ~is_active(ctx), return; end
P = ctx.import; A = ctx.attempts;
L{end+1} = '';
L{end+1} = sprintf('Import (roots of an earlier run as seeds): source %s', P.source_dir);
L{end+1} = sprintf('  source stamp %s; source software %s', P.stamp, P.source_software);
L{end+1} = sprintf('  source model version %s', P.source_model_version);
if isempty(P.differing_leaves)
    L{end+1} = '  scientific settings differing from the source: none';
else
    L{end+1} = sprintf('  scientific settings differing from the source (information): %s', ...
        strjoin(P.differing_leaves, ', '));
end
L{end+1} = sprintf('  guess mesh %s; mirror images %d; source edges %s; plan %.1f s (hashing %.1f s, %.3g GB)', ...
    P.guess_mesh, P.mirror, P.source_edges_mode, P.plan_s, P.hash_s, P.hash_bytes / 1e9);
if isfield(P, 'source_edges_auto')
    if isempty(P.source_edges_auto.why)
        L{end+1} = '  source edges auto -> same_or_stronger (same lattice, solver, validity, step rule and model version)';
    else
        L{end+1} = sprintf('  source edges auto -> none: %s', strjoin(P.source_edges_auto.why, '; '));
    end
end
if isfield(P, 'source_stopped')
    K = P.source_stopped;
    L{end+1} = sprintf(['  source KILLED (import.source_stopped): last call started %s; %d call(s) without an ' ...
        'end line ignored%s; %d root(s) whose call has no end line ignored; %d attempt files'], ...
        tern(isempty(K.last_started), '(none)', K.last_started), numel(K.open_calls), ...
        tern(isempty(K.open_calls), '', [' (' strjoin(K.open_calls, ', ') ')']), numel(K.ignored_roots), ...
        K.attempt_files);
    L{end+1} = sprintf(['  checkpoint used: queue.mat at call count %g (source-edge knowledge from it), ' ...
        'attempts.mat at %g; identity: SHA-256 of the planned attempt files (the source may be resumed)'], ...
        K.checkpoint_queue, K.checkpoint_attempts);
end
if isfield(P, 'selection')
    G = P.selection;
    if isempty(G.root_ids_file), src_ids = sprintf('%d id(s)', numel(G.root_ids));
    else, src_ids = sprintf('%d id(s) from %s', numel(G.root_ids), G.root_ids_file);
    end
    L{end+1} = sprintf(['  selection: import.root_ids %s; import.families {%s}: selected %d of %d source ' ...
        'roots, skipped %d; selected labels {%s}'], src_ids, strjoin(G.families, ', '), G.n_selected, ...
        G.n_source, G.n_skipped, strjoin(G.labels, ', '));
end
if isfield(P, 'fence')
    L{end+1} = sprintf('  fence (continuation.families ''selected''): {%s}', strjoin(P.fence, ', '));
end
imp = A(strcmp({A.kind}, 'import'));
for kind = {'direct', 'mirror'}
    R = P.rows(strcmp({P.rows.kind}, kind{1}));
    if isempty(R) && strcmp(kind{1}, 'mirror'), continue; end
    ii = [R.i];
    a = imp(ismember(arrayfun(@row_number, imp), ii));
    acc = a(strcmp({a.outcome}, 'accepted'));
    L{end+1} = sprintf(['  %s: planned %d; skipped no_node_near %d, outside_box %d, mirror_twin_planned %d; ' ...
        'called %d; accepted %d (new root %d, matched %d); invalid %d; no_return %d'], kind{1}, ...
        sum(strcmp({R.status}, 'planned')), sum(strcmp({R.status}, 'no_node_near')), ...
        sum(strcmp({R.status}, 'outside_box')), sum(strcmp({R.status}, 'mirror_twin_planned')), numel(a), ...
        numel(acc), sum(strcmp({acc.store_status}, 'new_root')), sum(~strcmp({acc.store_status}, 'new_root')), ...
        sum(strcmp({a.outcome}, 'invalid')), sum(strcmp({a.outcome}, 'no_return'))); %#ok<AGROW>
end
if ~isempty(imp)
    w = [imp.wall_s]; m = [imp.returned_mesh];
    L{end+1} = sprintf('  re-solve wall time per call: median %.3g s, 90th percentile %.3g s, max %.3g s', ...
        median(w), prctile_(w, 90), max(w));
    edges = [0 1000 2000 5000 10000 20000 Inf];
    for b = 1:numel(edges) - 1
        s = w(m > edges(b) & m <= edges(b+1));
        if isempty(s), continue; end
        L{end+1} = sprintf('    returned mesh (%g, %g]: %d calls, median %.3g s, p90 %.3g s, max %.3g s', ...
            edges(b), edges(b+1), numel(s), median(s), prctile_(s, 90), max(s)); %#ok<AGROW>
    end
    nm = sum(isnan(m));
    if nm > 0, L{end+1} = sprintf('    no returned mesh: %d calls', nm); end
    rf = imp(~strcmp({imp.outcome}, 'accepted'));
    for k = 1:numel(rf)
        i = row_number(rf(k));
        L{end+1} = sprintf('    refused: %s %s -> %s (%s, %s) %s', rf(k).attempt_id, P.rows(i).source_root_id, ...
            P.rows(i).key, P.rows(i).source_label, P.rows(i).kind, rf(k).failing); %#ok<AGROW>
    end
    acc = imp(strcmp({imp.outcome}, 'accepted'));
    mm = 0;
    for k = 1:numel(acc)
        i = row_number(acc(k));
        if ~strcmp(acc(k).observed_family, P.rows(i).source_label), mm = mm + 1; end
    end
    L{end+1} = sprintf('  label mismatches (information; labels are descriptions): %d of %d accepted', mm, numel(acc));
end
pl = P.rows(strcmp({P.rows.status}, 'planned'));
if ~isempty(pl)
    L{end+1} = sprintf('  largest map distance %.6g new-lattice steps; exact %d of %d planned', ...
        max([pl.map_distance]), sum([pl.exact]), numel(pl));
end
c = P.source_node_counts;
[u, ~, j] = unique(c);
L{end+1} = sprintf(['  source roots per node: %s; source nodes at or above its multiseed.max_roots_per_node (%g): ' ...
    '%d of %d'], strjoin(arrayfun(@(a, b) sprintf('%d:%d', a, b), u, accumarray(j(:), 1).', ...
    'UniformOutput', false), ' '), P.source_max_roots_per_node, P.source_nodes_at_cap, numel(c));
n = [];
if isfield(ctx, 'nroots'), n = ctx.nroots(ctx.nroots > 0); end
if ~isempty(n)
    [u, ~, j] = unique(n(:).');
    L{end+1} = sprintf('  roots per node in this run: %s', strjoin(arrayfun(@(a, b) sprintf('%d:%d', a, b), u, ...
        accumarray(j(:), 1).', 'UniformOutput', false), ' '));
end
if isstruct(ctx.Q) && isfield(ctx.Q, 'source_skips')
    sk = ctx.Q.source_skips;
    cats = {'pass_a', 'pass_b', 'lineage_same_or_stronger', 'lineage_weaker'};
    L{end+1} = sprintf('  source-edge skips (%s): %s', P.source_edges_mode, strjoin(cellfun(@(c) ...
        sprintf('%s %d', c, sum(strcmp({sk.category}, c))), cats, 'UniformOutput', false), ', '));
end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function v = prctile_(x, p)
x = sort(x(:));
if isempty(x), v = NaN; return; end
k = max(1, min(numel(x), ceil(p / 100 * numel(x))));
v = x(k);
end
