function varargout = sa_store_io(action, varargin)
%SA_STORE_IO  Save/load wrapper around the copied root-store kernel.
%
%   SA_STORE_IO('save', store, file)
%       Atomic write of a (monolithic, legacy) store (sa_write).
%   SA_STORE_IO('save_shard', store, file, key, attempt_count)
%       Atomic write of one node's shard: the kernel
%       store of node `key`, stamped with the run's attempt count at the time
%       of writing (a checkpoint is compared at its own count).
%       Profile deduplication (storage only; format 'sa_store_io/shard-2'): a discovery's
%       returned_profile / comparison_profile that is bit-identical to its
%       root's representative_profile (same class, size and IEEE bits) is
%       written as [] with a flag (profile_dedup{r}(1|2, q)); every reader
%       below restores it, so the in-memory kernel store is isequaln to the
%       one saved. The kernel, its fields and every verdict are untouched.
%       Older 'sa_store_io/shard-1' files load unchanged.
%   S = SA_STORE_IO('read_shard', file)
%       The shard file's variables with S.store DECODED and NOT
%       re-validated: the form a checkpoint is compared in (sa_search) and
%       the one analysis code must read (a raw load() of a shard-2 file
%       shows [] for every deduplicated profile).
%   store = SA_STORE_IO('load', file, M, rederive)
%   [store, key, count, prof] = SA_STORE_IO('load_shard', file, M, rederive [, profiles])
%       Load and re-validate. The kernel validates a store only inside 'add',
%       so the loaded store is REBUILT: every discovery is replayed, in
%       discovery order, into an empty store with the saved settings, and the
%       rebuilt store must equal the loaded one exactly. Every
%       discovery's verdict is RE-DERIVED, never read: rederive is a function
%       [valid, failing] = rederive(discovery) that re-evaluates the stored
%       attempt (sa_validity 'rederive'); without it the load is refused.
%       Any disagreement is an error (fail closed).
%       Profile check: load_shard also RECOMPUTES every
%       discovery's returned mesh, returned profile and comparison profile
%       from its attempt file's sol under the run's identity.candidate_mesh
%       (sa_identity 'candidate' + the kernel's interpolation) and compares
%       bits; a mismatch is SA:STORE_IO:E_PROFILE. Version scope as for verdicts
%       re-derivation: a record written under another software or MATLAB version or
%       thread mode is not failed; its largest absolute difference is
%       reported in prof (attempt_id, same_version, max_abs_diff, equal).
%       profiles: a function from 'run_profiles'; when omitted it is built
%       from the run folder holding the shard (<run>/roots/<file>), and the
%       check is skipped (prof.checked = false) only for a file outside a
%       run folder. read_shard stays raw (analysis only, never verified).
%   name = SA_STORE_IO('shard_file', key)    roots/<name>.mat for a node key
%   f = SA_STORE_IO('shard_files', run_dir)
%       The shard files of a state folder (dir entries of roots/*.mat)
%       without sa_write's temporary siblings. An interrupted save leaves
%       <name>.mat.tmp_<pid> (MATLAB R2025b's save keeps a name that does
%       not end in .mat as given), which *.mat does not match; the filter
%       is a guard for a save that appends .mat (<name>.mat.tmp_<pid>.mat,
%       sa_write's fallback branch): a name holding '.tmp_' is skipped (a
%       shard name never holds a dot). Used by resume prepare,
%       load_run_shards, seeded_atlas_roots, seeded_atlas_plot,
%       seeded_atlas_compact and the pilot report, which never read such a
%       file; resume prepare lists every *.tmp_* file as a stray temporary
%       file (resume_<n>.mat, journal note).
%   fn = SA_STORE_IO('run_rederive', run_dir, M)
%       The verdict re-derivation function of a run folder: loads the
%       discovery's attempt file, rebuilds the options from the effective
%       settings and the attempt's kind, re-evaluates (sa_validity
%       'rederive'); an integrity failure is an error (SA:INTEGRITY); for a
%       same-version record the discovery's stored validity must also equal
%       the re-derived one.
%   shards = SA_STORE_IO('load_run_shards', run_dir, M)
%       containers.Map node key -> re-validated kernel store, every shard of
%       the run folder (verdict re-derivation and the profile check).
%   fn = SA_STORE_IO('run_profiles', run_dir, M)   the profile recomputation
%   store = SA_STORE_IO('light', store)
%       The in-memory form the executor keeps between calls (memory only):
%       every discovery's returned_mesh and returned_profile (the n_dense
%       grid) set to []. Everything else is kept: representative and
%       comparison profiles, ids, routes, verdicts, comparisons. The kernel
%       refuses a light store (malformed discovery), so it can never be
%       added to or saved by mistake; shard FILES are always written from
%       the full store and are unchanged.
%   store = SA_STORE_IO('full', store, M, S, K, sol_of)
%       The full store again: each light discovery's returned_mesh and
%       returned_profile are recomputed from its attempt's sol (sol_of: a
%       function attempt_id -> sol, e.g. sa_search 'get_sol') with
%       sa_identity('returned'), the expressions sa_identity('candidate')
%       built them with, so the bits are the original ones (the check
%       load_shard's E_PROFILE makes). Discoveries already full are left
%       as they are. Recomputing from the attempt file, never reading the
%       node's shard file, is deliberate: a kill between an attempt write
%       and its shard write leaves the shard file one discovery behind.
%   tf = SA_STORE_IO('is_light', store)   true when a discovery is light
%   F = SA_STORE_IO('full_map', result)
%       result.shards (a seeded_atlas_run result, light stores) as full
%       stores: containers.Map node key -> 'full' of each store, the sols
%       read from the run folder's attempt files, the store settings from
%       its settings_effective.mat (memory: the whole store in full form)

switch action
    case 'save',       save_store(varargin{:});
    case 'save_shard', save_shard(varargin{:});
    case 'load',       varargout{1} = load_store(varargin{:});
    case 'load_shard', [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = load_shard(varargin{:});
    case 'run_profiles', varargout{1} = run_profiles(varargin{:});
    case 'read_shard', varargout{1} = read_shard(varargin{:});
    case 'shard_file', varargout{1} = shard_file(varargin{:});
    case 'shard_files', varargout{1} = shard_files(varargin{:});
    case 'run_rederive', varargout{1} = run_rederive(varargin{:});
    case 'load_run_shards', varargout{1} = load_run_shards(varargin{:});
    case 'light',      varargout{1} = light_store(varargin{:});
    case 'full',       varargout{1} = full_store(varargin{:});
    case 'is_light',   varargout{1} = is_light(varargin{:});
    case 'full_map',   varargout{1} = full_map(varargin{:});
    otherwise, error('SA:STORE_IO', 'Unknown sa_store_io action %s.', action);
end
end

function store = light_store(store)
% memory only: the n_dense returned mesh and profile of every discovery
for r = 1:numel(store.roots)
    for q = 1:numel(store.roots(r).discoveries)
        store.roots(r).discoveries(q).returned_mesh = [];
        store.roots(r).discoveries(q).returned_profile = [];
    end
end
end

function store = full_store(store, M, S, K, sol_of)
for r = 1:numel(store.roots)
    for q = 1:numel(store.roots(r).discoveries)
        if ~isempty(store.roots(r).discoveries(q).returned_mesh), continue; end
        [mesh, prof] = sa_identity('returned', M, sol_of(store.roots(r).discoveries(q).attempt_id), S, K);
        store.roots(r).discoveries(q).returned_mesh = mesh;
        store.roots(r).discoveries(q).returned_profile = prof;
    end
end
end

function F = full_map(res)
M = sa_model(res.settings.model);
D = sa_layout('state', res.run_dir);
E = load(fullfile(D, 'settings_effective.mat'), 'store_settings');
F = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = res.shards.keys
    F(k{1}) = full_store(res.shards(k{1}), M, res.settings, E.store_settings, @(id) attempt_sol(D, id));
end
end

function sol = attempt_sol(D, id)
Z = load(fullfile(D, 'attempts', [id '.mat']), 'sol');
sol = Z.sol;
end

function tf = is_light(store)
tf = false;
for r = 1:numel(store.roots)
    if any(arrayfun(@(d) isempty(d.returned_mesh), store.roots(r).discoveries)), tf = true; return; end
end
end

function f = shard_file(key)
f = [regexprep(key, '[^A-Za-z0-9\-]', '_') '.mat'];
end

function f = shard_files(run_dir)
% roots/*.mat without sa_write's temporary siblings (a guard: R2025b writes
% <name>.mat.tmp_<pid>, which *.mat does not match; a save that appended
% .mat would write <name>.mat.tmp_<pid>.mat)
f = dir(fullfile(run_dir, 'roots', '*.mat'));
f = f(cellfun(@isempty, strfind({f.name}, '.tmp_')));
end

function save_store(store, file)
sa_write('mat', file, struct('store', store, 'format', 'sa_store_io/1'));
end

function save_shard(store, file, key, count)
[enc, flags] = dedup_encode(store);
sa_write('mat', file, struct('store', enc, 'format', 'sa_store_io/shard-2', ...
    'key', key, 'attempt_count', count, 'profile_dedup', {flags}));
end

function S = read_shard(file)
S = load(file);
if ~isfield(S, 'format') || ~isfield(S, 'store') || ~ischar(S.format)
    error('SA:STORE_IO:E_FORMAT', '%s is not an sa_store_io shard file.', file);
end
switch S.format
    case 'sa_store_io/shard-1'                  % older format: as written
    case 'sa_store_io/shard-2'
        if ~isfield(S, 'profile_dedup')
            error('SA:STORE_IO:E_FORMAT', '%s: shard-2 file without profile_dedup.', file);
        end
        S.store = dedup_decode(S.store, S.profile_dedup, file);
    otherwise
        error('SA:STORE_IO:E_FORMAT', '%s is not an sa_store_io shard file (%s).', file, S.format);
end
end

% ---- profile deduplication (serialization only) ---------------------------
function [store, flags] = dedup_encode(store)
names = {'returned_profile', 'comparison_profile'};
flags = cell(1, numel(store.roots));
for r = 1:numel(store.roots)
    rep = store.roots(r).representative_profile;
    f = false(2, numel(store.roots(r).discoveries));
    for q = 1:numel(store.roots(r).discoveries)
        for k = 1:2
            if same_bits(store.roots(r).discoveries(q).(names{k}), rep)
                store.roots(r).discoveries(q).(names{k}) = [];
                f(k, q) = true;
            end
        end
    end
    flags{r} = f;
end
end

function store = dedup_decode(store, flags, file)
names = {'returned_profile', 'comparison_profile'};
if ~isfield(store, 'roots') || ~iscell(flags) || numel(flags) ~= numel(store.roots)
    error('SA:STORE_IO:E_DEDUP', '%s: profile_dedup does not match the roots.', file);
end
for r = 1:numel(store.roots)
    f = flags{r};
    n = numel(store.roots(r).discoveries);
    if ~islogical(f) || ~isequal(size(f), [2 n])
        error('SA:STORE_IO:E_DEDUP', '%s: profile_dedup{%d} is not a 2 x %d logical.', file, r, n);
    end
    rep = store.roots(r).representative_profile;
    for q = 1:n
        for k = 1:2
            v = store.roots(r).discoveries(q).(names{k});
            if f(k, q)
                if ~isa(v, 'double') || ~isequal(size(v), [0 0]) || isempty(rep)
                    error('SA:STORE_IO:E_DEDUP', '%s: root %d discovery %d %s is flagged but stored.', ...
                        file, r, q, names{k});
                end
                store.roots(r).discoveries(q).(names{k}) = rep;
            elseif isempty(v)
                error('SA:STORE_IO:E_DEDUP', '%s: root %d discovery %d %s is empty but not flagged.', ...
                    file, r, q, names{k});
            end
        end
    end
end
end

function tf = same_bits(a, b)
% bit identity (stricter than isequal: -0 ~= +0, so a restore is bit-exact)
tf = isa(a, 'double') && isa(b, 'double') && isreal(a) && isreal(b) && ...
    ~issparse(a) && ~issparse(b) && ~isempty(a) && isequal(size(a), size(b)) && ...
    isequal(typecast(a(:), 'uint64'), typecast(b(:), 'uint64'));
end

function store = load_store(file, M, rederive)
if nargin < 3, rederive = []; end
S = load(file);
if ~isfield(S, 'format') || ~strcmp(S.format, 'sa_store_io/1') || ~isfield(S, 'store')
    error('SA:STORE_IO:E_FORMAT', '%s is not an sa_store_io/1 file.', file);
end
store = rebuild(S.store, file, M, rederive);
end

function [store, key, count, prof] = load_shard(file, M, rederive, profiles)
if nargin < 3, rederive = []; end
if nargin < 4
    profiles = [];
    rd = fileparts(fileparts(file));
    if exist(fullfile(rd, 'settings_effective.mat'), 'file') == 2 && exist(fullfile(rd, 'attempts'), 'dir') == 7
        profiles = run_profiles(rd, M);
    end
end
S = read_shard(file);                           % shard-1 as written; shard-2 decoded
key = S.key; count = S.attempt_count;
store = rebuild(S.store, file, M, rederive);
if ~all(arrayfun(@(r) strcmp(r.target_key, key), store.roots))
    error('SA:STORE_IO:E_SHARD', '%s holds a root of another node.', file);
end
prof = check_profiles(store, file, profiles);
end

% ---- profile check -------------------------------------------------------------
function prof = check_profiles(store, file, profiles)
prof = struct('checked', ~isempty(profiles), 'rows', struct('attempt_id', {}, 'same_version', {}, ...
    'max_abs_diff', {}, 'equal', {}));
if isempty(profiles), return; end
for r = 1:numel(store.roots)
    for q = 1:numel(store.roots(r).discoveries)
        d = store.roots(r).discoveries(q);
        P = profiles(d);
        eq = same_bits(d.returned_mesh, P.mesh) && same_bits(d.returned_profile, P.returned) && ...
            same_bits(d.comparison_profile, P.comparison);
        mx = max([maxdiff(d.returned_mesh, P.mesh), maxdiff(d.returned_profile, P.returned), ...
            maxdiff(d.comparison_profile, P.comparison)]);
        prof.rows(end+1) = struct('attempt_id', d.attempt_id, 'same_version', P.same_version, ...
            'max_abs_diff', mx, 'equal', eq);
        if ~eq && P.same_version
            error('SA:STORE_IO:E_PROFILE', ['%s: discovery %d (attempt %s) profiles differ from their ' ...
                'recomputation from the attempt''s sol (largest difference %g).'], file, ...
                d.discovery_number, d.attempt_id, mx);
        end
    end
end
end

function m = maxdiff(a, b)
if ~isequal(size(a), size(b)), m = Inf; return; end
m = max([0; abs(double(a(:)) - double(b(:)))]);
if any(isnan(a(:)) ~= isnan(b(:))), m = Inf; end
end

function fn = run_profiles(run_dir, M)
E = load(fullfile(run_dir, 'settings_effective.mat'), 'settings', 'store_settings');
prov = sa_settings('provenance');
fn = @(d) profile_one(run_dir, M, E.settings, E.store_settings, prov, d);
end

function P = profile_one(run_dir, M, S, K, prov, d)
rec = load(fullfile(run_dir, 'attempts', [d.attempt_id '.mat']), 'sol', 'target_row', 'provenance');
if strcmp(S.identity.candidate_mesh, 'comparison')
    mesh = K.comparison_mesh;
else
    mesh = linspace(rec.sol.x(1), rec.sol.x(end), S.validity.n_dense);
end
prof = deval(rec.sol, mesh, [M.identity_fields.row]);
cand = struct('target_key', rec.target_row.key, 'coordinate', d.coordinate, 'mesh', mesh, ...
    'profile', prof, 'seed_family', d.seed_family, 'observed_family', d.observed_family, ...
    'route', d.route, 'attempt_id', d.attempt_id, 'valid', true, 'validity', d.validity);
st = sa_root_store('add', sa_root_store('empty', K), cand);
P = struct('mesh', mesh, 'returned', prof, ...
    'comparison', st.roots(1).discoveries(1).comparison_profile, ...
    'same_version', isfield(rec, 'provenance') && isequal(rec.provenance, prov));
end

function store = rebuild(store, file, M, rederive) %#ok<INUSL>
if isempty(rederive) || ~isa(rederive, 'function_handle')
    error('SA:STORE_IO:E_REDERIVE', ['AM-1: a reload re-derives every discovery''s verdict ' ...
        'from its attempt record; supply the rederive function.']);
end
rebuilt = sa_root_store('empty', store.settings);
disc = [];
for r = 1:numel(store.roots)
    for q = 1:numel(store.roots(r).discoveries)
        d = store.roots(r).discoveries(q);
        disc = [disc; struct('n', d.discovery_number, 'root', r, 'd', d)]; %#ok<AGROW>
    end
end
if ~isempty(disc)
    [~, ord] = sort([disc.n]);
    disc = disc(ord);
end
for k = 1:numel(disc)
    d = disc(k).d;
    [valid, failing] = rederive(d);
    if ~valid
        error('SA:STORE_IO:E_VALIDITY', ['Discovery %d (attempt %s) fails the re-derived ' ...
            'validity verdict on reload: %s.'], d.discovery_number, d.attempt_id, strjoin(failing, ', '));
    end
    root = store.roots(disc(k).root);
    cand = struct('target_key', root.target_key, 'coordinate', d.coordinate, ...
        'mesh', d.returned_mesh, 'profile', d.returned_profile, ...
        'seed_family', d.seed_family, 'observed_family', d.observed_family, ...
        'route', d.route, 'attempt_id', d.attempt_id, 'valid', true, ...
        'validity', d.validity);
    [rebuilt, ~] = sa_root_store('add', rebuilt, cand);
end
if ~isequaln(rebuilt, store)
    error('SA:STORE_IO:E_REPLAY', ['%s: the stored root store differs from its ' ...
        'replay through the kernel.'], file);
end
end

function fn = run_rederive(run_dir, M)
E = load(fullfile(run_dir, 'settings_effective.mat'), 'settings', 'theta');
E.settings = sa_settings('legacy_fill', E.settings);   % neutral values for older folders
prov = sa_settings('provenance');
fn = @(d) rederive_one(run_dir, M, E.settings, E.theta, prov, d);
end

function [valid, failing] = rederive_one(run_dir, M, S, theta, prov, d)
rec = load(fullfile(run_dir, 'attempts', [d.attempt_id '.mat']));
if isfield(rec, 'unit'), kind = rec.unit.kind; else, kind = rec.attempt.kind; end
ctx = struct('M', M, 'theta', theta, 'S', S);
opts = sa_search('opts_for', ctx, sa_search('solver_for', S, kind));
[V, info] = sa_validity('rederive', M, rec, theta, S, opts, prov);
if ~info.integrity_ok
    error('SA:INTEGRITY', 'Re-derivation of %s: %s', d.attempt_id, info.reason);
end
valid = V.valid; failing = V.failing;
if info.same_version && ~isequaln(d.validity, V)
    valid = false;
    failing = [failing, {'stored discovery validity differs from its re-derivation'}];
end
end

function shards = load_run_shards(run_dir, M)
fn = run_rederive(run_dir, M);
pf = run_profiles(run_dir, M);
shards = containers.Map('KeyType', 'char', 'ValueType', 'any');
f = shard_files(run_dir);
for k = 1:numel(f)
    [st, key] = load_shard(fullfile(f(k).folder, f(k).name), M, fn, pf);
    shards(key) = st;
end
end
