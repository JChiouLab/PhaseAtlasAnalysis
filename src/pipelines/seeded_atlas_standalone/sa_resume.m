function varargout = sa_resume(action, varargin)
%SA_RESUME  Restart contract.
%
%   replay = SA_RESUME('prepare', run_dir, S, current)
%       run_dir  the run folder (the top); its state folder (sa_layout) is
%                read and written
%       S        the resolved settings of the resuming call (resume = true)
%       current  struct of the current run's identities: provenance,
%                code_hashes, library_files, model_version, parameters
%       0. A compacted folder (internal/compaction.mat) is refused,
%          SA:RESUME:E_COMPACTED.
%       1. Preconditions: the parameters equal the stored ones (isequaln,
%          SA:RESUME:E_PARAMETERS, when current carries them); every
%          scientific settings leaf equals the stored
%          effective settings (only budget.*, checkpoint.*, resume,
%          resume_code_override, resolved_at and run_dir -- the folder's
%          identity, recorded in resume_<n>.mat -- may differ); software
%          version, MATLAB version, architecture and thread mode equal
%          (never waived); model version (a toolbox change named as such)
%          and library SHA-256 equal; the import source's stamp equals
%          the stored one (SA:RESUME:E_IMPORT; checked when the
%          caller carries current.import_check = true (the stamp is
%          then computed here) or current.import_stamp; for a killed source
%          imported with import.source_stopped the stamp is recomputed over
%          its planned attempt files, so resuming that source does not
%          block this run); for a focus run (when current
%          carries focus_inputs or focus_error) the field scale and sheet
%          (family label) lookup derived again from the source equal the stored
%          focus_inputs (SA:RESUME:E_FOCUS; also when the source's
%          settings_effective.mat cannot be read; a focus folder started
%          with a library by an older version is refused E_FOCUS before
%          E_SETTINGS and E_LIBRARY); code hashes equal, unless
%          S.resume_code_override holds a reason, which is journaled with
%          both hash sets. Otherwise SA:RESUME:E_* (a changed setting or a
%          repair is a new run).
%       2. Reconcile: a journal 'start' without 'end' or file becomes an
%          'interrupted' record (keeps its id; the unit is re-queued at its
%          original position under a new id with retry_of); a file without
%          'end' is complete and its 'end' line is appended; stray *.tmp_*
%          files (sa_write's temporary siblings of an interrupted save) are
%          listed, never read: the attempt and shard listings skip them;
%          ids must be contiguous. The
%          appended 'end' line carries the stored row's ' cut_by_cap' suffix,
%          as a live call writes it.
%       3. The replay plan: attempt files in id order, and the stamps of the
%          stored checkpoints (queue.mat, attempts.mat, every shard), each to
%          be compared at its own attempt count (sa_search).
%       4. resume_<n>.mat: budget, checkpoint, time, pid, attempt count,
%          override use, interrupted ids, stray temporary files, the import
%          source folder used.
%       replay.import: the stored import plan (effective.import), loaded and
%       never rebuilt (sa_import('blank') for a run without one).
%       Both settings structs are compared after sa_settings
%       'fill_focus_leaves' (the selection and fence leaves at their
%       neutral values where a struct made before them lacks them: a
%       folder's stored settings, or a caller passing them, such as
%       seeded_atlas_compact).
%   Z = SA_RESUME('status', run_dir, S, current)
%       What a run with settings S would do with run_dir, from the same
%       preconditions as 'prepare', without reconciliation and WITHOUT
%       WRITING. When current carries M and theta, S first goes through
%       sa_lattice('derive_step') (the run's own derivation). Exactly one
%       state, the first that matches, in this order:
%         new               run_dir missing or empty
%         legacy            a flat folder written by an older version: not
%                           resumable; it still loads
%         refused           a precondition differs (Z.id E_PARAMETERS,
%                           E_SETTINGS, E_VERSION, E_MODEL, E_LIBRARY,
%                           E_IMPORT, E_FOCUS, E_CODE, in that order;
%                           E_IMPORT and / or E_FOCUS alone leave a
%                           completed run completed, with a
%                           note; Z.items lists every
%                           differing item), or a journal integrity note
%                           follows the last end line (Z.id E_INTEGRITY),
%                           or a non-empty folder that holds no run
%         completed         queue_final.mat says completed and its attempt
%                           count equals the attempt files (a compacted
%                           folder: the count recorded before compaction);
%                           Z.compacted = true marks a compacted run, whose
%                           resume prepare refuses (E_COMPACTED)
%         budget_exhausted  stopped at budget_calls and the calls used are
%                           at least S.budget.max_calls
%         resumable         any other stop (budget, wall, disk, or a kill
%                           that left queue_final.mat stale or missing)
%       Z: state, id, message, items, compacted, layout, state_dir,
%       attempt_count, stop_reason, settings (S after the derivation).
%   cur = SA_RESUME('current', M, parameters, library_manifest [, S])
%       The current identities for 'status' (M, theta, provenance, code
%       hashes, library files, model version, parameters; with S also
%       import_check = true: status compares the import source's identity,
%       computing only the stamp the stored plan needs); zero solves. With
%       S.mode 'focus' no library is read: library_files is
%       empty, seed_modes is the source's, and focus_inputs (sa_import
%       'focus_inputs') or focus_error is carried for E_FOCUS.
%   L = SA_RESUME('scientific_diff', A, B)   differing scientific leaves
%   H = SA_RESUME('code_hashes', M)          SHA-256 of every sa_*.m, the
%                                            entry point and module file;
%                                            and of private/, solver/,
%                                            solver/private/, models/*/private/

switch action
    case 'prepare',         varargout{1} = prepare(varargin{:});
    case 'status',          varargout{1} = status(varargin{:});
    case 'current',         varargout{1} = current(varargin{:});
    case 'scientific_diff', varargout{1} = scientific_diff(varargin{:});
    case 'code_hashes',     varargout{1} = code_hashes(varargin{:});
    otherwise, error('SA:RESUME', 'Unknown sa_resume action %s.', action);
end
end

% =========================================================================

function H = code_hashes(M)
% Also private/, solver/, solver/private/ (the
% generated traced bvp4c copy and its helpers) and models/*/private/ (the
% forwarders)
pkg = fileparts(mfilename('fullpath'));
f = [dir(fullfile(pkg, 'sa_*.m')); dir(fullfile(pkg, 'seeded_atlas_run.m'))];
modf = which(['sa_model_' M.id]);
if ~isempty(modf)
    f = [f; dir(fullfile(fileparts(modf), '*.m'))];
end
f = [f; dir(fullfile(pkg, 'private', '*.m')); dir(fullfile(pkg, 'solver', '*.m')); ...
    dir(fullfile(pkg, 'solver', 'private', '*.m')); dir(fullfile(pkg, 'models', '*', 'private', '*.m'))];
H = struct('file', {}, 'sha256', {});
for k = 1:numel(f)
    p = fullfile(f(k).folder, f(k).name);
    H(end+1) = struct('file', strrep(p, [pkg filesep], ''), 'sha256', sa_sha256(p)); %#ok<AGROW>
end
[~, o] = sort({H.file});
H = H(o);
end

function t = toolbox_change(a, b)
% '' unless both model versions carry a solver: component and its toolbox
% bvp4c.m or deval.m SHA-256 prefix differs: then the change, e.g. 'bvp4c 597839cef613 -> 0123456789ab'
t = '';
ba = tag(a, 'toolbox_bvp4c'); bb = tag(b, 'toolbox_bvp4c');
da = tag(a, 'toolbox_deval'); db = tag(b, 'toolbox_deval');
if isempty(ba) || isempty(bb) || (strcmp(ba, bb) && strcmp(da, db)), return; end
t = sprintf('bvp4c %s -> %s', ba, bb);
if ~strcmp(da, db), t = sprintf('%s, deval %s -> %s', t, da, db); end
end

function v = tag(s, name)
v = regexp(s, [name '=([0-9a-f]+)'], 'tokens', 'once');
if isempty(v), v = ''; else, v = v{1}; end
end

function L = scientific_diff(A, B)
% free entries: a top-level name frees the whole group; a dotted path
% (import.source_dir) frees that one leaf (leaf_diff consults the list)
free = sa_settings('resume_free');
L = {};
fa = fieldnames(A); fb = fieldnames(B);
names = union(fa, fb);
for k = 1:numel(names)
    n = names{k};
    if any(strcmp(n, free)), continue; end
    if ~isfield(A, n) || ~isfield(B, n)
        L{end+1} = n; %#ok<AGROW>
    elseif isstruct(A.(n)) && isstruct(B.(n)) && isscalar(A.(n)) && isscalar(B.(n))
        sub = leaf_diff(A.(n), B.(n), n, free);
        L = [L, sub]; %#ok<AGROW>
    elseif ~isequaln(A.(n), B.(n)) || ~strcmp(class(A.(n)), class(B.(n)))
        L{end+1} = n; %#ok<AGROW>
    end
end
end

function L = leaf_diff(a, b, path, free)
if nargin < 4, free = {}; end
L = {};
names = union(fieldnames(a), fieldnames(b));
for k = 1:numel(names)
    n = names{k}; p = [path '.' n];
    if any(strcmp(p, free)), continue; end
    if ~isfield(a, n) || ~isfield(b, n)
        L{end+1} = p; %#ok<AGROW>
    elseif isstruct(a.(n)) && isstruct(b.(n)) && isscalar(a.(n)) && isscalar(b.(n))
        L = [L, leaf_diff(a.(n), b.(n), p, free)]; %#ok<AGROW>
    elseif ~isequaln(a.(n), b.(n)) || ~strcmp(class(a.(n)), class(b.(n)))
        L{end+1} = p; %#ok<AGROW>
    end
end
end

function [bad, item] = import_diff(E0, cur, S)
% E_IMPORT: the stored import (effective.import: plan,
% stamp) against the current source. Checked only when the caller asks for
% it: cur.import_check = true (the current stamp is computed here, once,
% and only the kind the stored plan needs) or cur.import_stamp (a stamp the
% caller computed; '' = no source set or unreadable). A killed source's
% plan (field source_stopped) is checked by its planned attempt files
% instead (sa_import 'current_stamp' with the stored plan; the folder stamp
% is not computed for it).
bad = false; item = '';
if ~isfield(cur, 'import_stamp') && ~(isfield(cur, 'import_check') && isequal(cur.import_check, true)), return; end
st = '';
if isfield(E0, 'import') && isstruct(E0.import) && isfield(E0.import, 'active') && E0.import.active
    st = E0.import.stamp;
end
killed = ~isempty(st) && isfield(E0.import, 'source_stopped');
if killed
    cs = sa_import('current_stamp', S, E0.import);
elseif isfield(cur, 'import_stamp')
    cs = cur.import_stamp;
else
    cs = sa_import('current_stamp', S);
end
cur.import_stamp = cs;
if ~strcmp(st, cur.import_stamp)
    bad = true;
    if isempty(st)
        item = 'import: the stored run has no import but import.source_dir is set';
    elseif killed && isempty(cs)
        item = sprintf(['import: a planned source attempt file of the killed source (identity %s) is missing, ' ...
            'or the source is not found at import.source_dir'], st(1:8));
    elseif killed
        item = sprintf(['import: a planned source attempt file of the killed source changed (identity %s -> %s; ' ...
            'attempt files are never rewritten, and resuming that source does not change them)'], st(1:8), cs(1:8));
    elseif isempty(cur.import_stamp)
        item = sprintf('import: the source run (stamp %s) is not found at import.source_dir', st(1:8));
    else
        item = sprintf(['import: the source run changed since this run started (stamp %s -> %s; ' ...
            'a source resumed or compacted after the follow-up started blocks its resume)'], st(1:8), ...
            cur.import_stamp(1:min(8, end)));
    end
end
end

% =========================================================================

function P = prepare(top, S, cur)
eid = 'SA:RESUME';
run_dir = sa_layout('state', top);
f = fullfile(run_dir, 'settings_effective.mat');
if exist(f, 'file') ~= 2
    error([eid ':E_NO_RUN'], 'resume: %s holds no settings_effective.mat.', run_dir);
end
if exist(fullfile(run_dir, 'compaction.mat'), 'file') == 2
    Zc = load(fullfile(run_dir, 'compaction.mat'));
    if isfield(Zc, 'compacted_from')
        st = sprintf('it was compacted after it stopped (%s) without completing', Zc.stop_reason);
    else
        st = 'it is complete';
    end
    error([eid ':E_COMPACTED'], ['resume refused: %s was compacted on %s: %d files removed; a ' ...
        'compacted run cannot be resumed (%s; choose a new run_dir for more work).'], ...
        top, Zc.compacted_at, Zc.n_removed, st);
end
E0 = load(f);
% ---- 1. preconditions ----------------------------------------------------
[legacy, item] = focus_legacy(E0, cur);
if legacy                                        % first: its library files would be named otherwise
    error([eid ':E_FOCUS'], 'resume refused: %s.', item);
end
if isfield(cur, 'parameters') && isfield(E0, 'parameters') && ~isequaln(E0.parameters, cur.parameters)
    error([eid ':E_PARAMETERS'], ['resume refused: the parameters differ from the stored run (%s); ' ...
        'other parameters are a new run: choose a new run_dir.'], ...
        strjoin(param_diff(E0.parameters, cur.parameters), ', '));
end
S = focus_seed_modes(S, E0, cur);
d = scientific_diff(sa_settings('fill_focus_leaves', E0.settings), sa_settings('fill_focus_leaves', S));
if ~isempty(d)
    error([eid ':E_SETTINGS'], ['resume refused: scientific setting(s) changed (%s); a ' ...
        'changed setting is a new run.'], strjoin(d, ', '));
end
if ~isfield(E0, 'provenance') || ~isequal(E0.provenance, cur.provenance)
    error([eid ':E_VERSION'], ['resume refused: software version, MATLAB version, architecture or ' ...
        'thread mode differ from the stored run (never waived).']);
end
if ~strcmp(E0.model.version, cur.model_version)
    t = toolbox_change(E0.model.version, cur.model_version);
    if isempty(t)
        error([eid ':E_MODEL'], 'resume refused: the model version differs.');
    end
    error([eid ':E_MODEL'], ['resume refused: the MATLAB toolbox changed since this run started (%s); a run ' ...
        'cannot continue on another bvp4c: start a new run_dir.'], t);
end
if ~isequal({E0.library_files.sha256}, {cur.library_files.sha256})
    error([eid ':E_LIBRARY'], 'resume refused: a library file''s SHA-256 differs.');
end
[bad, item] = import_diff(E0, cur, S);
if bad
    error([eid ':E_IMPORT'], 'resume refused: %s.', item);
end
[bad, item] = focus_diff(E0, cur);
if bad
    error([eid ':E_FOCUS'], 'resume refused: %s.', item);
end
override = '';
if ~isequal(E0.code_hashes, cur.code_hashes)
    if ischar(S.resume_code_override)
        override = S.resume_code_override;
        sa_journal('note', run_dir, sprintf('override code-hash change; reason: %s; stored: %s; current: %s', ...
            override, hash_text(E0.code_hashes), hash_text(cur.code_hashes)));
    else
        error([eid ':E_CODE'], ['resume refused: code hashes differ from the stored run ' ...
            '(a repair is a new run; resume_code_override needs a written reason).']);
    end
end
% ---- 2. reconcile ----------------------------------------------------------
J = sa_journal('read', run_dir);
files = attempt_files(run_dir);
ids = sort(cellfun(@(n) sscanf(n, 'A%d.mat'), {files.name}));
started = cellfun(@(s) sscanf(s, 'A%d'), {J.starts.id});
ended = cellfun(@(s) sscanf(s, 'A%d'), {J.ends.id});
if ~isempty(setdiff(ids, started))
    error([eid ':E_CORRUPT'], 'resume refused: attempt file(s) without a journal start: %s.', ...
        mat2str(setdiff(ids, started)));
end
interrupted = [];
for k = 1:numel(started)
    id = started(k);
    if any(ended == id), continue; end
    name = sprintf('A%06d', id);
    af = fullfile(run_dir, 'attempts', [name '.mat']);
    if exist(af, 'file') == 2
        try
            Z = load(af, 'attempt');
        catch ME
            error([eid ':E_CORRUPT'], 'resume refused: %s cannot be read (%s).', af, ME.message);
        end
        sa_journal('end', run_dir, name, sa_search('end_text', Z.attempt));   % cut_by_cap suffix
    else
        A = sa_report('blank');
        A.attempt_id = name; A.kind = 'interrupted'; A.outcome = 'interrupted';
        A.unit_key = J.starts(k).text; A.file = fullfile('attempts', [name '.mat']);
        sa_write('mat', af, struct('attempt', A, 'ukey', J.starts(k).text, 'sol', [], ...
            'provenance', cur.provenance));
        sa_journal('end', run_dir, name, 'interrupted');
        interrupted(end+1) = id; %#ok<AGROW>
    end
end
files = attempt_files(run_dir);
names = sort({files.name});
ids = cellfun(@(n) sscanf(n, 'A%d.mat'), names);
if ~isequal(ids, 1:numel(ids))
    error([eid ':E_CORRUPT'], 'resume refused: attempt ids are not contiguous.');
end
for k = 1:numel(names)
    try
        Z = load(fullfile(run_dir, 'attempts', names{k}), 'attempt');
        assert(strcmp(Z.attempt.attempt_id, names{k}(1:end-4)));
    catch ME
        error([eid ':E_CORRUPT'], 'resume refused: %s is corrupt (%s).', names{k}, ME.message);
    end
end
tmp = [dir(fullfile(run_dir, '*.tmp_*')); dir(fullfile(run_dir, 'attempts', '*.tmp_*')); ...
    dir(fullfile(run_dir, 'roots', '*.tmp_*'))];
tmpl = arrayfun(@(t) fullfile(t.folder, t.name), tmp, 'UniformOutput', false);
% ---- 3. replay plan ----------------------------------------------------------
P = struct('active', ~isempty(names), ...
    'files', {cellfun(@(n) fullfile(run_dir, 'attempts', n), names, 'UniformOutput', false)}, ...
    'pos', 1, 'queue_stamp', stamp(fullfile(run_dir, 'queue.mat')), ...
    'attempts_stamp', stamp(fullfile(run_dir, 'attempts.mat')), ...
    'queue_matched', false, 'attempts_matched', false, ...
    'shard_stamp', containers.Map('KeyType', 'char', 'ValueType', 'double'), ...
    'shard_matched', containers.Map('KeyType', 'char', 'ValueType', 'logical'), ...
    'override', override, 'interrupted', interrupted, 'tmp', {tmpl});
% the stored import plan: loaded, never rebuilt on resume
if isfield(E0, 'import'), P.import = E0.import; else, P.import = sa_import('blank'); end
sh = sa_store_io('shard_files', run_dir);   % without the stray temporary files listed above
for k = 1:numel(sh)
    Z = load(fullfile(sh(k).folder, sh(k).name), 'key', 'attempt_count');
    P.shard_stamp(Z.key) = Z.attempt_count;
end
if P.queue_stamp > numel(names) || P.attempts_stamp > numel(names) || ...
        any(cell2mat(P.shard_stamp.values) > numel(names))
    error([eid ':E_CORRUPT'], 'resume refused: a checkpoint is stamped beyond the attempt files.');
end
% ---- 4. resume_<n>.mat -------------------------------------------------------
n = numel(dir(fullfile(run_dir, 'resume_*.mat'))) + 1;
sa_write('mat', fullfile(run_dir, sprintf('resume_%d.mat', n)), struct( ...
    'resume_number', n, 'budget', S.budget, 'checkpoint', S.checkpoint, ...
    'time', datestr(now, 'yyyy-mm-ddTHH:MM:SS'), 'pid', feature('getpid'), ... %#ok<DATST,TNOW1>
    'attempt_count', numel(names), 'override', override, 'interrupted', interrupted, ...
    'stray_tmp', {tmpl}, 'provenance', cur.provenance, 'run_dir', top, ...
    'run_dir_stored', E0.settings.run_dir, 'import_source_dir', import_dir(S), ...
    'budget_stored', E0.settings.budget, 'checkpoint_stored', E0.settings.checkpoint));
sa_journal('note', run_dir, sprintf(['resume %d: %d attempt file(s), %d interrupted, %d stray ' ...
    'temporary file(s), override %d'], n, numel(names), numel(interrupted), numel(tmpl), ~isempty(override)));
end

% =========================================================================
% status (zero writes)

function cur = current(M, parameters, manifest, S)
% S (optional): the resolved settings; their seed_modes select
% the library files (a mode-n file enters the list only when a shape needs it)
theta = M.normalize(parameters);
if nargin >= 4 && ~isempty(S) && isfield(S, 'mode') && strcmp(S.mode, 'focus')
    cur = current_focus(M, theta, parameters, S);    % focus: no library
    return
end
if nargin >= 4 && ~isempty(S)
    [M, Sa, A] = sa_seed_modes('activate', M, S, manifest);
    lib = sa_seed_modes('load', M, manifest, theta, parameters, A);
    seed_modes = Sa.seed_modes;                       % as seeded_atlas_run stores them
else
    lib = M.library('load', manifest, theta, parameters);
    seed_modes = {};
end
cur = struct('M', M, 'theta', theta, 'provenance', sa_settings('provenance'), ...
    'code_hashes', code_hashes(M), 'library_files', lib.files, 'model_version', M.version, ...
    'parameters', parameters);
if ~isempty(seed_modes), cur.seed_modes = seed_modes; end
if nargin >= 4 && ~isempty(S)
    cur.import_check = true;          % the import stamp: computed by status / prepare when needed
end
end

function cur = current_focus(M, theta, parameters, S)
% A focus run reads no library: library_files is empty; the
% source's stored field scale and sheet (family label) lookup are derived (zero solves)
% for the E_FOCUS comparison, and seed_modes is the source's (as the run
% stores it). An unreadable source is recorded (focus_error), not thrown:
% status reports it.
cur = struct('M', M, 'theta', theta, 'provenance', sa_settings('provenance'), ...
    'code_hashes', code_hashes(M), 'library_files', struct('type', {}, 'path', {}, 'bytes', {}, 'sha256', {}), ...
    'model_version', M.version, 'parameters', parameters);
try
    [FI, ~, Sf] = sa_import('focus_inputs', S, M);
    cur.seed_modes = Sf.seed_modes;
    cur.focus_inputs = FI;
catch ME
    cur.focus_error = sprintf('%s: %s', ME.identifier, ME.message);
end
cur.import_check = true;
end

function [tf, item] = focus_legacy(E0, cur)
% a focus folder without focus_inputs (started with a half-cell library by
% an older version), checked for a focus caller before the other
% preconditions: its library files and seed_modes would otherwise be named
tf = false; item = '';
if ~isfield(cur, 'focus_inputs') && ~isfield(cur, 'focus_error'), return; end
if isfield(E0, 'focus_inputs'), return; end
if ~(isfield(E0, 'settings') && isfield(E0.settings, 'mode') && strcmp(E0.settings.mode, 'focus')), return; end
tf = true;
item = ['focus: this folder was started with an uncoupled library (by an earlier build of this software); a focus ' ...
    'run now runs without one'];
end

function S = focus_seed_modes(S, E0, cur)
% a focus caller whose source cannot be read has not derived the source's
% seed_modes: compare with the stored ones (E_FOCUS names the source)
if isfield(cur, 'focus_error') && isfield(E0, 'settings') && isfield(E0.settings, 'seed_modes')
    S.seed_modes = E0.settings.seed_modes;
end
end

function [bad, item] = focus_diff(E0, cur)
% E_FOCUS: a focus run's stored focus_inputs against the
% current derivation from the source (cur.focus_inputs), or the reason it
% failed (cur.focus_error). Checked only when the caller carries either.
bad = false; item = '';
if ~isfield(cur, 'focus_inputs') && ~isfield(cur, 'focus_error'), return; end
if isfield(cur, 'focus_error')
    bad = true;
    item = sprintf('focus: the source''s stored field scale and family lookup cannot be read (%s)', cur.focus_error);
    return
end
if ~isfield(E0, 'focus_inputs')
    [bad, item] = focus_legacy(E0, cur);         % not a focus folder: E_SETTINGS names the mode
    return
end
[bad, why] = sa_import('focus_compare', E0.focus_inputs, cur.focus_inputs);
if bad
    item = sprintf(['focus: the field scale and family lookup derived from the source differ from the stored ' ...
        'ones (%s)'], strjoin(why, ', '));
end
end

function Z = status(top, S, cur)
Z = struct('state', '', 'id', '', 'message', '', 'items', {{}}, 'compacted', false, ...
    'layout', sa_layout('layout', top), 'state_dir', sa_layout('state', top), ...
    'attempt_count', 0, 'stop_reason', '', 'settings', S);
if isfield(cur, 'M') && isfield(cur, 'theta')
    S = sa_lattice('derive_step', cur.M, cur.theta, S);
    Z.settings = S;
end
if isfield(cur, 'seed_modes'), S.seed_modes = cur.seed_modes; Z.settings = S; end   % 'all' resolved
% 1. new
if exist(top, 'dir') ~= 7 || numel(dir(top)) <= 2
    Z.state = 'new'; Z.message = sprintf('%s is new: the run starts there.', top); return
end
% 2. legacy
if strcmp(Z.layout, 'flat')
    Z.state = 'legacy';
    Z.message = sprintf(['%s was written by an early version of this package (flat layout); this software cannot resume it; ' ...
        'choose a new run_dir. It still loads: seeded_atlas_plot / seeded_atlas_trace with opts.out.'], top);
    return
end
D = Z.state_dir;
f = fullfile(D, 'settings_effective.mat');
if exist(f, 'file') ~= 2
    Z.state = 'refused'; Z.id = 'SA:RESUME:E_NO_RUN';
    Z.message = sprintf('%s exists, is not empty and holds no run: choose a new run_dir.', top);
    return
end
% 3. compacted
Z.compacted = exist(fullfile(D, 'compaction.mat'), 'file') == 2;
% 4. preconditions, in prepare's order (E_PARAMETERS first)
E0 = load(f);
S = focus_seed_modes(S, E0, cur);
items = {}; ids = {};
[legacy, item] = focus_legacy(E0, cur);
if legacy                                        % a focus folder started with a library: that reason alone
    Z.state = 'refused'; Z.id = 'SA:RESUME:E_FOCUS'; Z.items = {item};
    Z.message = sprintf('%s differs from the stored run: %s; choose a new run_dir.', top, item);
    return
end
if isfield(cur, 'parameters') && isfield(E0, 'parameters') && ~isequaln(E0.parameters, cur.parameters)
    ids{end+1} = 'E_PARAMETERS';
    items{end+1} = ['parameters: ' strjoin(param_diff(E0.parameters, cur.parameters), ', ')];
end
d = scientific_diff(sa_settings('fill_focus_leaves', E0.settings), sa_settings('fill_focus_leaves', S));
if ~isempty(d), ids{end+1} = 'E_SETTINGS'; items{end+1} = ['settings: ' strjoin(d, ', ')]; end
if isfield(cur, 'provenance') && (~isfield(E0, 'provenance') || ~isequal(E0.provenance, cur.provenance))
    ids{end+1} = 'E_VERSION'; items{end+1} = 'software version, MATLAB version, architecture or thread mode';
end
if isfield(cur, 'model_version') && ~strcmp(E0.model.version, cur.model_version)
    t = toolbox_change(E0.model.version, cur.model_version);
    ids{end+1} = 'E_MODEL';
    if isempty(t)
        items{end+1} = 'model version';
    else
        items{end+1} = sprintf('model version: the MATLAB toolbox changed since this run started (%s)', t);
    end
end
if isfield(cur, 'library_files') && ~isequal({E0.library_files.sha256}, {cur.library_files.sha256})
    ids{end+1} = 'E_LIBRARY'; items{end+1} = 'a library file''s SHA-256';
end
[bad, imp_item] = import_diff(E0, cur, S);
if bad, ids{end+1} = 'E_IMPORT'; items{end+1} = imp_item; end
[bad, foc_item] = focus_diff(E0, cur);
if bad, ids{end+1} = 'E_FOCUS'; items{end+1} = foc_item; end
if isfield(cur, 'code_hashes') && ~isequal(E0.code_hashes, cur.code_hashes) && ...
        ~(isfield(S, 'resume_code_override') && ischar(S.resume_code_override))
    ids{end+1} = 'E_CODE'; items{end+1} = sprintf('code hashes (%s)', strjoin(hash_diff(E0.code_hashes, ...
        cur.code_hashes), ', '));
end
% a changed source refuses only a run that would resume; a completed
% follow-up keeps its state (plot and traces), with a note
imp_note = '';
if isequal(ids, {'E_IMPORT'}) && completed_now(D)
    imp_note = items{1}; ids = {}; items = {};
elseif ~isempty(ids) && all(ismember(ids, {'E_IMPORT', 'E_FOCUS'})) && completed_now(D)
    imp_note = strjoin(items, '; '); ids = {}; items = {};    % focus: a focus source moved or changed
end
if ~isempty(ids)
    Z.state = 'refused'; Z.id = ['SA:RESUME:' ids{1}]; Z.items = items;
    Z.message = sprintf('%s differs from the stored run: %s; choose a new run_dir.', top, strjoin(items, '; '));
    return
end
% 5. integrity note after the last end line
J = sa_journal('read', D);
le = find(startsWith(J.lines, 'end '), 1, 'last'); if isempty(le), le = 0; end
bad = find(startsWith(J.lines, 'note integrity'));
bad = bad(bad > le);
if ~isempty(bad)
    Z.state = 'refused'; Z.id = 'SA:RESUME:E_INTEGRITY'; Z.items = J.lines(bad);
    Z.message = sprintf('%s stopped on an integrity error (%s); a resume would meet it again.', top, ...
        J.lines{bad(end)});
    return
end
% 6-8. completed / budget_exhausted / resumable
nfiles = numel(attempt_files(D));
if Z.compacted
    Zc = load(fullfile(D, 'compaction.mat'), 'attempt_count_before');
    nfiles = Zc.attempt_count_before;
end
Z.attempt_count = nfiles;
qf = fullfile(D, 'queue_final.mat');
if exist(qf, 'file') == 2
    Q = load(qf, 'attempt_count', 'stop_reason');
    Z.stop_reason = Q.stop_reason;
    if Q.attempt_count == nfiles
        if ~strcmp(Q.stop_reason, 'completed') && Z.compacted
            % compacted after a clean stop (seeded_atlas_compact opts.stopped_ok)
            Z.state = 'refused'; Z.id = 'SA:RESUME:E_COMPACTED';
            Z.message = sprintf(['%s was compacted after it stopped (%s) without completing: it can no longer ' ...
                'be resumed (choose a new run_dir for more work; plot, traces and import still work).'], top, ...
                Q.stop_reason);
            return
        end
        if strcmp(Q.stop_reason, 'completed')
            Z.state = 'completed';
            Z.message = sprintf('%s is complete (%d calls): no solve; plot and traces only.', top, nfiles);
            if ~isempty(imp_note), Z.message = sprintf('%s Note: %s.', Z.message, imp_note); end
            if Z.compacted, Z.message = [Z.message ' It is compacted: a resume is refused (E_COMPACTED).']; end
            return
        end
        if strcmp(Q.stop_reason, 'budget_calls') && nfiles >= S.budget.max_calls
            Z.state = 'budget_exhausted';
            Z.message = sprintf(['%s stopped at budget.max_calls with %d calls used; raise budget.max_calls ' ...
                'to continue (no solve now; plot and traces on what exists).'], top, nfiles);
            return
        end
    end
end
if Z.compacted
    Z.state = 'refused'; Z.id = 'SA:RESUME:E_COMPACTED';
    Z.message = sprintf('%s is compacted but not complete: it cannot be resumed.', top);
    return
end
Z.state = 'resumable';
Z.message = sprintf('%s stopped (%s) after %d calls: the run resumes (replay, then live).', top, ...
    tern(isempty(Z.stop_reason), 'killed or no final queue', Z.stop_reason), nfiles);
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

function L = hash_diff(A, B)
fa = {A.file}; fb = {B.file};
L = setxor(fa, fb);
for k = 1:numel(fa)
    j = find(strcmp(fb, fa{k}), 1);
    if ~isempty(j) && ~strcmp(A(k).sha256, B(j).sha256), L{end+1} = fa{k}; end %#ok<AGROW>
end
if isempty(L), L = {'(order)'}; end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function tf = completed_now(D)
% queue_final.mat says completed at the current attempt-file count (the
% 'completed' state below; a compacted folder counts its files before
% compaction)
tf = false;
qf = fullfile(D, 'queue_final.mat');
if exist(qf, 'file') ~= 2, return; end
Q = load(qf, 'attempt_count', 'stop_reason');
n = numel(attempt_files(D));
if exist(fullfile(D, 'compaction.mat'), 'file') == 2
    Zc = load(fullfile(D, 'compaction.mat'), 'attempt_count_before');
    n = Zc.attempt_count_before;
end
tf = strcmp(Q.stop_reason, 'completed') && Q.attempt_count == n;
end

function d = import_dir(S)
d = '';
if isfield(S, 'import') && isstruct(S.import) && isfield(S.import, 'source_dir'), d = S.import.source_dir; end
end

function f = attempt_files(run_dir)
% attempts/A######.mat only. An interrupted save leaves A######.mat.tmp_<pid>
% (R2025b's save keeps the name as given), which A*.mat does not match; the
% name filter guards a save that appends .mat (A######.mat.tmp_<pid>.mat).
% Either way the file is listed (stray temporary files), never read.
f = dir(fullfile(run_dir, 'attempts', 'A*.mat'));
f = f(~cellfun(@isempty, regexp({f.name}, '^A\d{6}\.mat$', 'once')));
end

function c = stamp(f)
c = 0;
if exist(f, 'file') == 2
    Z = load(f, 'attempt_count');
    if isfield(Z, 'attempt_count'), c = Z.attempt_count; end
end
end

function t = hash_text(H)
t = strjoin(arrayfun(@(h) sprintf('%s:%s', h.file, h.sha256(1:12)), H, 'UniformOutput', false), ',');
end
