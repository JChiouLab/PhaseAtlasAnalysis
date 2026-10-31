function result = seeded_atlas_run(parameters, library_manifest, settings, corrector)
%SEEDED_ATLAS_RUN  Standalone seeded coupled atlas: cold, continuation, rescue,
%   restart.
%
%   result = SEEDED_ATLAS_RUN(parameters, library_manifest, settings)
%   result = SEEDED_ATLAS_RUN(..., corrector)
%
%   parameters        the model's dimensional parameters (module 'normalize')
%   library_manifest  where the uncoupled seed library lives (module 'load');
%                     for the original model struct('directory', <outputs dir>)
%   settings          see sa_settings; settings.model and settings.run_dir are
%                     required
%   corrector         OPTIONAL injected corrector with the copied corrector's
%                     signature [sol, info] = f(theta, guess, target, opts).
%                     Tests inject stubs. When omitted the module's real
%                     corrector is called (real solves need numerical authority).
%
%   Path: parameters + library + settings -> theta, lattice, seedability
%   census, declared targets -> the scheduler sa_search (cold pass,
%   refinement, then continuation / rescue passes when routes enable them)
%   -> every call through one executor: validity boundary (sa_validity),
%   per-node root-store shards (copied root-store kernel), journal and atomic
%   attempt files. Run folder (sa_layout): the human-facing files on top
%   (README.txt, summary.txt, roots_index.csv, attempts.csv,
%   settings_effective.txt, pilot_report.txt; plots/ and traces/ for
%   seeded_atlas_plot and seeded_atlas_trace), the state below in
%   internal/:
%     settings_effective.mat  written once: resolved settings, software /
%                             model / MATLAB versions, thread mode, input file
%                             identities (library SHA-256), SHA-256 of every
%                             sa_*.m, the entry point and every module file,
%                             field scale, census, runtime evidence incl. the
%                             runtime.clean flag
%     journal.txt             append-only start / end lines
%     attempts/A######.mat    one per call, atomic, never rewritten
%     roots/<node>.mat        one kernel store per node (rewritten after its
%                             node changes, stamped with the attempt count)
%     queue.mat, attempts.mat derived, rebuildable, every
%                             checkpoint.every_calls and at stop (attempts.csv,
%                             summary.txt, roots_index.csv on top: every
%                             report.files_interval_s seconds and at stop)
%     identity_pairs.csv      at the end
%     resume_<n>.mat          one per resume
%     settings_user.mat       the settings struct exactly as passed (new
%                             runs; sa_settings('from_run') reads it)
%     settings_effective.txt  (top) one line per leaf, value and source
%                             (package / model / mode / user / derived)
%     roots_index.csv         (top) one row per root (sa_report 'roots_index')
%     pilot_report.txt (top), pilot/*.csv   mode pilot: sa_report('pilot') at
%                             every stop (zero solves)
%   settings.run_dir, the stored settings and summary.txt hold the folder
%   the user set (the top); every state reader and writer receives the state
%   folder sa_layout('state', run_dir) (ctx.state_dir).
%   lattice.points_per_axis: the per-axis step is derived (sa_lattice
%   'derive_step') right after theta and before the folder is prepared, so a
%   refused derivation leaves no folder.
%   Modes (sa_modes): smoke | pilot | production; layered defaults
%   package -> model (M.setting_defaults) -> mode -> user (sa_settings).
%   settings.resume = true continues an existing run folder (sa_resume):
%   preconditions, reconciliation, zero-solve replay of every stored attempt
%   with re-derived verdicts and checkpoint comparison, then live.
%   Without resume an existing non-empty run_dir is refused.
%   Screen: stage and status lines through sa_progress (interval
%   report.progress_interval_s; display only).
%   Memory: result.shards and result.store hold the LIGHT
%   in-memory stores (sa_store_io 'light': every discovery's returned_mesh
%   and returned_profile are []); everything else in them is as before. The
%   full store of a node is its shard file (sa_store_io 'read_shard'), or
%   sa_store_io('full', store, M, S, K, sol_of) from the attempt files;
%   sa_store_io('full_map', result) gives every node's full store.
%   The human-facing files (attempts.csv, summary.txt, roots_index.csv)
%   are written every report.files_interval_s seconds (checked after each
%   live call) and at the end; the checkpoint files every
%   checkpoint.every_calls calls and at the end.

if nargin < 4, corrector = []; end
t_run = tic;
if ~(isstruct(settings) && isscalar(settings) && isfield(settings, 'model'))
    error('SA:RUN', 'settings must be a scalar struct with a model field.');
end
M = sa_model(settings.model);
[S, SRC] = sa_settings('resolve', settings, M);
theta = M.normalize(parameters);                 % writes nothing
[S, SRC] = sa_lattice('derive_step', M, theta, S, SRC);
[S, SRC] = sa_settings('scale_counts', M, theta, S, SRC);      % point counts above the reference L
[M, S, A] = sa_seed_modes('activate', M, S, library_manifest); % seed_modes: active types (zero solves)
prepare_run_dir(S.run_dir, S.resume);
sa_progress('begin', S.report.progress_interval_s);
progress_end = onCleanup(@() sa_progress('end')); %#ok<NASGU>
sa_progress('stage', 'inputs (library, lattice, seedability census, targets)', ...
    struct('calls', 0, 'max_calls', S.budget.max_calls, 'roots', 0));
D = sa_layout('state', S.run_dir);               % the state folder, computed once

% ---- inputs ----------------------------------------------------------------
lib = sa_seed_modes('load', M, library_manifest, theta, parameters, A);
L = sa_modes('cap_lattice', sa_lattice('build', M, theta, S), S);   % lattice cap (Inf = none)
C = sa_seedability(M, lib, theta, L);
T = sa_targets(L, C, S);
sa_seed_modes('print', M, C, T, S);              % the pre-run count, before the first call
types = unique([{T.list.left_type}, {T.list.right_type}]);
lib = M.library('profiles', lib, types);
fs = M.field_scale(theta, lib, S);
K = sa_identity('store_settings', M, theta, S, fs);

prov = sa_settings('provenance');
runtime = runtime_evidence(M, corrector);
clean_line = clean_text(runtime);
fprintf('%s\n', clean_line);                     % printed before the search starts
hashes = sa_resume('code_hashes', M);

if S.resume
    cur = struct('provenance', prov, 'code_hashes', hashes, 'library_files', lib.files, ...
        'model_version', M.version, 'parameters', parameters);
    replay = sa_resume('prepare', S.run_dir, S, cur);
    sa_progress('stage', sprintf('resume: %d stored call(s) are replayed and checked first (no solves)', ...
        numel(replay.files)), struct('calls', 0, 'max_calls', S.budget.max_calls, 'roots', 0));
else
    effective = struct();
    effective.settings = S;
    effective.settings_sources = SRC;
    effective.parameters = parameters;
    effective.theta = theta;
    effective.model = struct('id', M.id, 'version', M.version, 'copies', M.copies, ...
        'n_states', M.n_states, 'fields', M.fields, 'identity_fields', M.identity_fields, ...
        'pools', M.pools, 'library_types', {M.library_types}, ...
        'required_checks', {[sa_validity('shared_required'), M.required_checks]});
    effective.library_manifest = library_manifest;
    effective.library_files = lib.files;
    effective.field_scale = fs;
    effective.store_settings = K;
    effective.lattice = rmfield(L, {'half_targets'});
    effective.census = struct('names', {C.names}, 'counts', C.counts, ...
        'n_in_range', sum(C.in_range), 'n_outside', sum(~C.in_range), ...
        'worst_offset', C.worst_offset, 'rule', C.rule, ...
        'seedable', C.seedable, 'in_range', C.in_range, ...  % for the pilot projection
        'seed_counts', sa_targets('seed_counts', C));        % per_seed cold targets per family
    effective.targets = T;
    effective.matlab_version = version;
    effective.runtime = runtime;
    effective.started_at = datestr(now, 'yyyy-mm-ddTHH:MM:SS'); %#ok<DATST,TNOW1>
    effective.provenance = prov;
    effective.code_hashes = hashes;
    sa_write('mat', fullfile(D, 'settings_effective.mat'), effective);
    sa_write('mat', fullfile(D, 'settings_user.mat'), struct('settings', settings));
    sa_write('text', fullfile(S.run_dir, 'settings_effective.txt'), [strjoin([ ...
        {sprintf('settings_effective (%s; mode %s): leaf, value, source', S.software_version, S.mode)}, ...
        sa_settings('effective_lines', S, SRC)], newline) newline]);
    replay = struct('active', false, 'files', {{}}, 'pos', 1, 'queue_stamp', 0, ...
        'attempts_stamp', 0, 'queue_matched', false, 'attempts_matched', false, ...
        'shard_stamp', containers.Map('KeyType', 'char', 'ValueType', 'double'), ...
        'shard_matched', containers.Map('KeyType', 'char', 'ValueType', 'logical'), ...
        'override', '', 'interrupted', [], 'tmp', {{}});
end

% ---- the scheduler -----------------------------------------------------------
ctx = struct('S', S, 'M', M, 'theta', theta, 'K', K, 'T', T, 'C', C, 'L', L, ...
    'lib', lib, 'corrector', corrector, 'prov', prov, 'field_scale', fs, ...
    'attempts', repmat(sa_report('blank'), 1, 0), 'calls', 0, 'stop_reason', 'running', ...
    'stop', false, 't_run', t_run, 'refinement', sa_report('blank_refinement'), 'pairs', [], ...
    'shards', containers.Map('KeyType', 'char', 'ValueType', 'any'), 'dirty', {{}}, ...
    'Q', [], 'Q_at_last', [], 'ring_calls', 0, 'ring_calls_at_last', 0, ...
    'replay', replay, 'pending_retry', [], 'clean_line', clean_line, ...
    'guess_nodes', zeros(1, 0), 'guess_fp', {cell(1, 0)}, 'state_dir', D, ...
    'nroots', zeros(L.n_nodes, 1), ...                                 % roots per lattice node
    'node_of', containers.Map(L.key(:).', num2cell(1:L.n_nodes)), ...   % key -> node
    't_reports', tic, 'rep', []);                                      % report interval, render cache
ctx = sa_search('run', ctx);
sa_progress('stage', sprintf('final reports (stop: %s)', ctx.stop_reason), ...
    struct('calls', ctx.calls, 'roots', root_count(ctx)));
ctx = sa_search('finish', ctx);
if strcmp(S.mode, 'pilot')
    sa_progress('stage', 'pilot report', struct('calls', ctx.calls, 'roots', root_count(ctx)));
    sa_report('pilot', S.run_dir);                % at every stop reason
end

roots = sa_report('roots', ctx);
if isempty(roots), rs = []; else, rs = [roots.root]; end
result = struct('run_dir', S.run_dir, 'stop_reason', ctx.stop_reason, ...
    'calls', ctx.calls, 'n_targets', numel(T.list), ...
    'n_roots', numel(roots), 'attempts', ctx.attempts, ...
    'refinement', ctx.refinement, 'pairs', ctx.pairs, 'store', struct('roots', rs), ...
    'root_ids', {{roots.root_id}}, 'shards', ctx.shards, 'Q', ctx.Q, ...
    'census', struct('names', {C.names}, 'counts', C.counts, 'n_in_range', sum(C.in_range), ...
    'n_outside', sum(~C.in_range), 'worst_offset', C.worst_offset, 'rule', C.rule), ...
    'targets', T, 'field_scale', fs, 'runtime_clean', runtime.clean, 'L', L, 'settings', S, ...
    'replay_interrupted', replay.interrupted, 'guess_nodes', ctx.guess_nodes, ...
    'wall_s', toc(t_run));
end

% =========================================================================

function n = root_count(ctx)
% The roots held now, for the progress lines (the result's n_roots at the end)
n = sum(ctx.nroots);
end

function prepare_run_dir(d, resume)
if resume
    D = sa_layout('state', d);
    if exist(fullfile(D, 'settings_effective.mat'), 'file') ~= 2
        error('SA:RESUME:E_NO_RUN', 'resume: %s is not a run folder.', d);
    end
    if exist(fullfile(D, 'roots'), 'dir') ~= 7, mkdir(fullfile(D, 'roots')); end
    return
end
if exist(d, 'dir') == 7
    listing = dir(d);
    if numel(listing) > 2
        error('SA:RUN', ['run_dir %s exists and is not empty; a new run needs a new ' ...
            'directory (or settings.resume = true to continue it).'], d);
    end
else
    [ok, msg] = mkdir(d);
    if ~ok, error('SA:RUN', 'Cannot create run_dir %s: %s', d, msg); end
end
sa_layout('create', d);
end

function R = runtime_evidence(M, corrector)
% What actually ran: the MATLAB path, the working folder and the file each
% entry point and corrector resolved to (no campaign path).
R = struct();
R.path_entries = strsplit(path, pathsep);
R.pwd = pwd;
R.which_seeded_atlas_run = which('seeded_atlas_run');
R.corrector_injected = ~isempty(corrector);
if R.corrector_injected
    R.corrector_name = func2str(corrector);
    info = functions(corrector);
    R.corrector_file = '';
    if isfield(info, 'file'), R.corrector_file = info.file; end
else
    R.corrector_name = '<module default>';
    R.corrector_file = '';
end
% every verbatim-copied corrector file of the module (the default corrector)
R.copies_which = struct('file', {}, 'which', {});
for k = 1:numel(M.copies)
    [~, nm] = fileparts(M.copies(k).file);
    R.copies_which(end+1) = struct('file', M.copies(k).file, 'which', which(nm));
end
% Where bvp4c and deval resolve for every file of
% the module folder and for the package root (the private forwarders), and
% sa_deval's own deval (the toolbox)
R.solver_which = solver_which(M);
% runtime.clean = false, with reasons, when pwd, a path entry or a
% resolved package file lies under data/experiments/ or tests/ of this
% repository. Recorded and printed before the search starts; not a refusal.
pkg = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(pkg)));
bad = @(p) startsWith(p, [repo filesep]) && (contains(p, [filesep 'data' filesep 'experiments']) || ...
    contains([p filesep], [filesep 'tests' filesep]));
reasons = {};
if bad(R.pwd), reasons{end+1} = ['pwd ' R.pwd]; end
pe = R.path_entries(cellfun(bad, R.path_entries));
for k = 1:numel(pe), reasons{end+1} = ['path ' pe{k}]; end %#ok<AGROW>
fl = dir(fullfile(pkg, 'sa_*.m'));
names = [{'seeded_atlas_run'}, strrep({fl.name}, '.m', '')];
modf = which(['sa_model_' M.id]);
if ~isempty(modf)
    fl = dir(fullfile(fileparts(modf), '*.m'));
    names = [names, strrep({fl.name}, '.m', '')];
end
for k = 1:numel(names)
    w = which(names{k});
    if ~isempty(w) && bad(w), reasons{end+1} = ['resolved ' w]; end %#ok<AGROW>
end
R.clean = isempty(reasons);
R.clean_reasons = reasons;
end

function W = solver_which(M)
callers = {'sa_validity', 'sa_deval'};
modf = which(['sa_model_' M.id]);
if ~isempty(modf)
    fl = dir(fullfile(fileparts(modf), '*.m'));
    callers = [callers, strrep({fl.name}, '.m', '')];
end
pairs = [repmat({'deval'}, numel(callers), 1), callers(:); repmat({'bvp4c'}, numel(callers) - 1, 1), ...
    callers([1, 3:end]).'];
W = struct('name', {}, 'caller', {}, 'which', {});
for k = 1:size(pairs, 1)
    try
        if isempty(which(pairs{k, 2}))
            w = '<caller not on the path>';
        else
            w = which(pairs{k, 1}, 'in', pairs{k, 2});
        end
    catch ME
        w = ['<which failed: ' ME.message '>'];
    end
    W(end+1) = struct('name', pairs{k, 1}, 'caller', pairs{k, 2}, 'which', w); %#ok<AGROW>
end
end

function t = clean_text(R)
if R.clean
    t = 'runtime.clean: true';
else
    t = sprintf('runtime.clean: false (%d reason(s)): %s', numel(R.clean_reasons), ...
        strjoin(R.clean_reasons, '; '));
end
end
