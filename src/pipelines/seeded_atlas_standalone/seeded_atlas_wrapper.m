function W = seeded_atlas_wrapper(cfg)
%SEEDED_ATLAS_WRAPPER  From the uncoupled library to coupled results, in one call.
%
%   W = SEEDED_ATLAS_WRAPPER(cfg) is what run_seeded_atlas.m (repository
%   root) runs after its edit block. cfg fields (see run_seeded_atlas.m and
%   the user guide):
%     model               'original_polarity' | 'guzzo_noxa'
%     par                 the dimensional parameters (both pipelines)
%     res                 the uncoupled seed library grid (noxa: also par.res)
%     library_dir         '' = the uncoupled pipeline's own folder
%                         data/<uncoupled run name>/outputs; built if missing
%     library_confirmation  noxa only: a pre-seeding confirmation file ('')
%     run_dir             '' = data/standalone_runs/<sa_run_name(par, model, mode,
%                         grid)>, the grid being the one the run uses (the
%                         resolved lattice.points_per_axis, else lattice.step,
%                         after points_per_axis, settings_overrides and
%                         from_pilot), e.g. ..._MP3p444_n51_92260295
%     mode                'smoke' | 'pilot' | 'production'; production's cold
%                         targets are per_seed (one cold call per library
%                         seed); settings_overrides.targets.policy =
%                         'all_seedable' calls every seedable node instead
%     points_per_axis     the coupled grid: an odd point count per pool axis
%                         (scalar or one per pool), passed as
%                         settings lattice.points_per_axis; [] = the fixed
%                         lattice.step (the package default, or
%                         settings_overrides.lattice.step). Setting it here and
%                         in settings_overrides.lattice.points_per_axis is
%                         refused (one setting, one place)
%     from_pilot          production: a pilot run folder whose settings to keep
%     budget              struct (max_calls, max_wall_s, ...): a user leaf; it
%                         overrides every mode's budget
%     settings_overrides  struct of further settings leaves (user guide,
%                         tuning table)
%     traces              struct array of trace entries (seeded_atlas_trace)
%     make_plot           true: plots/roots_3d.fig after the run (seeded_atlas_plot:
%                         zero solves, the root files read one at a time, so
%                         memory holds one root file plus the plotted roots,
%                         whatever the run's size)
%     compact_after       false: true runs seeded_atlas_compact after a
%                         completed run (removes the files of calls that found
%                         no root; resume and re-verification of failed calls
%                         are then lost)
%   Test hooks (never in the edit block): corrector (injected corrector),
%   integrator (injected trace integrator), status_only (true: print and
%   return after the status, no build, run, plot or trace), build_opts
%   (passed to atlas_build_libraries), library_manifest (a ready manifest:
%   no library step at all, e.g. the stub model's in-memory library).
%
%   Order. 1. Resolve and PRINT model, parameters, library folder (exists /
%   will be built), run folder and its status, plot path and trace folders;
%   nothing is written before this. 2. Library: the uncoupled folders go on
%   the path, atlas_build_libraries builds or reuses, then exactly the
%   folders the wrapper added are removed again (the uncoupled dependency
%   folder and its models/guzzo_noxa/; the shared noxa core may stay).
%   3. Settings: sa_settings('template', model, mode) (or 'from_run' of
%   from_pilot), then run_dir, budget and the overrides. 4. Status,
%   sa_resume('status') (zero writes), and act on it: new -> run;
%   resumable -> resume = true, run; completed / budget_exhausted -> no
%   solve; legacy / refused -> error, zero calls, the folder untouched.
%   A folder whose status is legacy or refused at step 1 stops there, before
%   the library step: no library is built for a run that would be refused.
%   seeded_atlas_run runs with the repository root as the current folder
%   (the caller's folder is restored afterwards, also on an error), so
%   runtime.clean does not depend on where the wrapper is launched from;
%   when run_dir or the library folder is a relative path, or the package
%   folder is not on the path (found only as the current folder), the
%   current folder is kept and a line says so.
%   5. Plot (when make_plot), then each trace entry (only new or changed
%   ones integrate), then compaction (when compact_after and completed),
%   then the result paths are printed.
%   On a completed folder nothing is replayed or solved: steps 4-5 are the
%   status read, the plot and the traces. The status still checks the
%   resume preconditions, so a completed folder whose stored code hashes
%   differ from the current package's is refused at step 1 (no plot, no
%   trace), unless settings_overrides sets resume = true and
%   resume_code_override = '<reason>' (the override is read, not
%   journaled, when nothing resumes). seeded_atlas_plot(run_dir) and
%   seeded_atlas_trace(run_dir, entry) give the same files without the
%   status check.
%   W: the resolved paths, status (before the run), status_after, result
%   (seeded_atlas_run's, when a run happened), plot, traces, path_at_run,
%   lattice (the grid printed in the header: rule, steps, sizes, n_nodes).
%   Screen: the header gives the grid (count or step, the derived per-axis
%   steps, the node count); stage lines through sa_progress for the library
%   build, the run's stages, the plot, each trace and the compaction.

cfg = fill(cfg);
pkg = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(pkg)));
deps = fullfile(repo, 'src', 'solvers', 'atlas_surface_pipeline_dependencies');
deps_noxa = fullfile(deps, 'models', 'guzzo_noxa');
core = fullfile(repo, 'src', 'models', 'guzzo_noxa');
W = struct('model', cfg.model, 'library_dir', '', 'run_dir', '', 'plot', '', 'trace_dirs', {{}}, ...
    'status', [], 'status_after', [], 'result', [], 'plot_out', [], 'traces', [], 'path_at_run', '', ...
    'compaction', [], 'lattice', []);

% ---- 1. resolve and print ---------------------------------------------------
direct = ~isempty(cfg.library_manifest);
added = {};
if ~direct, added = add_paths(cfg.model, deps, deps_noxa, core); end
cleanup = onCleanup(@() remove_paths(added));
par = cfg.par;
if strcmp(cfg.model, 'guzzo_noxa'), par.res = cfg.res; end
lib_par = cfg.par;
if strcmp(cfg.model, 'guzzo_noxa') && isfield(lib_par, 'res'), lib_par = rmfield(lib_par, 'res'); end
library_dir = cfg.library_dir;
if direct
    manifest = cfg.library_manifest; lib_ok = true;
    if isfield(manifest, 'directory'), library_dir = manifest.directory; end
else
    if isempty(library_dir)
        library_dir = atlas_build_libraries('default_dir', cfg.model, lib_par, cfg.res);
    end
    lib_ok = atlas_build_libraries('present', cfg.model, library_dir);
    manifest = struct('directory', library_dir);
end
if strcmp(cfg.model, 'guzzo_noxa') && ~isempty(cfg.library_confirmation)
    manifest.confirmation = cfg.library_confirmation;
end
M = sa_model(cfg.model);
run_dir = cfg.run_dir;
if isempty(run_dir)
    % the default name shows the grid of the settings the run resolves (the
    % folder is not part of the grid: any text stands in for it here)
    S0 = sa_settings('resolve', user_settings(cfg, fullfile(repo, 'data', 'standalone_runs')), M);
    run_dir = fullfile(repo, 'data', 'standalone_runs', sa_run_name(par, struct('model', cfg.model, ...
        'mode', cfg.mode, 'lattice', S0.lattice)));
end
U = user_settings(cfg, run_dir);
[S, SRCw] = sa_settings('resolve', U, M);
S = sa_settings('scale_counts', M, M.normalize(par), S, SRCw);   % milestone 8 (L above the reference)
sa_progress('begin', S.report.progress_interval_s);
progress_end = onCleanup(@() sa_progress('end')); %#ok<NASGU>
if lib_ok
    cur = sa_resume('current', M, par, manifest, S);
else
    cur = struct('M', M, 'theta', M.normalize(par), 'provenance', sa_settings('provenance'), ...
        'code_hashes', sa_resume('code_hashes', M), 'model_version', M.version, 'parameters', par);
end
Z = sa_resume('status', run_dir, S, cur);
W.lattice = grid_of(M, cur.theta, Z.settings);
W.library_dir = library_dir; W.run_dir = run_dir; W.plot = fullfile(run_dir, 'plots', 'roots_3d.fig');
W.trace_dirs = arrayfun(@(t) fullfile(run_dir, 'traces', t.name), cfg.traces(:).', 'UniformOutput', false);
fprintf('\n================ seeded coupled atlas (%s) ================\n', sa_settings('version'));
fprintf(' model            %s   mode %s\n', cfg.model, cfg.mode);
fprintf(' parameters       %s\n', param_text(par));
fprintf(' library folder   %s   (%s)\n', library_dir, tern(lib_ok, 'exists', 'will be built'));
if isfield(manifest, 'confirmation')
    fprintf(' confirmation     %s\n', manifest.confirmation);
end
fprintf(' grid             %s\n', W.lattice.text);
fprintf(' run folder       %s\n', run_dir);
fprintf(' run status       %s: %s\n', Z.state, Z.message);
fprintf(' plot             %s\n', W.plot);
for k = 1:numel(W.trace_dirs), fprintf(' trace            %s\n', W.trace_dirs{k}); end
fprintf('============================================================\n');
W.status = Z;
if cfg.status_only, return; end
if any(strcmp(Z.state, {'legacy', 'refused'}))   % stop before the library step: nothing is built
    refuse(Z);
end

% ---- 2. library -----------------------------------------------------------------
if ~direct
    if ~lib_ok, sa_progress('stage', 'library build (Stages 1-3 of the uncoupled pipeline)', struct()); end
    manifest = atlas_build_libraries(cfg.model, lib_par, cfg.res, library_dir, cfg.build_opts);
    if strcmp(cfg.model, 'guzzo_noxa') && ~isempty(cfg.library_confirmation)
        manifest.confirmation = cfg.library_confirmation;
    end
end
remove_paths(added); added = {}; %#ok<NASGU>
cleanup = []; %#ok<NASGU>

% ---- 3-4. settings and status ---------------------------------------------------------
cur = sa_resume('current', M, par, manifest, S);
Z = sa_resume('status', run_dir, S, cur);
W.status = Z;
switch Z.state
    case {'new', 'resumable'}
        if strcmp(Z.state, 'resumable'), U.resume = true; end
        fprintf('\n%s: running (%s)\n', run_dir, Z.state);
        W.path_at_run = path;
        back = enter_root(repo, pkg, run_dir, manifest); %#ok<NASGU>
        W.result = seeded_atlas_run(par, manifest, U, cfg.corrector);
        back = []; %#ok<NASGU>                      % the caller's folder again
        fprintf('run: stop %s, %d calls, %d roots\n', W.result.stop_reason, W.result.calls, W.result.n_roots);
    case 'completed'
        fprintf('\n%s is complete: no solve.\n', run_dir);
    case 'budget_exhausted'
        fprintf('\n%s\n', Z.message);
    case {'legacy', 'refused'}
        refuse(Z);
    otherwise
        error('SA:WRAPPER:E_STATE', 'unknown run status %s.', Z.state);
end
W.status_after = sa_resume('status', run_dir, S, cur);
cr = struct();                                   % the run's final counts for the stage lines below
if ~isempty(W.result), cr = struct('calls', W.result.calls, 'roots', W.result.n_roots); end

% ---- 5. plot, traces, compaction -----------------------------------------------------------
if cfg.make_plot
    sa_progress('stage', 'plot (plots/roots_3d.fig)', cr);
    W.plot_out = seeded_atlas_plot(run_dir);
end
T = struct('name', {}, 'dir', {}, 'file', {}, 'ran', {}, 'reason', {}, 'stop', {}, 'frames', {}, ...
    'wall_s', {}, 'frame0_gap', {}, 'at_resolution', {});
for k = 1:numel(cfg.traces)
    e = drop_empty(cfg.traces(k));
    sa_progress('stage', sprintf('trace %d of %d (%s)', k, numel(cfg.traces), entry_name(e)), cr);
    T(end+1) = seeded_atlas_trace(run_dir, e, struct('integrator', cfg.integrator)); %#ok<AGROW>
end
W.traces = T;
if cfg.compact_after
    if strcmp(W.status_after.state, 'completed') && ~W.status_after.compacted
        sa_progress('stage', 'compaction', cr);
        W.compaction = seeded_atlas_compact(run_dir);
    elseif W.status_after.compacted
        fprintf('compact_after: %s is already compacted\n', run_dir);
    else
        fprintf('compact_after: %s is %s, not completed: not compacted\n', run_dir, W.status_after.state);
    end
end
fprintf('\nresults\n  run folder  %s\n  summary     %s\n  roots       %s\n', run_dir, ...
    fullfile(run_dir, 'summary.txt'), fullfile(run_dir, 'roots_index.csv'));
if cfg.make_plot, fprintf('  plot        %s\n', W.plot); end
for k = 1:numel(T), fprintf('  trace       %s\n', T(k).dir); end
end

% =========================================================================

function cfg = fill(cfg)
d = struct('model', 'original_polarity', 'par', struct(), 'res', 50, 'library_dir', '', ...
    'library_confirmation', '', 'run_dir', '', 'mode', 'pilot', 'points_per_axis', [], 'from_pilot', '', ...
    'budget', struct(), ...
    'settings_overrides', struct(), 'traces', struct('name', {}, 'root', {}), 'make_plot', true, ...
    'compact_after', false, 'corrector', [], 'integrator', [], 'status_only', false, 'build_opts', struct(), ...
    'library_manifest', []);
extra = setdiff(fieldnames(cfg), fieldnames(d));
if ~isempty(extra)
    error('SA:WRAPPER:E_FIELD', 'unknown wrapper field(s): %s.', strjoin(extra(:).', ', '));
end
for f = fieldnames(d).'
    if ~isfield(cfg, f{1}), cfg.(f{1}) = d.(f{1}); end
end
if isempty(cfg.traces), cfg.traces = d.traces; end
end

function U = user_settings(cfg, run_dir)
if ~isempty(cfg.from_pilot)
    U = sa_settings('from_run', cfg.from_pilot, cfg.mode);
    U.model = cfg.model;
else
    U = sa_settings('template', cfg.model, cfg.mode);
end
U.run_dir = run_dir;
so = cfg.settings_overrides;
if ~isempty(cfg.points_per_axis)
    if isfield(so, 'lattice') && isstruct(so.lattice) && isfield(so.lattice, 'points_per_axis')
        error('SA:WRAPPER:E_FIELD', ['the grid is set twice: points_per_axis in the edit block and ' ...
            'settings_overrides.lattice.points_per_axis; set it in one place.']);
    end
    U.lattice.points_per_axis = cfg.points_per_axis;
end
U = merge(U, so);
b = cfg.budget;                     % the edit block's budget wins over every other budget leaf
for f = fieldnames(b).'
    U.budget.(f{1}) = b.(f{1});
end
end

function a = merge(a, b)
for f = fieldnames(b).'
    if isstruct(b.(f{1})) && isscalar(b.(f{1})) && isfield(a, f{1}) && isstruct(a.(f{1}))
        a.(f{1}) = merge(a.(f{1}), b.(f{1}));
    else
        a.(f{1}) = b.(f{1});
    end
end
end

function G = grid_of(M, theta, S)
% The coupled grid the run uses (zero solves): rule, per-axis steps, sizes,
% node count, and one line of text for the header.
G = struct('rule', '', 'steps', [], 'sizes', [], 'n_nodes', NaN, 'text', '');
if ~isempty(S.lattice.points_per_axis)
    G.rule = sprintf('points_per_axis %s', mat2str(S.lattice.points_per_axis));
else
    G.rule = sprintf('fixed step %s (points_per_axis [])', mat2str(S.lattice.step));
end
try
    L = sa_modes('cap_lattice', sa_lattice('build', M, theta, S), S);
    G.steps = sa_lattice('steps', L);
    G.sizes = L.index_max - L.index_min + 1;
    G.n_nodes = L.n_nodes;
    G.text = sprintf('%s -> per-axis steps %s, %s = %d nodes', G.rule, mat2str(G.steps, 6), ...
        strjoin(arrayfun(@(n) sprintf('%d', n), G.sizes, 'UniformOutput', false), ' x '), G.n_nodes);
catch ME
    G.text = sprintf('%s (grid not derived: %s)', G.rule, ME.message);
end
end

function n = entry_name(e)
n = '';
if isfield(e, 'name') && ischar(e.name), n = e.name; end
end

function e = drop_empty(e)
for f = fieldnames(e).'
    if isempty(e.(f{1})), e = rmfield(e, f{1}); end
end
end

function added = add_paths(model, deps, deps_noxa, core)
want = {deps};
if strcmp(model, 'guzzo_noxa'), want = {deps, deps_noxa, core}; end
added = {};
for k = 1:numel(want)
    if ~contains([pathsep path pathsep], [pathsep want{k} pathsep])
        addpath(want{k}); added{end+1} = want{k}; %#ok<AGROW>
    end
end
end

function remove_paths(added)
% the uncoupled folders the wrapper added; the shared noxa core stays (the
% coupled noxa module needs it and re-adds it in any case)
for k = 1:numel(added)
    if endsWith(added{k}, fullfile('src', 'models', 'guzzo_noxa')), continue; end
    if contains([pathsep path pathsep], [pathsep added{k} pathsep]), rmpath(added{k}); end
end
end

function t = param_text(p)
f = fieldnames(p);
parts = cellfun(@(n) sprintf('%s=%.17g', n, p.(n)), f, 'UniformOutput', false);
t = strjoin(parts(:).', ' ');
end

function back = enter_root(repo, pkg, run_dir, manifest)
% cd to the repository root for the run; back (onCleanup) returns to the
% caller's folder. Kept in place when a path the run reads is relative, or
% when the package is reachable only through the current folder.
back = [];
p = {run_dir};
for f = {'directory', 'confirmation'}
    if isstruct(manifest) && isfield(manifest, f{1}) && (ischar(manifest.(f{1})) || isStringScalar(manifest.(f{1}))) ...
            && strlength(manifest.(f{1})) > 0
        p{end+1} = char(manifest.(f{1})); %#ok<AGROW>
    end
end
if ~all(cellfun(@(x) startsWith(x, filesep), p))
    fprintf('run folder or library folder is a relative path: the run starts in the current folder %s\n', pwd);
    return
end
if ~any(strcmp(strsplit(path, pathsep), pkg))
    fprintf('the package folder is not on the path: the run starts in the current folder %s\n', pwd);
    return
end
here = pwd;
if strcmp(here, repo), return; end
cd(repo);
back = onCleanup(@() cd(here));
end

function refuse(Z)
if strcmp(Z.state, 'legacy')
    error('SA:WRAPPER:E_LEGACY', '%s', Z.message);
end
error('SA:WRAPPER:E_REFUSED', 'run folder refused (%s; zero calls, the folder is untouched): %s', Z.id, Z.message);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
