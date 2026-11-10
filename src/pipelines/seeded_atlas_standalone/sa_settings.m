function varargout = sa_settings(action, varargin)
%SA_SETTINGS  Documented settings defaults, validation and resolution.
%
%   [S, SRC] = SA_SETTINGS('resolve', user, M)
%       Fill documented defaults, reject unknown or contradictory fields and
%       return the fully resolved ("effective") settings struct. `user` is the
%       caller's settings struct; `M` is the model module (sa_model), which
%       supplies the schema and defaults of the model-owned group
%       `model_options` -- shared settings name no species and no model
%       constant. LAYERED defaults,
%       package -> model (optional module member M.setting_defaults) ->
%       mode (sa_modes('layer', mode)) -> user, each layer filling only the
%       leaves the next does not set; SRC has the shape of S with the source
%       of every leaf ('package' | 'model' | 'mode' | 'user' | 'derived').
%       The model and mode layers must set disjoint leaves (fail closed).
%
%   T = SA_SETTINGS('template', model_id, mode)
%       The editable user struct of a mode: model, mode, run_dir '' (to be
%       set); production and focus add budget.max_calls / max_wall_s = NaN,
%       which the resolver refuses until the user sets them (focus also
%       needs import.source_dir).
%   [U, dropped] = SA_SETTINGS('from_run', run_dir, new_mode)
%       The user-set leaves of an earlier run (settings_user.mat), with every
%       mode-layer leaf (sa_modes('mode_leaves')), pilot.*, mode, run_dir,
%       resume and resume_code_override dropped (listed in `dropped`) and the
%       mode set to new_mode; run_dir is left '' (the pilot-to-production
%       path; a pilot root never seeds production).
%   L = SA_SETTINGS('effective_lines', S, SRC)   one line per leaf: value, source
%   P = SA_SETTINGS('leaves', S)                 every leaf path of a struct
%
%   D = SA_SETTINGS('defaults', M)
%       The documented default struct without any user input.
%
%   SA_SETTINGS('describe')
%       Prints the documented default table (one line per leaf).
%
%   GROUPS AND DOCUMENTED DEFAULTS (shared, model-free)
%     model              (char, REQUIRED)  model-module id, resolved by sa_model
%     mode               'smoke'           'smoke' | 'pilot' | 'production' |
%                                          'focus' (selected roots of an
%                                          earlier run as seeds, with a
%                                          family fence; see FOCUS MODE
%                                          below) (a scientific leaf: a
%                                          resume with another mode is refused)
%     run_dir            (char, REQUIRED)  output directory for this run; must
%                                          not exist yet or be empty
%     lattice.step       0.02              H; coordinates are anchor + H*k;
%                                          scalar (every axis) or one per pool
%     lattice.points_per_axis []           an odd count n >= 3, scalar or one
%                                          per pool: each axis gets its own
%                                          step, derived from the pool's
%                                          admissible range (sa_lattice
%                                          'derive_step'); set lattice.step or
%                                          lattice.points_per_axis, not both;
%                                          not with a moved anchor
%     lattice.anchor     [0 0]             one entry per model pool
%     lattice.margin     1e-4              admissibility margin on half-pool
%                                          targets
%     lattice.max_abs_coordinate Inf       per-axis
%                                          cap on |coordinate| in lattice
%                                          coordinate units (scalar or one per
%                                          pool); Inf = the admissibility
%                                          boundary (the default)
%     targets.policy     'central_nine'    'central_nine' (the node nearest
%                                          d = 0 and its 8 neighbours) |
%                                          'explicit' | 'all_seedable' |
%                                          'per_seed' (one cold call per
%                                          library seed; production mode) |
%                                          'nearest_anchor' | 'pilot_sample'
%                                          (central_nine plus the evenly
%                                          spaced sample of pilot.*) |
%                                          'per_seed_centre' (per_seed's
%                                          targets, cold calls ordered nearest
%                                          d = 0 first; a run option, never a
%                                          default; milestone 8)
%     targets.list       zeros(0,n_pools)  explicit integer lattice indices
%                                          (each selected family, if seedable)
%     targets.per_family struct()          'nearest_anchor': count per family,
%                                          nodes ordered by index distance to
%                                          the anchor, ties by lattice order
%     families           {}                {} = every seedable family; names are
%                                          checked against the library census
%     seed_modes         {'0','1+','2-'}   the LEFT-drawn half shapes seeded on
%                                          both halves ('0', 'n+', 'n-', n >= 1)
%                                          or 'all' (every shape the library
%                                          folder offers); families are their
%                                          ordered pairs (sa_seed_modes; the
%                                          module's seed_shape maps shapes to
%                                          files); milestone 8
%     scaling.reference_length  Inf        model layer (fission yeast 14,
%                                          noxa 6): above it the model scales
%                                          its point counts and mesh limits
%                                          with L (module scale_counts;
%                                          source 'derived'); at or below it
%                                          nothing changes
%     solver.reltol      1e-6
%     solver.abstol      1e-8
%     solver.nmax        20000
%     solver.mesh_reject 20000             explicit, NOT the corrector's 3000
%     validity.n_dense                20001  uniform dense grid (positivity,
%                                            residual, store profile)
%     validity.interval_samples       [0.25 0.5 0.75]  fractions inside every
%                                            returned mesh interval
%     validity.residual_flag_factor   10     flag scaled sup > factor*reltol
%     validity.refine_enabled         true   refinement re-solves after
%                                            the cold pass (budget permitting)
%     validity.refine_reltol          1e-8   refinement re-solve
%     validity.refine_abstol          1e-10
%     validity.refine_expected_distance 1e-4
%     identity.similarity_tolerance   1e-3   provisional (not calibrated)
%     identity.near_threshold         [1e-4 1e-2]
%     identity.metric_id              'scaled_rms_v1'
%     identity.n_comparison_nodes     2001   comparison mesh node count
%     identity.field_scale_floor      1e-3
%     identity.relsup_nodes           4001   dedup_roots.py common grid for the
%                                            relative-sup report distance
%     budget.max_calls   35                 corrector calls, refinement included
%     budget.max_wall_s  28800              checked before every call
%     budget.wall_s      Inf                clean per-process stop
%                                           (budget_wall) under the outer alarm
%     budget.ring_calls  Inf                cap on the
%                                           ring phase's calls (stop budget_ring)
%   ROUTE GROUPS (PACKAGE defaults are cold only; the mode layer of every
%   mode turns continuation and rescue
%   on, see sa_modes)
%     routes.cold / .continuation / .rescue   true / false / false
%     continuation.step_rule        'adaptive'  'adaptive' | 'myxo_restart'
%     continuation.max_halvings     4           h; adaptive edges refused when
%                                               refusals > h (<= 2h+1 calls)
%     continuation.starts           'anchor_and_cold'  | 'anchor'
%     continuation.max_edges_per_node 2         per lineage and node
%     continuation.max_lineages     10000       non-anchor lineage slots
%                                               (milestone 10; was 16)
%     continuation.max_calls_per_lineage 10000  non-anchor lineages
%                                               (milestone 10; was 100)
%     continuation.box              []          [] = whole lattice, or [b1 b2]:
%                                               nodes with |k_c| <= b_c
%     continuation.region           'box'       'box' (continuation.box) |
%                                               'cross' (the lattice lines
%                                               through the central node; box
%                                               must be [])
%     continuation.guess_mesh       'parent'    (sa_guess_mesh) mesh
%                                               of every guess built from a
%                                               stored iterate; 'parent' =
%                                               its own mesh | 'uniform' |
%                                               'subsample'
%     continuation.guess_nodes      1000        N for 'uniform' / 'subsample'
%     continuation.families         {}          the family fence: {} = none;
%                                               a cell row of observed family
%                                               labels (mode names such as
%                                               '1+/0'); or 'selected' (with
%                                               an import: the labels of the
%                                               selected source roots, plus
%                                               their mirror labels when
%                                               import.mirror is on; the
%                                               focus default; refused,
%                                               SA:IMPORT:E_FENCE, when it
%                                               resolves to no label). A root whose
%                                               observed family (its
%                                               representative call's label)
%                                               is outside the fence is
%                                               stored, reported and counted,
%                                               but never made expandable by
%                                               a lineage, never starts a
%                                               lineage (the anchor lineage
%                                               included), never starts the
%                                               per-root pass and is never a
%                                               rescue parent or secant
%                                               basis. Steps aimed at fenced
%                                               families may land on other
%                                               families; those roots are
%                                               kept. Not the top-level
%                                               families (library families
%                                               for cold seeding)
%     rescue.rungs                  {'diagonal','secant','deep'} (that order)
%     rescue.max_calls_per_node     16
%     rescue.ring                   true        rescue
%                                               also in the ring band; false =
%                                               continuation only outside the
%                                               library range
%     pilot.sample_per_axis         7           pilot_sample nodes per axis
%                                               (user-set only in pilot mode)
%     report.large_change_factor    3           edge_large_change report-only factor
%                                               (may change on resume)
%     report.replay_s_per_record    0.65        replay seconds per
%                                               stored record on a resume, used
%                                               only by the pilot report (live
%                                               overhead after a resume; replay
%                                               term of the projection); model
%                                               layer values in the model
%                                               factories (report-only)
%     report.solver_warnings        false       true prints bvp4c's warnings
%                                               on the screen again (debugging);
%                                               every call records its warning
%                                               either way (screen only)
%     report.progress_interval_s    300         seconds between progress
%                                               status lines (Inf = stage lines
%                                               only; screen only)
%     report.files_interval_s       600         seconds between rewrites of
%                                               attempts.csv, summary.txt and
%                                               roots_index.csv (checked after
%                                               each live call; always at the
%                                               end; Inf = at the end only;
%                                               the checkpoint files keep
%                                               checkpoint.every_calls;
%                                               report-only, may change on
%                                               resume)
%     search.max_passes             3           continuation+rescue passes
%     search.ring_last              true        ring nodes (outside the
%                                               library's native range) only
%                                               after the interior
%     multiseed.max_roots_per_node  10000       no lineage edge into a node
%                                               holding this many distinct
%                                               roots (milestone 10; was 4)
%     checkpoint.every_calls        25          checkpoint files (attempts.mat,
%                                               queue.mat, process_time.mat)
%                                               rewritten
%     checkpoint.min_free_gib       8           stop disk_floor below it
%                                               (checked at start and at
%                                               every checkpoint)
%     resume                        false       true: continue run_dir
%     resume_code_override          false       false, or the reason (text)
%                                               for a lead-approved resume
%                                               across a code-hash change
%     identity.candidate_mesh       'dense'     | 'comparison'
%   MESH CAP, PASS AND THINNING GROUPS (PACKAGE defaults are neutral; the
%   defaults in force sit in the model layer --
%   cold.mesh_cap, root_pass.enabled -- and the mode layer -- rescue.rungs,
%   rescue.ring; see sa_modes and the model factories)
%     cold.mesh_cap                 Inf         cold calls only: NMax =
%                                               mesh_reject = cap (positive
%                                               integer or Inf = off); model
%                                               layer 10000 (both models);
%                                               must not exceed
%                                               solver.mesh_reject
%     root_pass.enabled             false       per-root continuation pass
%                                               (sa_root_pass); model layer:
%                                               fission yeast true, noxa false
%     root_pass.max_halvings        1           h of a pass edge (<= 2h+1
%                                               calls); <= continuation.max_halvings
%     root_pass.max_rounds          10000       rounds per search pass
%                                               (milestone 10; was 8)
%     root_pass.max_calls_per_round Inf         reached: the pass ends for the
%                                               current search pass
%                                               (status root_pass_cap)
%     targets.thin.families         {}          cold-target thinning lever
%                                               (never a default): families
%                                               whose (family, node) pairs are
%                                               kept only on a stride
%     targets.thin.stride           1           positive integer
%     targets.thin.offset           []          one integer per pool in
%                                               [0, stride-1] ([] = zeros)
%   IMPORT GROUP (milestone 10; inert while import.source_dir is ''):
%     import.source_dir             ''          a finished run folder whose
%                                               roots seed this run: each is
%                                               re-solved once at the nearest
%                                               node (sa_import); a free leaf
%                                               on resume (its identity is
%                                               the source stamp, E_IMPORT)
%     import.guess_mesh             'source'    'source' (the stored x, y on
%                                               its own mesh) | 'rule'
%                                               (sa_guess_mesh 'from')
%     import.max_map_steps          1           a source root farther than
%                                               this (in new-lattice steps)
%                                               from every node is skipped
%     import.mirror                 true        also import each source
%                                               root's mirror image; a model
%                                               without a reflection member
%                                               resolves a package-default
%                                               true to false (printed); a
%                                               user-set true is refused
%     import.source_edges           'none'      'none' | 'same_or_stronger' |
%                                               'all': skip edges the source
%                                               refused (same grid, solver
%                                               settings and model version
%                                               only) | 'auto' (the focus
%                                               default): 'same_or_stronger'
%                                               when those conditions hold,
%                                               else 'none' with a printed
%                                               note
%     import.root_ids               {}          a cell row of source root ids
%                                               ('k:1,4/R000001', as in the
%                                               source's roots_index.csv), or
%                                               the path of a CSV file with a
%                                               root_id column
%     import.families               {}          a cell row of observed family
%                                               labels of source roots
%                                               ('1+/0'); a source root is
%                                               imported when it matches
%                                               either list; both empty =
%                                               every root. An id or label
%                                               the source does not hold is
%                                               refused (SA:IMPORT:E_SELECTION,
%                                               listing them). Mirror rows
%                                               are made for selected roots
%                                               only
%     import.source_stopped         false       true: accept a source that
%                                               was killed (no lock or pid
%                                               file exists, so a running
%                                               source cannot be told from a
%                                               killed one: set it only when
%                                               the source will not write
%                                               again). Its calls without an
%                                               end line are ignored; a
%                                               journal written in the last
%                                               10 minutes is refused; its
%                                               identity is the SHA-256 of
%                                               each planned row's attempt
%                                               file (the source may be
%                                               resumed later); source-edge
%                                               knowledge comes from its
%                                               last checkpoint queue.mat
%   FOCUS MODE: the production mode layer with routes.cold false (import
%   only); import.source_dir is required; budgets are required from the
%   user; import.source_edges defaults to 'auto' and continuation.families
%   to 'selected' (source 'mode'; a value set by the user wins, {} = no
%   fence).
%     model_options      (from M.model_options_defaults; for the original model
%                         hss_branch, guess mesh, corrector gate tolerances and
%                         the family label's constants (shape_flat_tol,
%                         shape_turn_frac, shape_edge_frac, shape_nodes; see
%                         sa_mode_name); main_variable, the '+'-sum of the
%                         model's fields that the family label and the
%                         plots' third coordinate look at, see
%                         sa_main_variable)
%
%   The solver itself is NOT a setting: seeded_atlas_run takes an optional
%   injected corrector (tests use stubs); otherwise the module's real one.
%
%   DERIVED MEMBERS filled by the resolver and rejected in user input:
%     software_version, resolved_at, model_version, settings_format

switch action
    case 'resolve',  [varargout{1}, varargout{2}] = resolve(varargin{:});
    case 'template', varargout{1} = template(varargin{:});
    case 'from_run', [varargout{1}, varargout{2}] = from_run(varargin{:});
    case 'effective_lines', varargout{1} = effective_lines(varargin{:});
    case 'leaves',   varargout{1} = leaves(varargin{:});
    case 'defaults', varargout{1} = defaults(varargin{:});
    case 'describe', describe(varargin{:});
    case 'version',  varargout{1} = software_version();
    case 'provenance', varargout{1} = provenance();
    case 'resume_free', varargout{1} = resume_free_leaves();
    case 'legacy_fill', varargout{1} = legacy_fill(varargin{:});
    case 'fill_focus_leaves', varargout{1} = fill_focus_leaves(varargin{:});
    case 'scale_counts', [varargout{1}, varargout{2}] = scale_counts(varargin{:});
    otherwise
        error('SA:SETTINGS:E_ACTION', 'Unknown sa_settings action %s.', action);
end
end

% =========================================================================

function v = software_version()
v = 'seeded_atlas_standalone/1.0.0-rc2';
end

function P = provenance()
% The versions that verdict re-derivation and the resume preconditions compare;
% arch: no resume across machine architectures.
n = maxNumCompThreads;
if n == 1, tm = 'single'; else, tm = sprintf('multi:%d', n); end
P = struct('software_version', software_version(), 'matlab_version', version, ...
    'thread_mode', tm, 'arch', computer('arch'));
end

function L = resume_free_leaves()
% Leaves that may differ on resume; everything else is a
% scientific leaf and must equal the stored effective settings. A dotted
% entry frees one leaf of a group (import.source_dir: the source folder may
% move; its identity is the source stamp, sa_resume E_IMPORT).
L = {'budget', 'checkpoint', 'resume', 'resume_code_override', 'resolved_at', 'run_dir', ...
    'report', ...                         % report-only, may change on resume
    'import.source_dir'};
end

function v = settings_format()
v = 'sa_settings/1';
end

function D = defaults(M)
if nargin < 1, M = []; end
D = struct();
D.model            = '';
D.mode             = 'smoke';
D.run_dir          = '';
D.lattice  = struct('step', 0.02, 'anchor', [], 'margin', 1e-4, 'max_abs_coordinate', Inf, ...
                    'points_per_axis', []);
D.targets  = struct('policy', 'central_nine', 'list', zeros(0,0), ...
                    'per_family', struct(), ...
                    'thin', struct('families', {{}}, 'stride', 1, 'offset', []));   % lever (off)
D.families = {};
D.seed_modes = {'0', '1+', '2-'};                                      % milestone 8: today's shapes
D.scaling  = struct('reference_length', Inf);                          % off at package level
D.solver   = struct('reltol', 1e-6, 'abstol', 1e-8, 'nmax', 20000, ...
                    'mesh_reject', 20000);
D.validity = struct('n_dense', 20001, ...
                    'interval_samples', [0.25 0.5 0.75], ...
                    'residual_flag_factor', 10, ...
                    'refine_enabled', true, ...
                    'refine_reltol', 1e-8, 'refine_abstol', 1e-10, ...
                    'refine_expected_distance', 1e-4);
D.identity = struct('similarity_tolerance', 1e-3, ...
                    'near_threshold', [1e-4 1e-2], ...
                    'metric_id', 'scaled_rms_v1', ...
                    'n_comparison_nodes', 2001, ...
                    'field_scale_floor', 1e-3, ...
                    'relsup_nodes', 4001, ...
                    'candidate_mesh', 'dense');
D.budget   = struct('max_calls', 35, 'max_wall_s', 28800, 'wall_s', Inf, ...
                    'ring_calls', Inf);
D.routes   = struct('cold', true, 'continuation', false, 'rescue', false);
D.continuation = struct('step_rule', 'adaptive', 'max_halvings', 4, ...
                    'starts', 'anchor_and_cold', 'max_edges_per_node', 2, ...
                    'max_lineages', 10000, 'max_calls_per_lineage', 10000, 'box', [], ...
                    'guess_mesh', 'parent', 'guess_nodes', 1000, 'region', 'box', ...
                    'families', {{}});                                   % the family fence (off)
D.rescue   = struct('rungs', {{'diagonal', 'secant', 'deep'}}, 'max_calls_per_node', 16, ...
                    'ring', true);
D.pilot    = struct('sample_per_axis', 7);
D.cold     = struct('mesh_cap', Inf);                                  % off at package level
D.root_pass = struct('enabled', false, 'max_halvings', 1, 'max_rounds', 10000, ...
                    'max_calls_per_round', Inf);                         % off at package level
D.report   = struct('large_change_factor', 3, 'replay_s_per_record', 0.65, 'solver_warnings', false, ...
                    'progress_interval_s', 300, ...                      % screen only
                    'files_interval_s', 600);                            % human-facing files
D.search   = struct('max_passes', 3, 'ring_last', true);
D.multiseed = struct('max_roots_per_node', 10000);
D.import   = struct('source_dir', '', 'guess_mesh', 'source', 'max_map_steps', 1, ...
                    'mirror', true, 'source_edges', 'none', ...
                    'root_ids', {{}}, 'families', {{}}, 'source_stopped', false);   % (off)
D.checkpoint = struct('every_calls', 25, 'min_free_gib', 8);
D.resume   = false;
D.resume_code_override = false;
if isstruct(M) && isfield(M, 'model_options_defaults')
    D.model_options = M.model_options_defaults;
    D.lattice.anchor = zeros(1, numel(M.pools));
else
    D.model_options = struct();
end
end

function S = legacy_fill(S)
% Fill the leaves missing from an older stored effective settings struct
% with their NEUTRAL values (cap
% Inf, pass off, no thinning, import off), so loaders treat a legacy folder
% as it ran. Never written back; leaves present are left untouched. A
% missing root_pass group takes today's defaults (max_rounds 10000): inert,
% since such a folder (before milestone 5) ran with the pass off.
if ~isstruct(S), return; end
D = defaults();
if isfield(S, 'lattice') && isstruct(S.lattice) && ~isfield(S.lattice, 'points_per_axis')
    S.lattice.points_per_axis = [];
end
if ~isfield(S, 'cold') || ~isstruct(S.cold), S.cold = D.cold; end
if ~isfield(S.cold, 'mesh_cap'), S.cold.mesh_cap = Inf; end
if ~isfield(S, 'seed_modes'), S.seed_modes = D.seed_modes; end         % before 0.12: today's shapes
if ~isfield(S, 'scaling') || ~isstruct(S.scaling), S.scaling = D.scaling; end
if ~isfield(S, 'root_pass') || ~isstruct(S.root_pass), S.root_pass = struct(); end
for f = fieldnames(D.root_pass).'
    if ~isfield(S.root_pass, f{1}), S.root_pass.(f{1}) = D.root_pass.(f{1}); end
end
S.root_pass.enabled = logical(S.root_pass.enabled);
if ~isfield(S, 'import') || ~isstruct(S.import), S.import = D.import; end   % before 1.0.0-rc2: no import
S = fill_focus_leaves(S);
if isfield(S, 'targets') && isstruct(S.targets)
    if ~isfield(S.targets, 'thin') || ~isstruct(S.targets.thin), S.targets.thin = D.targets.thin; end
    for f = fieldnames(D.targets.thin).'
        if ~isfield(S.targets.thin, f{1}), S.targets.thin.(f{1}) = D.targets.thin.(f{1}); end
    end
end
end

function S = fill_focus_leaves(S)
% The selection, stopped-source and fence leaves (focus mode) at their
% neutral values, where a stored struct lacks them (a folder made before
% they existed ran without them); leaves present are left untouched. Used
% by legacy_fill and by the resume comparison (sa_resume), so such a
% folder is compared as it ran.
if ~isstruct(S), return; end
if isfield(S, 'import') && isstruct(S.import)
    if ~isfield(S.import, 'root_ids'), S.import.root_ids = {}; end
    if ~isfield(S.import, 'families'), S.import.families = {}; end
    if ~isfield(S.import, 'source_stopped'), S.import.source_stopped = false; end
end
if isfield(S, 'continuation') && isstruct(S.continuation) && ~isfield(S.continuation, 'families')
    S.continuation.families = {};
end
end

function check_label_list(v, name, id)
% a cell row of distinct nonempty texts ({} allowed)
if ~(iscell(v) && (isempty(v) || isrow(v)) && all(cellfun(@(c) ischar(c) && isrow(c) && ~isempty(c), v)) ...
        && numel(unique(v)) == numel(v))
    error(id, '%s must be a cell row of distinct nonempty texts ({} = none).', name);
end
end

function check_seed_modes(v, id)
% grammar only (zero solves); the library-dependent refusals (a missing or
% 'empty' mode file, a grid or par mismatch) come from the module's loader
sid = 'SA:SEED:E_SHAPE';
if ischar(v)
    if ~strcmp(v, 'all')
        error(sid, 'settings.seed_modes must be a cell row of shapes or ''all''.');
    end
    return
end
if ~(iscell(v) && isrow(v) && ~isempty(v) && all(cellfun(@(c) ischar(c) && isrow(c), v)))
    error(id, 'settings.seed_modes must be a nonempty cell row of shape names (or ''all'').');
end
if any(strcmp(v, 'all'))
    if numel(v) > 1
        error(sid, 'settings.seed_modes: ''all'' cannot be mixed with shape names.');
    end
    return
end
for k = 1:numel(v)
    if isempty(regexp(v{k}, '^(0|[1-9][0-9]*[+-])$', 'once'))
        error(sid, ['settings.seed_modes: unknown shape ''%s'' (a shape is ''0'' or ''n+'' / ' ...
            '''n-'' with n >= 1).'], v{k});
    end
end
if numel(unique(v)) ~= numel(v)
    error(sid, 'settings.seed_modes lists a shape twice.');
end
end

function [S, SRC] = scale_counts(M, theta, S, SRC)
% Point counts that grow with L (packet section 6): the module's scale_counts
% (the arithmetic lives in the model layer) above scaling.reference_length;
% nothing at or below it. User values win (the module scales only leaves
% whose source is package or model). The cold cap may only tighten.
if ~isfield(M, 'scale_counts'), return; end
[S, SRC] = M.scale_counts(theta, S, SRC);
cid = 'SA:SETTINGS:E_CONTRADICTION';
if S.solver.mesh_reject > S.solver.nmax
    error(cid, 'after length scaling solver.mesh_reject (%d) exceeds solver.nmax (%d).', ...
        S.solver.mesh_reject, S.solver.nmax);
end
if isfinite(S.cold.mesh_cap) && S.cold.mesh_cap > S.solver.mesh_reject
    error(cid, 'after length scaling cold.mesh_cap (%d) exceeds solver.mesh_reject (%d).', ...
        S.cold.mesh_cap, S.solver.mesh_reject);
end
end

function derived = derived_names()
derived = {'software_version', 'resolved_at', 'model_version', 'settings_format'};
end

function [S, SRC] = resolve(user, M)
id = 'SA:SETTINGS:E_FIELD';
if nargin < 2 || ~isstruct(M)
    error('SA:SETTINGS:E_MODEL', 'resolve needs the model module from sa_model.');
end
if ~isstruct(user) || ~isscalar(user)
    error(id, 'settings must be a scalar struct.');
end
D = defaults(M);

% ---- derived members may never be supplied -----------------------------
for k = derived_names()
    if isfield(user, k{1})
        error('SA:SETTINGS:E_DERIVED', ...
            ['settings.%s is derived by the pipeline and must not be supplied; ' ...
             'it is written into the effective settings.'], k{1});
    end
end

% ---- layers: package -> model -> mode -> user ---------------------------
mode = D.mode;
if isfield(user, 'mode'), mode = user.mode; end
require_text(mode, 'settings.mode', id);
require_member(mode, sa_modes('modes'), 'settings.mode', 'SA:SETTINGS:E_MODE');
ML = struct();
if isfield(M, 'setting_defaults') && ~isempty(M.setting_defaults)
    ML = M.setting_defaults;
end
MO = sa_modes('layer', mode);
lid = 'SA:SETTINGS:E_LAYER';
if ~(isstruct(ML) && isscalar(ML))
    error(lid, 'The model layer (setting_defaults) must be a scalar struct.');
end
bad = intersect(fieldnames(ML), [{'model', 'mode', 'run_dir', 'resume', 'resume_code_override', ...
    'budget', 'checkpoint', 'model_options'}, derived_names()]);
if ~isempty(bad)
    error(lid, 'The model layer may not set %s.', strjoin(bad(:).', ', '));
end
both = intersect(sa_modes('layer_leaves', ML), sa_modes('layer_leaves', MO));
if ~isempty(both)
    error(lid, 'The model layer and the %s mode layer both set %s (they must be disjoint).', ...
        mode, strjoin(both(:).', ', '));
end
SRC = mark_all(D, 'package');

% ---- recursive merge with unknown-field rejection ----------------------
S = merge(D, ML, 'settings(model layer)', id);   SRC = mark(SRC, ML, 'model');
S = merge(S, MO, 'settings(mode layer)', id);    SRC = mark(SRC, MO, 'mode');
S = merge(S, user, 'settings', id);              SRC = mark(SRC, user, 'user');
% focus mode: the defaults that depend on an import (the mode's own, unless
% the user set the leaf)
if strcmp(S.mode, 'focus')
    if strcmp(source_of(SRC, {'import', 'source_edges'}), 'package')
        S.import.source_edges = 'auto'; SRC.import.source_edges = 'mode';
    end
    if strcmp(source_of(SRC, {'continuation', 'families'}), 'package')
        S.continuation.families = 'selected'; SRC.continuation.families = 'mode';
    end
end

% ---- required members ---------------------------------------------------
require_text(S.model, 'settings.model', id);
require_text(S.run_dir, 'settings.run_dir', id);
require_text(S.mode, 'settings.mode', id);
if ~strcmp(S.model, M.id)
    error('SA:SETTINGS:E_MODEL', ...
        'settings.model is %s but the supplied module reports id %s.', S.model, M.id);
end
cid = 'SA:SETTINGS:E_CONTRADICTION';
% contradiction: production and focus need user-set budgets
if any(strcmp(S.mode, {'production', 'focus'})) && (~strcmp(SRC.budget.max_calls, 'user') || ...
        ~strcmp(SRC.budget.max_wall_s, 'user') || any(isnan(double(S.budget.max_calls))) || ...
        any(isnan(double(S.budget.max_wall_s))))
    error(cid, ['mode %s needs budget.max_calls and budget.max_wall_s set by the ' ...
        'user (there is no %s default).'], S.mode, S.mode);
end
% contradiction: focus needs the earlier run whose roots seed it
if strcmp(S.mode, 'focus') && ~(ischar(S.import.source_dir) && ~isempty(S.import.source_dir))
    error(cid, ['mode focus needs import.source_dir: the earlier run folder whose (selected) roots seed ' ...
        'this run (the wrapper''s import_from). Pick roots with import.root_ids or import.families; ' ...
        'without a source use smoke, pilot or production.']);
end
% contradiction: pilot.* leaves outside pilot mode
if isfield(user, 'pilot') && ~strcmp(S.mode, 'pilot')
    error(cid, 'pilot.* settings are given but mode is %s (they belong to pilot mode).', S.mode);
end

% ---- types and ranges ---------------------------------------------------
v = S.lattice.step;
if ~(isa(v, 'double') && isreal(v) && ~isempty(v) && (isscalar(v) || (isrow(v) && ...
        numel(v) == numel(M.pools))) && all(isfinite(v)) && all(v > 0))
    error(id, ['settings.lattice.step must be a finite positive double, scalar or one per ' ...
        'model pool (%d).'], numel(M.pools));
end
check_points_per_axis(S.lattice.points_per_axis, numel(M.pools), id);
require_pos(S.lattice.margin, 'settings.lattice.margin', id);
if ~isa(S.lattice.anchor,'double') || ~isreal(S.lattice.anchor) || ...
        ~isrow(S.lattice.anchor) || numel(S.lattice.anchor) ~= numel(M.pools) || ...
        any(~isfinite(S.lattice.anchor))
    error(id, ['settings.lattice.anchor must be a finite real row with one ' ...
        'entry per model pool (%d).'], numel(M.pools));
end
policies = {'central_nine','explicit','all_seedable','per_seed','nearest_anchor','pilot_sample', ...
    'per_seed_centre'};
v = S.lattice.max_abs_coordinate;
if ~(isa(v, 'double') && isreal(v) && ~isempty(v) && (isscalar(v) || (isrow(v) && ...
        numel(v) == numel(M.pools))) && ~any(isnan(v)) && all(v > 0))
    error(id, ['settings.lattice.max_abs_coordinate must be a positive double (Inf allowed), ' ...
        'scalar or one per model pool (%d).'], numel(M.pools));
end
require_member(S.targets.policy, policies, 'settings.targets.policy', id);
if ~isa(S.targets.list,'double') || ~isreal(S.targets.list) || ...
        (~isempty(S.targets.list) && (size(S.targets.list,2) ~= numel(M.pools) || ...
         any(~isfinite(S.targets.list(:))) || any(S.targets.list(:) ~= fix(S.targets.list(:)))))
    error(id, ['settings.targets.list must be an n-by-%d matrix of integer ' ...
        'lattice indices.'], numel(M.pools));
end
require_scalar_struct(S.targets.per_family, 'settings.targets.per_family', id);
pf = fieldnames(S.targets.per_family);
for k = 1:numel(pf)
    require_posint(S.targets.per_family.(pf{k}), ...
        sprintf('settings.targets.per_family.%s', pf{k}), id);
end
if ~iscell(S.families) || (~isempty(S.families) && ~isrow(S.families))
    error(id, 'settings.families must be a cell row of family names.');
end
check_seed_modes(S.seed_modes, id);
v = S.scaling.reference_length;
if ~(isa(v, 'double') && isscalar(v) && isreal(v) && ~isnan(v) && v > 0)
    error(id, 'settings.scaling.reference_length must be a positive double scalar (Inf = off).');
end
for k = 1:numel(S.families)
    require_text(S.families{k}, sprintf('settings.families{%d}', k), id);
end
require_pos(S.solver.reltol, 'settings.solver.reltol', id);
require_pos(S.solver.abstol, 'settings.solver.abstol', id);
require_posint(S.solver.nmax, 'settings.solver.nmax', id);
require_posint(S.solver.mesh_reject, 'settings.solver.mesh_reject', id);
require_posint(S.validity.n_dense, 'settings.validity.n_dense', id);
if ~isa(S.validity.interval_samples,'double') || ~isrow(S.validity.interval_samples) || ...
        isempty(S.validity.interval_samples) || ...
        any(S.validity.interval_samples <= 0) || any(S.validity.interval_samples >= 1)
    error(id, ['settings.validity.interval_samples must be a nonempty row of ' ...
        'fractions strictly inside (0,1).']);
end
require_pos(S.validity.residual_flag_factor, 'settings.validity.residual_flag_factor', id);
require_logical(S.validity.refine_enabled, 'settings.validity.refine_enabled', id);
require_pos(S.validity.refine_reltol, 'settings.validity.refine_reltol', id);
require_pos(S.validity.refine_abstol, 'settings.validity.refine_abstol', id);
require_pos(S.validity.refine_expected_distance, ...
    'settings.validity.refine_expected_distance', id);
require_pos(S.identity.similarity_tolerance, 'settings.identity.similarity_tolerance', id);
if ~isa(S.identity.near_threshold,'double') || ~isrow(S.identity.near_threshold) || ...
        numel(S.identity.near_threshold) ~= 2 || any(~isfinite(S.identity.near_threshold)) || ...
        any(S.identity.near_threshold <= 0)
    error(id, 'settings.identity.near_threshold must be a finite positive 1x2 row.');
end
require_member(S.identity.metric_id, {'scaled_rms_v1'}, ...
    'settings.identity.metric_id', id);
require_posint(S.identity.n_comparison_nodes, 'settings.identity.n_comparison_nodes', id);
require_pos(S.identity.field_scale_floor, 'settings.identity.field_scale_floor', id);
require_posint(S.identity.relsup_nodes, 'settings.identity.relsup_nodes', id);
require_posint(S.budget.max_calls, 'settings.budget.max_calls', id);
require_pos(S.budget.max_wall_s, 'settings.budget.max_wall_s', id);
require_pos_or_inf(S.budget.wall_s, 'settings.budget.wall_s', id);
require_nonneg_int_or_inf(S.budget.ring_calls, 'settings.budget.ring_calls', id);
require_member(S.identity.candidate_mesh, {'dense', 'comparison'}, ...
    'settings.identity.candidate_mesh', id);
require_scalar_struct(S.routes, 'settings.routes', id);
for f = {'cold', 'continuation', 'rescue'}
    require_logical(S.routes.(f{1}), sprintf('settings.routes.%s', f{1}), id);
end
require_member(S.continuation.step_rule, {'adaptive', 'myxo_restart'}, ...
    'settings.continuation.step_rule', id);
require_posint(S.continuation.max_halvings, 'settings.continuation.max_halvings', id);
require_member(S.continuation.starts, {'anchor', 'anchor_and_cold'}, ...
    'settings.continuation.starts', id);
require_posint(S.continuation.max_edges_per_node, 'settings.continuation.max_edges_per_node', id);
require_nonneg_int_or_inf(S.continuation.max_lineages, 'settings.continuation.max_lineages', id);
require_posint(S.continuation.max_calls_per_lineage, ...
    'settings.continuation.max_calls_per_lineage', id);
require_member(S.continuation.guess_mesh, {'parent', 'uniform', 'subsample'}, ...
    'settings.continuation.guess_mesh', id);
require_member(S.continuation.region, {'box', 'cross'}, 'settings.continuation.region', id);
v = S.continuation.families;
if ischar(v)
    if ~strcmp(v, 'selected')
        error(id, ['settings.continuation.families must be {} (no fence), a cell row of observed family ' ...
            'labels, or ''selected'' (the labels of the imported selection).']);
    end
else
    check_label_list(v, 'settings.continuation.families', id);
end
require_logical(S.rescue.ring, 'settings.rescue.ring', id);
require_posint(S.pilot.sample_per_axis, 'settings.pilot.sample_per_axis', id);
require_pos(S.report.large_change_factor, 'settings.report.large_change_factor', id);
require_pos(S.report.replay_s_per_record, 'settings.report.replay_s_per_record', id);
require_logical(S.report.solver_warnings, 'settings.report.solver_warnings', id);
require_pos_or_inf(S.report.progress_interval_s, 'settings.report.progress_interval_s', id);
require_pos_or_inf(S.report.files_interval_s, 'settings.report.files_interval_s', id);
require_posint(S.continuation.guess_nodes, 'settings.continuation.guess_nodes', id);
if S.continuation.guess_nodes < 3
    error(id, 'settings.continuation.guess_nodes must be at least 3.');
end
b = S.continuation.box;
if ~(isa(b, 'double') && isreal(b) && (isempty(b) || (isrow(b) && numel(b) == numel(M.pools) && ...
        all(isfinite(b)) && all(b >= 0) && all(b == fix(b)))))
    error(id, ['settings.continuation.box must be [] or a row of %d non-negative ' ...
        'integer half-widths (lattice steps).'], numel(M.pools));
end
r = S.rescue.rungs;
allowed = {'diagonal', 'secant', 'deep'};
if ~(iscell(r) && (isempty(r) || isrow(r)) && all(cellfun(@(c) ischar(c) && any(strcmp(c, allowed)), r)) ...
        && numel(unique(r)) == numel(r))
    error(id, 'settings.rescue.rungs must be a cell row of distinct names from: %s.', ...
        strjoin(allowed, ', '));
end
[~, ord] = ismember(r, allowed);
if any(diff(ord) < 0)
    error('SA:SETTINGS:E_CONTRADICTION', ['settings.rescue.rungs must keep the ladder order ' ...
        'diagonal, secant, deep.']);
end
require_posint(S.rescue.max_calls_per_node, 'settings.rescue.max_calls_per_node', id);
% ---- mesh cap, pass and thinning leaves -----------------------------------
require_scalar_struct(S.cold, 'settings.cold', id);
v = S.cold.mesh_cap;
if ~(isa(v, 'double') && isreal(v) && isscalar(v) && ~isnan(v) && v >= 1 && (isinf(v) || v == fix(v)))
    error(id, 'settings.cold.mesh_cap must be a positive integer double scalar or Inf (off).');
end
require_scalar_struct(S.root_pass, 'settings.root_pass', id);
require_logical(S.root_pass.enabled, 'settings.root_pass.enabled', id);
v = S.root_pass.max_halvings;
if ~(isa(v, 'double') && isreal(v) && isscalar(v) && isfinite(v) && v >= 0 && v == fix(v))
    error(id, 'settings.root_pass.max_halvings must be a non-negative integer double scalar.');
end
if v > S.continuation.max_halvings
    error(id, ['settings.root_pass.max_halvings (%d) must not exceed continuation.max_halvings (%d) ' ...
        '(a pass edge makes at most 2h+1 calls).'], v, S.continuation.max_halvings);
end
require_posint(S.root_pass.max_rounds, 'settings.root_pass.max_rounds', id);
v = S.root_pass.max_calls_per_round;
if ~(isa(v, 'double') && isreal(v) && isscalar(v) && ~isnan(v) && v >= 1 && (isinf(v) || v == fix(v)))
    error(id, 'settings.root_pass.max_calls_per_round must be a positive integer double scalar or Inf.');
end
require_scalar_struct(S.targets.thin, 'settings.targets.thin', id);
tf = S.targets.thin.families;
if ~(iscell(tf) && (isempty(tf) || isrow(tf)) && all(cellfun(@(c) ischar(c) && isrow(c) && ~isempty(c), tf)) ...
        && numel(unique(tf)) == numel(tf))
    error(id, 'settings.targets.thin.families must be a cell row of distinct family names.');
end
require_posint(S.targets.thin.stride, 'settings.targets.thin.stride', id);
o = S.targets.thin.offset;
if ~(isa(o, 'double') && isreal(o) && (isempty(o) || (isrow(o) && numel(o) == numel(M.pools) && ...
        all(isfinite(o)) && all(o == fix(o)) && all(o >= 0) && all(o <= S.targets.thin.stride - 1))))
    error(id, ['settings.targets.thin.offset must be [] or a row of %d integers in [0, stride-1] ' ...
        '(one per pool).'], numel(M.pools));
end
require_posint(S.search.max_passes, 'settings.search.max_passes', id);
require_logical(S.search.ring_last, 'settings.search.ring_last', id);
require_posint(S.multiseed.max_roots_per_node, 'settings.multiseed.max_roots_per_node', id);
% ---- import (milestone 10) ------------------------------------------------
require_scalar_struct(S.import, 'settings.import', id);
v = S.import.source_dir;
if ~(ischar(v) && (isempty(v) || isrow(v)))
    error(id, 'settings.import.source_dir must be character text ('''' = no import).');
end
require_member(S.import.guess_mesh, {'source', 'rule'}, 'settings.import.guess_mesh', id);
require_pos(S.import.max_map_steps, 'settings.import.max_map_steps', id);
require_logical(S.import.mirror, 'settings.import.mirror', id);
require_member(S.import.source_edges, {'none', 'same_or_stronger', 'all', 'auto'}, ...
    'settings.import.source_edges', id);
v = S.import.root_ids;
if ischar(v)
    if ~(isrow(v) && ~isempty(v))
        error(id, ['settings.import.root_ids must be a cell row of source root ids or the path of a CSV ' ...
            'file with a root_id column ({} = none).']);
    end
else
    check_label_list(v, 'settings.import.root_ids', id);
end
check_label_list(S.import.families, 'settings.import.families', id);
require_logical(S.import.source_stopped, 'settings.import.source_stopped', id);
require_posint(S.checkpoint.every_calls, 'settings.checkpoint.every_calls', id);
v = S.checkpoint.min_free_gib;
if ~(isa(v, 'double') && isreal(v) && isscalar(v) && isfinite(v) && v >= 0)
    error(id, 'settings.checkpoint.min_free_gib must be a finite non-negative double scalar.');
end
require_logical(S.resume, 'settings.resume', id);
o = S.resume_code_override;
if ~(isequal(o, false) || (ischar(o) && isrow(o) && ~isempty(strtrim(o))))
    error(id, ['settings.resume_code_override must be false or the nonempty reason ' ...
        '(text) for a lead-approved resume across a code-hash change.']);
end

% ---- contradictions ----------------------------------------------------
% contradiction: a production run under a smoke/pilot name
if any(strcmp(S.mode, {'smoke', 'pilot'})) && any(strcmp(S.targets.policy, {'all_seedable', 'per_seed', 'per_seed_centre'})) && ...
        strcmp(S.continuation.region, 'box') && isempty(S.continuation.box)
    error(cid, ['mode %s with targets.policy %s over the whole lattice is a production ' ...
        'run under another name; use mode production.'], S.mode, S.targets.policy);
end
if strcmp(S.continuation.region, 'cross') && ~isempty(S.continuation.box)
    error(cid, 'continuation.region is cross but continuation.box is set (%s).', ...
        mat2str(S.continuation.box));
end
if ~isempty(S.lattice.points_per_axis)
    % by source: the model layer's step never counts as set
    if strcmp(source_of(SRC, {'lattice', 'step'}), 'user')
        error(cid, ['set lattice.step or lattice.points_per_axis, not both (the step is derived ' ...
            'from the count).']);
    end
    if any(S.lattice.anchor ~= 0)
        error(cid, ['lattice.points_per_axis needs the anchor at the equal-halves point (anchor %s): ' ...
            'the derivation is symmetric about it; a moved anchor is a diagnostic, use lattice.step.'], ...
            mat2str(S.lattice.anchor));
    end
end
cap = S.lattice.max_abs_coordinate;
if isscalar(cap), cap = repmat(cap, 1, numel(M.pools)); end
if any(abs(S.lattice.anchor) > cap)
    error(cid, 'lattice.max_abs_coordinate %s excludes the anchor %s.', mat2str(cap), ...
        mat2str(S.lattice.anchor));
end
if strcmp(S.targets.policy, 'explicit') && isempty(S.targets.list)
    error(cid, 'targets.policy is explicit but targets.list is empty.');
end
if ~strcmp(S.targets.policy, 'explicit') && ~isempty(S.targets.list)
    error(cid, 'targets.list is supplied but targets.policy is %s.', S.targets.policy);
end
if strcmp(S.targets.policy, 'nearest_anchor') && isempty(fieldnames(S.targets.per_family))
    error(cid, 'targets.policy is nearest_anchor but targets.per_family is empty.');
end
if ~strcmp(S.targets.policy, 'nearest_anchor') && ~isempty(fieldnames(S.targets.per_family))
    error(cid, 'targets.per_family is supplied but targets.policy is %s.', S.targets.policy);
end
if S.solver.abstol >= S.solver.reltol
    error(cid, 'solver.abstol (%g) must be tighter than solver.reltol (%g).', ...
        S.solver.abstol, S.solver.reltol);
end
if S.validity.refine_abstol >= S.validity.refine_reltol
    error(cid, 'validity.refine_abstol must be tighter than validity.refine_reltol.');
end
if S.validity.refine_reltol >= S.solver.reltol
    error(cid, ['validity.refine_reltol (%g) must be tighter than solver.reltol ' ...
        '(%g): a refinement at the production tolerance measures nothing.'], ...
        S.validity.refine_reltol, S.solver.reltol);
end
if S.solver.mesh_reject > S.solver.nmax
    error(cid, ['solver.mesh_reject (%d) exceeds solver.nmax (%d): bvp4c cannot ' ...
        'return that mesh, so the rejection could never fire.'], ...
        S.solver.mesh_reject, S.solver.nmax);
end
% The cold cap only ever tightens the cold gates (cap <= mesh_reject <=
% nmax, the last by the check above)
if isfinite(S.cold.mesh_cap) && S.cold.mesh_cap > S.solver.mesh_reject
    error(cid, ['cold.mesh_cap (%d) exceeds solver.mesh_reject (%d): the cold cap may only tighten ' ...
        'the cold calls'' NMax and mesh gates.'], S.cold.mesh_cap, S.solver.mesh_reject);
end
if S.root_pass.enabled && ~S.routes.continuation
    error(cid, 'root_pass.enabled is true but routes.continuation is false (the pass is a continuation route).');
end
if ~isempty(S.targets.thin.families) && ~strcmp(S.targets.policy, 'all_seedable')
    error(cid, ['targets.thin.families is set but targets.policy is %s (thinning applies to all_seedable ' ...
        'only; set targets.policy = ''all_seedable'' to thin).'], S.targets.policy);
end
if ~isempty(S.targets.thin.families) && S.targets.thin.stride == 1
    error(cid, 'targets.thin.families is set with targets.thin.stride 1 (nothing would be thinned).');
end
if isempty(S.targets.thin.families) && S.targets.thin.stride > 1
    error(cid, 'targets.thin.stride is %d but targets.thin.families is empty.', S.targets.thin.stride);
end
if S.identity.near_threshold(1) >= S.identity.near_threshold(2)
    error(cid, 'identity.near_threshold must be strictly increasing.');
end
if S.identity.similarity_tolerance < S.identity.near_threshold(1) || ...
        S.identity.similarity_tolerance > S.identity.near_threshold(2)
    error(cid, ['identity.similarity_tolerance (%g) lies outside the ' ...
        'near-threshold band [%g, %g], so no pair could ever be flagged.'], ...
        S.identity.similarity_tolerance, S.identity.near_threshold(1), ...
        S.identity.near_threshold(2));
end
if numel(unique(S.families)) ~= numel(S.families)
    error(cid, 'settings.families lists a family twice.');
end
imp = ~isempty(S.import.source_dir);
if (S.routes.continuation || S.routes.rescue) && ~S.routes.cold && ~imp
    error(cid, ['routes.continuation/rescue need routes.cold or import: continuation starts from ' ...
        'an accepted cold or imported root (the anchor), so a run without either has no start.']);
end
if imp && strcmp(S.mode, 'pilot')
    error(cid, ['import.source_dir is set but mode is pilot: the pilot projection has no model of ' ...
        'imported roots; import in smoke or production mode.']);
end
if ~imp && ~strcmp(S.import.source_edges, 'none')
    error(cid, 'import.source_edges is %s but import.source_dir is empty.', S.import.source_edges);
end
if ~imp && (~isempty(S.import.root_ids) || ~isempty(S.import.families) || S.import.source_stopped)
    error(cid, ['import.root_ids, import.families or import.source_stopped is set but import.source_dir ' ...
        'is empty (they select roots of a source run).']);
end
if ~imp && ischar(S.continuation.families)
    error(cid, ['continuation.families is ''selected'' but import.source_dir is empty (the fence is made ' ...
        'from the labels of the imported selection; give a list of labels instead).']);
end
% import.mirror on a model without a reflection member (sa_model): a
% package-default true resolves to false, with a printed line when an
% import is set; a user-set true is refused
if S.import.mirror && ~isfield(M, 'reflection')
    if strcmp(source_of(SRC, {'import', 'mirror'}), 'user')
        error('SA:IMPORT:E_OPTION', ['import.mirror = true is set but model %s declares no reflection ' ...
            '(module member reflection): mirror images cannot be made.'], M.id);
    end
    S.import.mirror = false;
    SRC.import.mirror = 'derived';
    if imp
        fprintf('import.mirror: model %s declares no reflection; mirror images are off (false).\n', M.id);
    end
end
if ischar(S.resume_code_override) && ~S.resume
    error(cid, 'resume_code_override is given but resume is false.');
end

% ---- model-owned group -------------------------------------------------
if isfield(M, 'validate_model_options') && ~isempty(M.validate_model_options)
    S.model_options = M.validate_model_options(S.model_options);
end

% ---- derived members ---------------------------------------------------
S.software_version = software_version();
S.settings_format  = settings_format();
S.model_version    = M.version;
S.resolved_at      = datestr(now, 'yyyy-mm-ddTHH:MM:SS'); %#ok<DATST,TNOW1>
S = orderfields(S);
for k = derived_names(), SRC.(k{1}) = 'derived'; end
SRC = orderfields(SRC);
end

% =========================================================================
% sources, template, from_run, effective lines

function T = mark_all(D, tag)
T = D;
f = fieldnames(D);
for k = 1:numel(f)
    v = D.(f{k});
    if isstruct(v) && isscalar(v) && ~isempty(fieldnames(v))
        T.(f{k}) = mark_all(v, tag);
    else
        T.(f{k}) = tag;
    end
end
end

function T = mark(T, layer, tag)
if isempty(layer), return; end
f = fieldnames(layer);
for k = 1:numel(f)
    v = layer.(f{k});
    if isstruct(v) && isscalar(v) && ~isempty(fieldnames(v)) && isfield(T, f{k}) && isstruct(T.(f{k}))
        T.(f{k}) = mark(T.(f{k}), v, tag);
    else
        T.(f{k}) = tag;
    end
end
end

function T = template(model, mode)
sa_modes('layer', mode);                     % validates the mode name
T = struct('model', model, 'mode', mode, 'run_dir', '');
if any(strcmp(mode, {'production', 'focus'}))
    T.budget = struct('max_calls', NaN, 'max_wall_s', NaN);
end
end

function [U, dropped] = from_run(run_dir, new_mode)
f = fullfile(sa_layout('state', run_dir), 'settings_user.mat');
if exist(f, 'file') ~= 2
    error('SA:SETTINGS:E_FROM_RUN', '%s holds no settings_user.mat (not a milestone-4 run folder).', run_dir);
end
Z = load(f, 'settings');
U = Z.settings;
sa_modes('layer', new_mode);
drop = union(sa_modes('mode_leaves'), {'mode', 'run_dir', 'resume', 'resume_code_override'});
P = leaves(U);
dropped = P(ismember(P, drop) | startsWith(P, 'pilot.') | strcmp(P, 'pilot'));
for k = 1:numel(dropped)
    U = rm_path(U, strsplit(dropped{k}, '.'));
end
U.mode = new_mode;
U.run_dir = '';
end

function s = rm_path(s, parts)
if numel(parts) == 1
    if isfield(s, parts{1}), s = rmfield(s, parts{1}); end
    return
end
if ~isfield(s, parts{1}) || ~isstruct(s.(parts{1})), return; end
sub = rm_path(s.(parts{1}), parts(2:end));
if isempty(fieldnames(sub)), s = rmfield(s, parts{1}); else, s.(parts{1}) = sub; end
end

function P = leaves(s, path)
if nargin < 2, path = ''; end
P = {};
f = fieldnames(s);
for k = 1:numel(f)
    if isempty(path), p = f{k}; else, p = [path '.' f{k}]; end
    v = s.(f{k});
    if isstruct(v) && isscalar(v) && ~isempty(fieldnames(v))
        P = [P, leaves(v, p)]; %#ok<AGROW>
    else
        P{end+1} = p; %#ok<AGROW>
    end
end
end

function L = effective_lines(S, SRC)
P = leaves(S);
L = cell(1, numel(P));
for k = 1:numel(P)
    parts = strsplit(P{k}, '.');
    v = getfield(S, parts{:}); %#ok<GFLD>
    try
        src = getfield(SRC, parts{:}); %#ok<GFLD>
        if ~ischar(src), src = 'package'; end
    catch
        src = 'package';
    end
    L{k} = sprintf('%-44s %-40s %s', P{k}, value_text(v), src);
end
end

function t = value_text(v)
if ischar(v), t = ['''' v ''''];
elseif islogical(v) || isnumeric(v)
    if isempty(v), t = '[]';
    elseif numel(v) <= 8, t = mat2str(v, 17);
    else, t = sprintf('[%s %d elements]', class(v), numel(v));
    end
elseif iscell(v)
    if all(cellfun(@ischar, v)), t = ['{' strjoin(v, ', ') '}']; else, t = sprintf('{%d cell}', numel(v)); end
elseif isstruct(v), t = 'struct()';
elseif isa(v, 'function_handle'), t = func2str(v);
else, t = class(v);
end
end

% =========================================================================

function out = merge(tmpl, user, path, id)
out = tmpl;
if isempty(user), return; end
fn = fieldnames(user);
known = fieldnames(tmpl);
for k = 1:numel(fn)
    name = fn{k};
    if ~any(strcmp(name, known))
        error('SA:SETTINGS:E_UNKNOWN', ...
            ['unknown field %s.%s. Known fields here: %s. Unknown settings are ' ...
             'rejected rather than silently ignored.'], path, name, ...
            strjoin(known', ', '));
    end
    tv = tmpl.(name);
    uv = user.(name);
    if isstruct(tv) && isscalar(tv) && ~isempty(fieldnames(tv)) && isstruct(uv)
        out.(name) = merge(tv, uv, [path '.' name], id);
    else
        out.(name) = uv;
    end
end
end

function describe(M)
if nargin < 1, M = []; end
D = defaults(M);
fprintf('sa_settings documented defaults (%s)\n', software_version());
walk(D, 'settings');
end

function walk(s, path)
fn = fieldnames(s);
for k = 1:numel(fn)
    v = s.(fn{k});
    p = [path '.' fn{k}];
    if isstruct(v) && isscalar(v) && ~isempty(fieldnames(v))
        walk(v, p);
    else
        fprintf('  %-44s %s\n', p, render(v));
    end
end
end

function t = render(v)
if ischar(v),                t = ['''' v ''''];
elseif islogical(v),         t = mat2str(v);
elseif isa(v,'function_handle'), t = func2str(v);
elseif iscell(v),            t = sprintf('{%d cell}', numel(v));
elseif isstruct(v),          t = 'struct()';
elseif isempty(v),           t = '[]';
else,                        t = mat2str(v);
end
end

function check_points_per_axis(v, np, id)
if isempty(v) && isa(v, 'double'), return; end
if ~(isa(v, 'double') && isreal(v) && (isscalar(v) || (isrow(v) && numel(v) == np)) && ...
        all(isfinite(v)) && all(v == fix(v)))
    error(id, ['settings.lattice.points_per_axis must be [] or integer counts, scalar or one per ' ...
        'model pool (%d).'], np);
end
k = find(v < 3, 1);
if ~isempty(k)
    error(id, 'lattice.points_per_axis = %d is below 3: an axis needs the anchor and one node on each side.', ...
        v(k));
end
k = find(mod(v, 2) == 0, 1);
if ~isempty(k)
    error(id, ['lattice.points_per_axis = %d is even; the anchor is the middle node, so the count must ' ...
        'be odd: use %d or %d. (This is not the uncoupled res, a closed grid over 95 %% of the range ' ...
        'that allows even counts; points_per_axis spans the open physical range.)'], v(k), v(k) - 1, v(k) + 1);
end
end

function t = source_of(SRC, parts)
t = '';
try
    t = getfield(SRC, parts{:}); %#ok<GFLD>
catch
end
if ~ischar(t), t = ''; end
end

function require_text(v, name, id)
if ~(ischar(v) && isrow(v) && ~isempty(v))
    error(id, '%s must be nonempty character text.', name);
end
end

function require_logical(v, name, id)
if ~(islogical(v) && isscalar(v))
    error(id, '%s must be a logical scalar.', name);
end
end

function require_pos_or_inf(v, name, id)
if ~(isa(v,'double') && isreal(v) && isscalar(v) && ~isnan(v) && v > 0)
    error(id, '%s must be a positive double scalar (Inf allowed).', name);
end
end

function require_nonneg_int_or_inf(v, name, id)
if ~(isa(v,'double') && isreal(v) && isscalar(v) && ~isnan(v) && v >= 0 && ...
        (isinf(v) || v == fix(v)))
    error(id, '%s must be a non-negative integer double scalar (Inf allowed).', name);
end
end

function require_pos(v, name, id)
if ~(isa(v,'double') && isreal(v) && isscalar(v) && isfinite(v) && v > 0)
    error(id, '%s must be a finite positive double scalar.', name);
end
end

function require_posint(v, name, id)
if ~(isa(v,'double') && isreal(v) && isscalar(v) && isfinite(v) && v >= 1 && v == fix(v))
    error(id, '%s must be a positive integer double scalar.', name);
end
end

function require_member(v, allowed, name, id)
require_text(v, name, id);
if ~any(strcmp(v, allowed))
    error(id, '%s must be one of: %s.', name, strjoin(allowed, ', '));
end
end

function require_scalar_struct(v, name, id)
if ~(isstruct(v) && isscalar(v))
    error(id, '%s must be a scalar struct.', name);
end
end
