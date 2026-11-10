function varargout = sa_model(action, varargin)
%SA_MODEL  Resolve and validate a model module: the only model boundary.
%
%   M = SA_MODEL(id)
%       Resolve the module whose factory function is named sa_model_<id>. If
%       the package holds models/<id>/, that folder is put on the path first;
%       otherwise the factory must already be on the path (test stubs). The
%       returned struct is validated against the member list below. No model
%       id, species name or model constant appears in this file.
%
%   SA_MODEL('validate', M)      validate an already-built module struct
%   names = SA_MODEL('members')  the required member list
%
%   REQUIRED MEMBERS
%     id, version               char: module identity and version text
%     unclipped                 logical true: the reaction/ODE applies no
%                               clipping (source-scanned by a test)
%     n_states                  double: rows of the solver state
%     fields                    struct array (name, row): the concentration
%                               rows the positivity rule scans
%     identity_fields           struct array (name, row): root-identity rows
%     pools                     struct array (name): one per conserved pool,
%                               i.e. one per lattice coordinate
%     library_types             cellstr: half types held by the seed library;
%                               families are the ordered pairs 'A_B'
%     domain                    fh(theta) -> [a b] half-domain range
%     normalize                 fh(parameters) -> theta
%     pool_targets              fh(theta, d) -> .values (all half targets),
%                               .labels, .solver_target
%     library                   fh(action, ...): 'load' (manifest, theta, p),
%                               'profiles' (lib, types), 'cell' (lib, theta, d)
%     pair_guess                fh(lib, theta, typeL, cellL, typeR, cellR, mo)
%                               -> [guess, why_unavailable]
%     ode                       fh(x, Y, sol, theta) -> F (the corrector's RHS)
%     solve                     fh(theta, guess, target, solver, mo, corrector)
%                               -> R (.sol .opts .warning_id .warning_msg
%                               .corrector_class .corrector_reason .wall_s
%                               .corrector_nodal_min)
%     convergence               fh(sol, opts, wid, wmsg) -> classifier struct
%     diagnostics               fh(sol, theta, target, mo) -> module checks
%     shifts                    fh(sol) -> .names .values .max
%     shape_label               fh(sol, theta, mo, fs) -> [char, info]: the
%                               family label, a mode name (sa_mode_name). It
%                               sets the attempt's observed_family, the
%                               per-root pass's unit family string (unit
%                               keys, replay identity) and the child's
%                               seed_family; it never chooses a guess, a
%                               node, an order or a cap. fs = field_scale's
%                               struct (anchor_hss_values: the flatness
%                               reference); info for proofs. A
%                               label needing no reference may take
%                               (sol, theta, mo) only
%     field_scale               fh(theta, lib, S) -> .values (column) .source
%     comparison_mesh           fh(theta, S) -> strictly increasing row
%     required_checks           cellstr: module check names that must pass
%     model_options_defaults    struct: defaults of the model-owned settings
%     validate_model_options    fh(model_options) -> validated model_options
%     copies                    struct array: verbatim-copied sources
%
%   OPTIONAL MEMBERS (agreed boundary changes)
%     setting_defaults          scalar struct: the model layer of the shared
%                               settings (sa_settings resolves package ->
%                               model -> mode -> user)
%     shape_label_profile       fh(xi, P, theta, mo, fs) -> [char, info]:
%                               shape_label's rule on a stored profile P
%                               (identity rows on the mesh xi); post tools
%                               only (a plot relabelled)
%     mode_families             struct .left / .right: n x 2 cells {half
%                               mode name, library type}: the library type
%                               a half's mode name corresponds to, per side
%                               (sa_mode_name 'library_family'; the per-root
%                               pass's family lookup)
%     seed_shape                fh(action, ...): the model's seed shapes
%                               (milestone 8; sa_seed_modes): 'resolve'
%                               (shape -> type id, base file, transform),
%                               'shape_of_type', 'order' (the type order),
%                               'mirror' (the drawn name of a right half),
%                               'all' (shapes a library folder offers),
%                               'mode_families' (the map from active types),
%                               'cold_rate', 'defaults'
%     scale_counts              fh(theta, S, SRC) -> [S, SRC]: point counts
%                               that grow with L above
%                               scaling.reference_length (milestone 8)
%     field_scale_settings      cellstr of dotted settings leaves (other than
%                               the parameters) that field_scale reads, e.g.
%                               'lattice.anchor': a focus run takes the
%                               source's stored field scale only when these
%                               leaves equal the source's (sa_import
%                               'focus_inputs'); absent = lattice.anchor and
%                               identity.field_scale_floor
%     reflection                struct('state_perm', p, 'identity_perm', pi)
%                               (milestone 10; sa_import's mirror images):
%                               the left/right half swap as a state row
%                               permutation p (the mirror image of y is
%                               y(p, :)) and the permutation pi it induces
%                               on the identity_fields rows; both are
%                               involutions, validated here. Not part of
%                               the model version.

switch action
    case 'validate', varargout{1} = validate(varargin{:});
    case 'members',  varargout{1} = members();
    otherwise,       varargout{1} = load_module(action);
end
end

% =========================================================================

function names = members()
names = {'id','version','unclipped','n_states','fields','identity_fields', ...
    'pools','library_types','domain','normalize','pool_targets','library', ...
    'pair_guess','ode','solve','convergence','diagnostics','shifts', ...
    'shape_label','field_scale','comparison_mesh','required_checks', ...
    'model_options_defaults','validate_model_options','copies'};
end

function M = load_module(id)
if ~(ischar(id) && isrow(id) && ~isempty(id))
    error('SA:MODEL:E_ID', 'A model id must be nonempty character text.');
end
if isempty(regexp(id, '^[a-z][a-z0-9_]*$', 'once'))
    error('SA:MODEL:E_ID', ...
        ['Model id %s is not a legal identifier (lower case, digits and ' ...
         'underscore, starting with a letter).'], id);
end
fname = ['sa_model_' id];
folder = fullfile(fileparts(mfilename('fullpath')), 'models', id);
if exist(folder, 'dir') == 7 && ~contains([pathsep path pathsep], [pathsep folder pathsep])
    addpath(folder);
end
% The solver folder (traced bvp4c copy, mesh
% trace recorder, sa_deval) beside the model folder. The module factories
% verify the generated copy (sa_bvp4c_setup 'component'); a missing or
% edited copy refuses the run at zero calls (SA:SOLVER:E_TRACED).
solver = fullfile(fileparts(mfilename('fullpath')), 'solver');
if exist(solver, 'dir') == 7 && ~contains([pathsep path pathsep], [pathsep solver pathsep])
    addpath(solver);
end
if isempty(which(fname))
    error('SA:MODEL:E_NOT_FOUND', ...
        ['No model module found for id %s: the factory function %s is not on ' ...
         'the MATLAB path.'], id, fname);
end
M = feval(fname);
M = validate(M);
if ~strcmp(M.id, id)
    error('SA:MODEL:E_ID', ...
        'Module %s reports id %s; the factory name and the id must agree.', fname, M.id);
end
end

function M = validate(M)
id = 'SA:MODEL:E_MEMBER';
if ~(isstruct(M) && isscalar(M))
    error(id, 'A model module must be a scalar struct.');
end
need = members();
have = fieldnames(M);
missing = setdiff(need, have);
if ~isempty(missing)
    error(id, ['Model module is missing required member(s): %s. The boundary ' ...
        'is closed: a module that cannot fill a member is redesigned, not ' ...
        'worked around.'], strjoin(sort(missing)', ', '));
end
extra = setdiff(have, [need, {'setting_defaults', 'shape_label_profile', 'mode_families', ...
    'seed_shape', 'scale_counts', 'reflection', 'field_scale_settings'}]);
if ~isempty(extra)
    error(id, ['Model module declares member(s) outside the boundary: %s. Add ' ...
        'them to sa_model only by an agreed boundary change.'], ...
        strjoin(sort(extra)', ', '));
end

require_text(M.id, 'id', id);
require_text(M.version, 'version', id);
if ~(islogical(M.unclipped) && isscalar(M.unclipped) && M.unclipped)
    error(id, ['Member unclipped must be the logical true: the acceptance ' ...
        'boundary is not defined for a clipped reaction (RC-10).']);
end
if ~(isa(M.n_states,'double') && isscalar(M.n_states) && M.n_states >= 1 && ...
        M.n_states == fix(M.n_states))
    error(id, 'Member n_states must be a positive integer.');
end
require_named_rows(M.fields, 'fields', M.n_states, id);
require_named_rows(M.identity_fields, 'identity_fields', M.n_states, id);
if ~(isstruct(M.pools) && isrow(M.pools) && ~isempty(M.pools) && ...
        all(isfield(M.pools, {'name'})))
    error(id, 'Member pools must be a nonempty row struct array with a name.');
end
require_cellstr(M.library_types, 'library_types', id);
require_cellstr(M.required_checks, 'required_checks', id);
for f = {'domain','normalize','pool_targets','library', ...
         'pair_guess','ode','solve','convergence','diagnostics','shifts', ...
         'shape_label','field_scale','comparison_mesh','validate_model_options'}
    if ~isa(M.(f{1}), 'function_handle')
        error(id, 'Member %s must be a function handle.', f{1});
    end
end
if ~(isstruct(M.model_options_defaults) && isscalar(M.model_options_defaults))
    error(id, 'Member model_options_defaults must be a scalar struct.');
end
if ~(isstruct(M.copies) && (isempty(M.copies) || isrow(M.copies)))
    error(id, 'Member copies must be a row struct array (possibly empty).');
end
if isfield(M, 'setting_defaults') && ~(isstruct(M.setting_defaults) && isscalar(M.setting_defaults))
    error(id, 'Optional member setting_defaults must be a scalar struct.');
end
for f = {'seed_shape', 'scale_counts'}
    if isfield(M, f{1}) && ~isa(M.(f{1}), 'function_handle')
        error(id, 'Optional member %s must be a function handle.', f{1});
    end
end
if isfield(M, 'reflection')
    check_reflection(M, id);
end
if isfield(M, 'field_scale_settings') && ~(iscellstr(M.field_scale_settings) && ...
        all(cellfun(@(c) ~isempty(regexp(c, '^[a-z_]+(\.[a-z_]+)+$', 'once')), M.field_scale_settings)))
    error(id, 'Optional member field_scale_settings must be a cellstr of dotted settings leaves.');
end
if isfield(M, 'shape_label_profile') && ~isa(M.shape_label_profile, 'function_handle')
    error(id, 'Optional member shape_label_profile must be a function handle.');
end
if isfield(M, 'mode_families')
    f = M.mode_families;
    ok = isstruct(f) && isscalar(f) && all(isfield(f, {'left', 'right'}));
    for side = {'left', 'right'}
        ok = ok && iscell(f.(side{1})) && size(f.(side{1}), 2) == 2 && ...
            all(ismember(f.(side{1})(:, 2), M.library_types));
    end
    if ~ok
        error(id, ['Optional member mode_families must be a struct with left and right: n x 2 cells of ' ...
            '{half mode name, library type}.']);
    end
end
end

function require_text(v, name, id)
if ~(ischar(v) && isrow(v) && ~isempty(v))
    error(id, 'Member %s must be nonempty character text.', name);
end
end

function check_reflection(M, id)
% Optional member reflection (milestone 10): state_perm p a permutation row
% of 1..n_states, identity_perm pi one of 1..numel(identity_fields), each an
% involution, and pi the permutation p induces on the identity rows: the
% identity row k maps to identity row pi(k) iff p(row_k) = row_pi(k).
R = M.reflection;
if ~(isstruct(R) && isscalar(R) && all(isfield(R, {'state_perm', 'identity_perm'})))
    error(id, 'Optional member reflection must be a struct with state_perm and identity_perm.');
end
p = R.state_perm; q = R.identity_perm;
ni = numel(M.identity_fields);
if ~is_perm(p, M.n_states) || ~is_perm(q, ni)
    error(id, ['Optional member reflection: state_perm must be a permutation row of 1..%d and identity_perm ' ...
        'one of 1..%d.'], M.n_states, ni);
end
if ~isequal(p(p), 1:M.n_states) || ~isequal(q(q), 1:ni)
    error(id, 'Optional member reflection: state_perm and identity_perm must be involutions.');
end
rows = [M.identity_fields.row];
if ~isequal(p(rows), rows(q))
    error(id, ['Optional member reflection: identity_perm must be the permutation state_perm induces on the ' ...
        'identity_fields rows.']);
end
end

function tf = is_perm(v, n)
tf = isa(v, 'double') && isrow(v) && numel(v) == n && isequal(sort(v), 1:n);
end

function require_cellstr(v, name, id)
if ~iscell(v) || isempty(v) || ~isrow(v) || ...
        ~all(cellfun(@(c) ischar(c) && isrow(c) && ~isempty(c), v))
    error(id, 'Member %s must be a nonempty cell row of nonempty char rows.', name);
end
end

function require_named_rows(v, name, n_states, id)
if ~(isstruct(v) && isrow(v) && ~isempty(v) && all(isfield(v, {'name','row'})))
    error(id, 'Member %s must be a nonempty row struct array with name and row.', name);
end
rows = [v.row];
if ~isa(rows,'double') || any(rows < 1) || any(rows > n_states) || ...
        any(rows ~= fix(rows)) || numel(unique(rows)) ~= numel(rows)
    error(id, ['Member %s must carry distinct integer row indices inside ' ...
        '1..n_states.'], name);
end
for k = 1:numel(v)
    require_text(v(k).name, sprintf('%s(%d).name', name, k), id);
end
end
