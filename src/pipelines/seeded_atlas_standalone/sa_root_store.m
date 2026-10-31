function varargout = sa_root_store(action, varargin)
%SA_ROOT_STORE Validated complete-link identity within one atlas run.
%   Arithmetic and interpolation only; no solver or accounting calls.

switch action
    case 'empty'
        varargout{1} = empty_store(varargin{1});
    case 'add'
        [varargout{1}, varargout{2}] = add_candidate(varargin{1}, varargin{2});
    otherwise
        error('SA:ROOT_STORE:E_CANDIDATE', 'Unknown root-store action.');
end
end

function store = empty_store(settings)
settings = validate_settings(settings, 'SA:ROOT_STORE:E_SETTINGS');
store = struct('schema','seeded_atlas/root-store/2', ...
    'settings',settings,'next_root_number',1, ...
    'next_discovery_number',1,'roots',empty_roots());
end

function [store, outcome] = add_candidate(store, candidate)
validate_store(store);
candidate = validate_candidate(candidate, store.settings);
if ~candidate.valid
    outcome = struct('status','refused_invalid','root_id','', ...
        'member_root_id','','representative_distance',NaN, ...
        'tolerance',store.settings.similarity_tolerance, ...
        'comparisons',empty_comparisons(),'validity',candidate.validity, ...
        'reason','candidate_marked_invalid');
    return
end

profile = interpolate_profile(candidate.mesh, candidate.profile, ...
    store.settings.comparison_mesh, 'SA:ROOT_STORE:E_INTERPOLATION_COVERAGE');
indices = same_target_indices(store.roots, candidate.target_key, candidate.coordinate);
comparisons = comparisons_for(profile, store.roots, indices, store.settings);
eligible = find([comparisons.eligible]);
if isempty(eligible)
    root_id = sprintf('R%06d', store.next_root_number);
    member_index = 0;
    representative_distance = 0;
else
    rows = comparisons(eligible);
    order_key = [[rows.representative_distance].', root_numbers({rows.root_id}).'];
    [~, order] = sortrows(order_key, [1 2]);
    chosen = rows(order(1));
    root_id = chosen.root_id;
    member_index = find(strcmp({store.roots.root_id}, root_id), 1);
    representative_distance = chosen.representative_distance;
end

discovery_number = store.next_discovery_number;
discovery = discovery_record(candidate, profile, discovery_number, root_id, ...
    representative_distance, store.settings.similarity_tolerance, comparisons);
store.next_discovery_number = discovery_number + 1;
if member_index == 0
    root = struct('root_id',root_id, ...
        'created_discovery_number',discovery_number, ...
        'target_key',candidate.target_key,'coordinate',candidate.coordinate, ...
        'comparison_mesh',store.settings.comparison_mesh, ...
        'representative_profile',profile, ...
        'representative_attempt_id',candidate.attempt_id, ...
        'discoveries',discovery);
    store.roots(end+1) = root;
    store.next_root_number = store.next_root_number + 1;
    status = 'new_root';
else
    store.roots(member_index).discoveries(end+1) = discovery;
    status = 'matched_existing';
end
outcome = struct('status',status,'root_id',root_id, ...
    'member_root_id',root_id,'representative_distance',representative_distance, ...
    'tolerance',store.settings.similarity_tolerance, ...
    'comparisons',comparisons,'validity',candidate.validity);
end

function settings = validate_settings(settings, errid)
require_scalar_struct(settings, errid, 'settings');
require_exact_fields(settings, {'comparison_mesh','field_scale', ...
    'similarity_tolerance','metric_id'}, errid, 'settings');
if ~is_real_double_row(settings.comparison_mesh) || ...
        numel(settings.comparison_mesh) < 2 || ...
        any(~isfinite(settings.comparison_mesh)) || ...
        any(diff(settings.comparison_mesh) <= 0)
    error(errid, 'comparison_mesh must be a finite strictly increasing double row.');
end
if ~isa(settings.field_scale,'double') || ~isreal(settings.field_scale) || ...
        ~iscolumn(settings.field_scale) || isempty(settings.field_scale) || ...
        any(~isfinite(settings.field_scale)) || any(settings.field_scale <= 0)
    error(errid, 'field_scale must be a finite positive double column.');
end
if ~isa(settings.similarity_tolerance,'double') || ...
        ~isreal(settings.similarity_tolerance) || ...
        ~isscalar(settings.similarity_tolerance) || ...
        ~isfinite(settings.similarity_tolerance) || settings.similarity_tolerance < 0
    error(errid, 'similarity_tolerance must be a finite nonnegative double scalar.');
end
if ~is_scalar_text(settings.metric_id) || ...
        ~strcmp(char(settings.metric_id),'scaled_rms_v1')
    error(errid, 'metric_id must be scaled_rms_v1.');
end
settings.metric_id = char(settings.metric_id);
end

function candidate = validate_candidate(candidate, settings)
id = 'SA:ROOT_STORE:E_CANDIDATE';
require_scalar_struct(candidate, id, 'candidate');
require_exact_fields(candidate, {'target_key','coordinate','mesh','profile', ...
    'seed_family','observed_family','route','attempt_id','valid','validity'}, ...
    id, 'candidate');
for name = {'target_key','seed_family','observed_family','route','attempt_id'}
    field = name{1};
    if ~is_scalar_text(candidate.(field)) || isempty(char(candidate.(field)))
        error(id, '%s must be nonempty scalar text.', field);
    end
    candidate.(field) = char(candidate.(field));
end
if ~islogical(candidate.valid) || ~isscalar(candidate.valid)
    error(id, 'valid must be a logical scalar.');
end
require_scalar_struct(candidate.validity, id, 'validity');
if ~is_real_double_row(candidate.coordinate) || isempty(candidate.coordinate) || ...
        any(~isfinite(candidate.coordinate))
    error(id, 'coordinate must be a finite nonempty double row.');
end
if ~is_real_double_row(candidate.mesh) || numel(candidate.mesh) < 2 || ...
        any(~isfinite(candidate.mesh)) || any(diff(candidate.mesh) <= 0)
    error(id, 'mesh must be a finite strictly increasing double row.');
end
if ~isa(candidate.profile,'double') || ~isreal(candidate.profile) || ...
        ~ismatrix(candidate.profile) || ...
        size(candidate.profile,1) ~= numel(settings.field_scale) || ...
        size(candidate.profile,2) ~= numel(candidate.mesh) || ...
        any(~isfinite(candidate.profile),'all')
    error(id, 'profile must be a finite real double matrix of declared shape.');
end
if settings.comparison_mesh(1) < candidate.mesh(1) || ...
        settings.comparison_mesh(end) > candidate.mesh(end)
    error('SA:ROOT_STORE:E_INTERPOLATION_COVERAGE', ...
        'candidate mesh does not cover the comparison mesh.');
end
end

function validate_store(store)
schema_id = 'SA:ROOT_STORE:E_STORE_SCHEMA';
require_scalar_struct(store, schema_id, 'store');
if ~isfield(store,'schema') || ~is_scalar_text(store.schema) || ...
        ~strcmp(char(store.schema),'seeded_atlas/root-store/2')
    error(schema_id, 'Unsupported root-store schema.');
end
require_exact_fields(store, {'schema','settings','next_root_number', ...
    'next_discovery_number','roots'}, 'SA:ROOT_STORE:E_STORE_ROOT', 'store');
normalized_settings = validate_settings(store.settings, 'SA:ROOT_STORE:E_STORE_SETTINGS');
if ~isequaln(normalized_settings,store.settings)
    error('SA:ROOT_STORE:E_STORE_SETTINGS', ...
        'Persisted settings must already use canonical character text.');
end
counter_id = 'SA:ROOT_STORE:E_STORE_COUNTER';
if ~is_uint(store.next_root_number) || ~is_uint(store.next_discovery_number)
    error(counter_id, 'Store counters must be positive integer doubles.');
end
if ~isstruct(store.roots) || ~isrow(store.roots)
    error('SA:ROOT_STORE:E_STORE_ROOT', 'roots must be a row struct array.');
end
root_fields = {'root_id','created_discovery_number','target_key','coordinate', ...
    'comparison_mesh','representative_profile','representative_attempt_id', ...
    'discoveries'};
discovery_fields = {'discovery_number','attempt_id','coordinate','seed_family', ...
    'observed_family','route','validity','returned_mesh','returned_profile', ...
    'comparison_profile','member_root_id','representative_distance', ...
    'tolerance','comparisons'};
comparison_fields = {'root_id','representative_distance', ...
    'complete_link_max_distance','eligible'};
require_exact_fields(store.roots,root_fields, ...
    'SA:ROOT_STORE:E_STORE_ROOT','roots');

all_numbers = zeros(1,0);
for r = 1:numel(store.roots)
    root = store.roots(r);
    require_exact_fields(root, root_fields, 'SA:ROOT_STORE:E_STORE_ROOT', 'root');
    expected_id = sprintf('R%06d',r);
    if ~is_char_row(root.root_id) || ~strcmp(root.root_id,expected_id) || ...
            ~is_uint(root.created_discovery_number) || ...
            ~is_char_row(root.target_key) || isempty(root.target_key) || ...
            ~is_real_double_row(root.coordinate) || isempty(root.coordinate) || ...
            any(~isfinite(root.coordinate)) || ...
            ~isequal(root.comparison_mesh,store.settings.comparison_mesh) || ...
            ~valid_profile(root.representative_profile,store.settings) || ...
            ~is_char_row(root.representative_attempt_id) || ...
            isempty(root.representative_attempt_id) || ...
            ~isstruct(root.discoveries) || ~isrow(root.discoveries) || ...
            isempty(root.discoveries)
        error('SA:ROOT_STORE:E_STORE_ROOT', 'Malformed root record.');
    end
    for q = 1:numel(root.discoveries)
        d = root.discoveries(q);
        require_exact_fields(d, discovery_fields, ...
            'SA:ROOT_STORE:E_STORE_DISCOVERY', 'discovery');
        validate_discovery_shape(d, root, store.settings, comparison_fields);
        all_numbers(end+1) = d.discovery_number; %#ok<AGROW>
    end
    numbers = [root.discoveries.discovery_number];
    [first_number, first_index] = min(numbers);
    first = root.discoveries(first_index);
    if root.created_discovery_number ~= first_number || ...
            ~isequal(root.coordinate,first.coordinate) || ...
            ~isequal(root.representative_profile,first.comparison_profile) || ...
            ~strcmp(char(root.representative_attempt_id),char(first.attempt_id))
        error('SA:ROOT_STORE:E_STORE_ROOT', 'Representative does not equal first discovery.');
    end
end
if store.next_root_number ~= numel(store.roots)+1
    error(counter_id, 'next_root_number is inconsistent.');
end
n_discoveries = numel(all_numbers);
if store.next_discovery_number ~= n_discoveries+1 || ...
        ~isequal(sort(all_numbers),1:n_discoveries)
    error(counter_id, 'Discovery numbering is inconsistent.');
end
for r = 1:numel(store.roots)
    if store.roots(r).created_discovery_number ~= ...
            min([store.roots(r).discoveries.discovery_number])
        error(counter_id, 'Root creation number is inconsistent.');
    end
    if r > 1 && store.roots(r).created_discovery_number <= ...
            store.roots(r-1).created_discovery_number
        error(counter_id, 'Root creation order is inconsistent.');
    end
end
replay_history(store);
end

function validate_discovery_shape(d, root, settings, comparison_fields)
id = 'SA:ROOT_STORE:E_STORE_DISCOVERY';
if ~is_uint(d.discovery_number) || ~is_char_row(d.attempt_id) || ...
        isempty(d.attempt_id) || ~isequal(d.coordinate,root.coordinate) || ...
        ~is_char_row(d.seed_family) || ~is_char_row(d.observed_family) || ...
        ~is_char_row(d.route) || ~isstruct(d.validity) || ~isscalar(d.validity) || ...
        ~is_real_double_row(d.returned_mesh) || numel(d.returned_mesh)<2 || ...
        any(~isfinite(d.returned_mesh)) || any(diff(d.returned_mesh)<=0) || ...
        ~valid_returned_profile(d.returned_profile,d.returned_mesh,settings) || ...
        ~valid_profile(d.comparison_profile,settings) || ...
        ~is_char_row(d.member_root_id) || ...
        ~strcmp(d.member_root_id,root.root_id) || ...
        ~isa(d.representative_distance,'double') || ...
        ~isscalar(d.representative_distance) || ...
        ~isfinite(d.representative_distance) || d.representative_distance<0 || ...
        ~isequal(d.tolerance,settings.similarity_tolerance) || ...
        ~isstruct(d.comparisons) || ~isrow(d.comparisons)
    error(id, 'Malformed discovery record.');
end
if settings.comparison_mesh(1)<d.returned_mesh(1) || ...
        settings.comparison_mesh(end)>d.returned_mesh(end)
    error('SA:ROOT_STORE:E_INTERPOLATION_COVERAGE', ...
        'Persisted returned mesh does not cover the comparison mesh.');
end
recomputed = interpolate_profile(d.returned_mesh,d.returned_profile, ...
    settings.comparison_mesh,'SA:ROOT_STORE:E_INTERPOLATION_COVERAGE');
if ~isequal(d.comparison_profile,recomputed)
    error(id, 'Stored comparison profile disagrees with interpolation.');
end
for k = 1:numel(d.comparisons)
    row = d.comparisons(k);
    require_exact_fields(row,comparison_fields,id,'comparison');
    if ~is_char_row(row.root_id) || ...
            ~isa(row.representative_distance,'double') || ...
            ~isscalar(row.representative_distance) || ...
            ~isfinite(row.representative_distance) || row.representative_distance<0 || ...
            ~isa(row.complete_link_max_distance,'double') || ...
            ~isscalar(row.complete_link_max_distance) || ...
            ~isfinite(row.complete_link_max_distance) || ...
            row.complete_link_max_distance<0 || ...
            ~islogical(row.eligible) || ~isscalar(row.eligible)
        error(id, 'Malformed comparison record.');
    end
end
end

function replay_history(store)
id = 'SA:ROOT_STORE:E_STORE_DISCOVERY';
n = store.next_discovery_number-1;
for number = 1:n
    [root_index, discovery] = discovery_by_number(store.roots, number);
    prior = find([store.roots.created_discovery_number] < number);
    same = prior(arrayfun(@(k) same_target(store.roots(k), ...
        store.roots(root_index).target_key,store.roots(root_index).coordinate),prior));
    expected = historical_comparisons(discovery.comparison_profile, ...
        store.roots,same,number,store.settings);
    if ~comparisons_equal(discovery.comparisons,expected)
        error(id, 'Stored historical comparison set disagrees with replay.');
    end
    eligible = find([expected.eligible]);
    created_now = store.roots(root_index).created_discovery_number == number;
    if isempty(eligible)
        if ~created_now || discovery.representative_distance ~= 0
            error(id, 'Discovery should create a new root.');
        end
    else
        rows = expected(eligible);
        key = [[rows.representative_distance].', root_numbers({rows.root_id}).'];
        [~,order] = sortrows(key,[1 2]);
        chosen = rows(order(1));
        if created_now || ~strcmp(char(discovery.member_root_id),chosen.root_id) || ...
                discovery.representative_distance ~= chosen.representative_distance
            error(id, 'Discovery member choice disagrees with replay.');
        end
    end
end
end

function rows = historical_comparisons(profile, roots, indices, number, settings)
rows = empty_comparisons();
for k = 1:numel(indices)
    root = roots(indices(k));
    prior_discoveries = root.discoveries([root.discoveries.discovery_number] < number);
    values = zeros(1,numel(prior_discoveries));
    for q = 1:numel(prior_discoveries)
        values(q) = scaled_rms(profile,prior_discoveries(q).comparison_profile, ...
            settings.field_scale);
    end
    representative_distance = scaled_rms(profile,root.representative_profile, ...
        settings.field_scale);
    maximum = max(values);
    rows(end+1) = struct('root_id',char(root.root_id), ...
        'representative_distance',representative_distance, ...
        'complete_link_max_distance',maximum, ...
        'eligible',maximum<=settings.similarity_tolerance); %#ok<AGROW>
end
end

function rows = comparisons_for(profile, roots, indices, settings)
rows = empty_comparisons();
for k = 1:numel(indices)
    root = roots(indices(k));
    values = zeros(1,numel(root.discoveries));
    for q = 1:numel(root.discoveries)
        values(q) = scaled_rms(profile,root.discoveries(q).comparison_profile, ...
            settings.field_scale);
    end
    representative_distance = scaled_rms(profile,root.representative_profile, ...
        settings.field_scale);
    maximum = max(values);
    rows(end+1) = struct('root_id',char(root.root_id), ...
        'representative_distance',representative_distance, ...
        'complete_link_max_distance',maximum, ...
        'eligible',maximum<=settings.similarity_tolerance); %#ok<AGROW>
end
end

function discovery = discovery_record(candidate, profile, number, root_id, distance, tolerance, comparisons)
discovery = struct('discovery_number',number,'attempt_id',candidate.attempt_id, ...
    'coordinate',candidate.coordinate,'seed_family',candidate.seed_family, ...
    'observed_family',candidate.observed_family,'route',candidate.route, ...
    'validity',candidate.validity,'returned_mesh',candidate.mesh, ...
    'returned_profile',candidate.profile,'comparison_profile',profile, ...
    'member_root_id',root_id,'representative_distance',distance, ...
    'tolerance',tolerance,'comparisons',comparisons);
end

function indices = same_target_indices(roots, key, coordinate)
indices = find(arrayfun(@(r) same_target(r,key,coordinate),roots));
end

function tf = same_target(root, key, coordinate)
tf = strcmp(char(root.target_key),key) && isequal(root.coordinate,coordinate);
end

function profile = interpolate_profile(mesh, values, comparison_mesh, errid)
if comparison_mesh(1)<mesh(1) || comparison_mesh(end)>mesh(end)
    error(errid, 'Returned mesh does not cover the comparison mesh.');
end
profile = zeros(size(values,1),numel(comparison_mesh));
for field = 1:size(values,1)
    profile(field,:) = interp1(mesh,values(field,:),comparison_mesh,'linear');
end
if any(~isfinite(profile),'all')
    error(errid, 'Interpolation returned nonfinite values.');
end
end

function value = scaled_rms(left, right, scale)
delta = (left-right)./scale;
value = sqrt(sum(delta(:).^2)/numel(delta));
end

function [root_index, discovery] = discovery_by_number(roots, number)
root_index = 0; discovery = struct();
for r = 1:numel(roots)
    q = find([roots(r).discoveries.discovery_number]==number,1);
    if ~isempty(q), root_index=r; discovery=roots(r).discoveries(q); return; end
end
error('SA:ROOT_STORE:E_STORE_COUNTER','Missing discovery number.');
end

function tf = comparisons_equal(left, right)
tf = isequal(left,right);
end

function numbers = root_numbers(ids)
numbers = zeros(1,numel(ids));
for k = 1:numel(ids), numbers(k)=sscanf(ids{k},'R%d'); end
end

function tf = valid_profile(value, settings)
tf = isa(value,'double') && isreal(value) && ismatrix(value) && ...
    size(value,1)==numel(settings.field_scale) && ...
    size(value,2)==numel(settings.comparison_mesh) && ...
    all(isfinite(value),'all');
end

function tf = valid_returned_profile(value, mesh, settings)
tf = isa(value,'double') && isreal(value) && ismatrix(value) && ...
    size(value,1)==numel(settings.field_scale) && ...
    size(value,2)==numel(mesh) && all(isfinite(value),'all');
end

function tf = is_real_double_row(value)
tf = isa(value,'double') && isreal(value) && isrow(value);
end

function tf = is_uint(value)
tf = isa(value,'double') && isreal(value) && isscalar(value) && ...
    isfinite(value) && value>=1 && value==fix(value);
end

function tf = is_scalar_text(value)
tf = (ischar(value)&&isrow(value)) || (isstring(value)&&isscalar(value));
end

function tf = is_char_row(value)
tf = ischar(value) && isrow(value);
end

function require_scalar_struct(value, id, label)
if ~isstruct(value) || ~isscalar(value), error(id,'%s must be a scalar struct.',label); end
end

function require_exact_fields(value, expected, id, label)
if ~isequal(sort(fieldnames(value)),sort(expected(:)))
    error(id,'%s has an unexpected field set.',label);
end
end

function rows = empty_roots()
rows = repmat(struct('root_id','','created_discovery_number',NaN, ...
    'target_key','','coordinate',[],'comparison_mesh',[], ...
    'representative_profile',[],'representative_attempt_id','', ...
    'discoveries',empty_discoveries()),1,0);
end

function rows = empty_discoveries()
rows = repmat(struct('discovery_number',NaN,'attempt_id','','coordinate',[], ...
    'seed_family','','observed_family','','route','','validity',struct(), ...
    'returned_mesh',[],'returned_profile',[],'comparison_profile',[], ...
    'member_root_id','','representative_distance',NaN,'tolerance',NaN, ...
    'comparisons',empty_comparisons()),1,0);
end

function rows = empty_comparisons()
rows = repmat(struct('root_id','','representative_distance',NaN, ...
    'complete_link_max_distance',NaN,'eligible',false),1,0);
end
