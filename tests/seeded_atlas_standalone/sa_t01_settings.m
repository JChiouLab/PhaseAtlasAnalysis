function rows = sa_t01_settings(E)
%SA_T01_SETTINGS  Settings: defaults resolve; unknown/contradictory input rejected.
T = 'T01 settings';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
M = sa_model('original_polarity');
base = struct('model', 'original_polarity', 'run_dir', fullfile(E.scratch, 'unused'));
S = sa_settings('resolve', base, M);
rows = sa_check(rows, T, 'defaults resolve; mesh_reject = nmax = 20000; anchor [0 0]', ...
    S.solver.mesh_reject == 20000 && S.solver.nmax == 20000 && isequal(S.lattice.anchor, [0 0]) ...
    && S.identity.similarity_tolerance == 1e-3 && S.identity.n_comparison_nodes == 2001);
rows = sa_check(rows, T, 'resolved settings carry no validity.fail_closed field', ...
    ~isfield(S.validity, 'fail_closed'));

bad = {
  'unknown top-level field',        setf(base, 'colour', 1),                         'SA:SETTINGS:E_UNKNOWN'
  'unknown nested field',           setf(base, 'solver', struct('foo', 1)),          'SA:SETTINGS:E_UNKNOWN'
  'unknown model option',           setf(base, 'model_options', struct('bar', 1)),   'SA:SETTINGS:E_UNKNOWN'
  'derived field supplied',         setf(base, 'software_version', 'x'),             'SA:SETTINGS:E_DERIVED'
  'unknown mode (M4: smoke, pilot, production)', setf(base, 'mode', 'fast'),     'SA:SETTINGS:E_MODE'
  'model id mismatch',              setf(base, 'model', 'other_model'),              'SA:SETTINGS:E_MODEL'
  'missing run_dir',                rmfield(base, 'run_dir'),                        'SA:SETTINGS:E_FIELD'
  'abstol >= reltol',               setf(base, 'solver', struct('abstol', 1e-5)),    'SA:SETTINGS:E_CONTRADICTION'
  'mesh_reject > nmax',             setf(base, 'solver', struct('mesh_reject', 30000)), 'SA:SETTINGS:E_CONTRADICTION'
  'explicit policy, empty list',    setf(base, 'targets', struct('policy', 'explicit')), 'SA:SETTINGS:E_CONTRADICTION'
  'list with central_nine',         setf(base, 'targets', struct('list', [0 0])),    'SA:SETTINGS:E_CONTRADICTION'
  'tolerance outside band',         setf(base, 'identity', struct('similarity_tolerance', 0.5)), 'SA:SETTINGS:E_CONTRADICTION'
  'refine reltol not tighter',      setf(base, 'validity', struct('refine_reltol', 1e-6, 'refine_abstol', 1e-9)), 'SA:SETTINGS:E_CONTRADICTION'
  'anchor of wrong length',         setf(base, 'lattice', struct('anchor', [0 0 0])), 'SA:SETTINGS:E_FIELD'
  'non-integer target index',       setf(base, 'targets', struct('policy', 'explicit', 'list', [0.5 0])), 'SA:SETTINGS:E_FIELD'
  'hss_branch other rule',          setf(base, 'model_options', struct('hss_branch', 'lowest_u')), 'SA:SETTINGS:E_FIELD'
  'duplicated family',              setf(base, 'families', {'HSS_HSS','HSS_HSS'}), 'SA:SETTINGS:E_CONTRADICTION'
  'negative step',                  setf(base, 'lattice', struct('step', -0.02)),    'SA:SETTINGS:E_FIELD'
  'validity.fail_closed (removed; fail-closed is not optional)', setf(base, 'validity', struct('fail_closed', false)), 'SA:SETTINGS:E_UNKNOWN'
  };
for k = 1:size(bad, 1)
    got = '';
    try
        sa_settings('resolve', bad{k, 2}, M);
    catch ME
        got = ME.identifier;
    end
    rows = sa_check(rows, T, ['rejects: ' bad{k, 1}], strcmp(got, bad{k, 3}), ...
        sprintf('expected %s, got %s', bad{k, 3}, got));
end

% parameter normalization rejects missing / unknown / non-positive fields
p = E.p14;
cases = {rmfield(p, 'M_P'), setf(p, 'L_half', 7), setf(p, 'Dm', -1), setf(p, 'M_P', 30)};
labels = {'missing parameter', 'unknown parameter', 'non-positive parameter', 'M_P >= M_BG'};
for k = 1:numel(cases)
    got = '';
    try, M.normalize(cases{k}); catch ME, got = ME.identifier; end
    rows = sa_check(rows, T, ['normalize rejects: ' labels{k}], strcmp(got, 'SAP:NORMALIZE'), got);
end
end

function s = setf(s, f, v)
s.(f) = v;
end
