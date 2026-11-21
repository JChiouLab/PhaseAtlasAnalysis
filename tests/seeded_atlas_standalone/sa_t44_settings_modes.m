function rows = sa_t44_settings_modes(E)
%SA_T44_SETTINGS_MODES  Layered
%   settings package -> model -> mode -> user with a source per leaf;
%   model/mode layers disjoint (both real modules, every mode; a planted
%   overlap is refused); the three new contradictions; unknown leaf and mode;
%   mode change on resume refused at 0 calls; template / from_run round
%   trip; settings_user.mat and settings_effective.txt in the run folder.
%   The stub fixture sets the rescue ladder explicitly (rescue.ring true is
%   sourced 'user'); a mode-sourced leaf is checked on budget.max_wall_s
%   instead.
T = 'T44 settings and modes';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
d = sa_m3_testdir(E, 't44');
% ---- layers and sources (real modules: resolve only, no library) -----------
Mf = sa_model('original_polarity'); Mn = sa_model('guzzo_noxa');
[Sf, Pf] = sa_settings('resolve', struct('model', 'original_polarity', 'run_dir', 'x'), Mf);
[Sn, Pn] = sa_settings('resolve', struct('model', 'guzzo_noxa', 'run_dir', 'x'), Mn);
ok = strcmp(Sf.mode, 'smoke') && Sf.lattice.step == 0.02 && strcmp(Pf.lattice.step, 'model') && ...
    strcmp(Sf.continuation.guess_mesh, 'subsample') && Sf.continuation.guess_nodes == 401 && ...
    strcmp(Pf.continuation.guess_nodes, 'model') && strcmp(Sf.identity.candidate_mesh, 'dense') && ...
    Sf.routes.continuation && strcmp(Pf.routes.continuation, 'mode') && Sf.budget.max_calls == 100 && ...
    strcmp(Pf.budget.max_calls, 'mode') && strcmp(Pf.run_dir, 'user') && ...
    strcmp(Pf.identity.similarity_tolerance, 'package') && strcmp(Pf.software_version, 'derived') && ...
    isequal(Sf.continuation.box, [1 1]) && strcmp(Sf.targets.policy, 'central_nine');
rows = sa_check(rows, T, 'fission yeast smoke: model layer (step 0.02, subsample 401, dense), mode layer (routes on, box [1 1], 100 calls), sources per leaf', ...
    ok, sprintf('step %g (%s), N %d (%s), routes.continuation %d (%s)', Sf.lattice.step, Pf.lattice.step, ...
    Sf.continuation.guess_nodes, Pf.continuation.guess_nodes, Sf.routes.continuation, Pf.routes.continuation));
ok = Sn.lattice.step == 0.15 && Sn.solver.reltol == 1e-8 && Sn.solver.abstol == 1e-10 && ...
    Sn.validity.refine_reltol == 1e-10 && Sn.validity.refine_abstol == 1e-12 && ...
    strcmp(Sn.continuation.guess_mesh, 'subsample') && strcmp(Sn.identity.candidate_mesh, 'comparison') && ...
    strcmp(Pn.solver.reltol, 'model');
rows = sa_check(rows, T, 'Guzzo_noXA: model layer = the tested values (the fixture leaves leave the test fixture)', ...
    ok, sprintf('N %d', Sn.continuation.guess_nodes));
U = struct('model', 'guzzo_noxa', 'run_dir', 'x', 'lattice', struct('step', 0.1), 'mode', 'pilot');
[S2, P2] = sa_settings('resolve', U, Mn);
rows = sa_check(rows, T, 'user overrides the model layer; pilot mode layer applied', S2.lattice.step == 0.1 && ...
    strcmp(P2.lattice.step, 'user') && strcmp(S2.continuation.region, 'cross') && ...
    strcmp(S2.targets.policy, 'pilot_sample') && S2.budget.max_calls == 500 && strcmp(P2.mode, 'user'));
% ---- disjointness ----------------------------------------------------------
dis = true; det = {};
for M = {Mf, Mn}
    for m = sa_modes('modes')
        x = intersect(sa_modes('layer_leaves', M{1}.setting_defaults), sa_modes('layer_leaves', sa_modes('layer', m{1})));
        dis = dis && isempty(x);
        det{end+1} = sprintf('%s/%s %d', M{1}.id, m{1}, numel(x)); %#ok<AGROW>
    end
end
Mbad = Mn; Mbad.setting_defaults.continuation.box = [2 2];
got = err(@() sa_settings('resolve', struct('model', 'guzzo_noxa', 'run_dir', 'x'), Mbad));
Mbad2 = Mn; Mbad2.setting_defaults.budget = struct('max_calls', 5);
got2 = err(@() sa_settings('resolve', struct('model', 'guzzo_noxa', 'run_dir', 'x'), Mbad2));
rows = sa_check(rows, T, 'model and mode layers set disjoint leaves (2 models x 3 modes); a planted overlap / budget in the model layer is refused', ...
    dis && strcmp(got, 'SA:SETTINGS:E_LAYER') && strcmp(got2, 'SA:SETTINGS:E_LAYER'), [strjoin(det, ', ') ' | ' got ' ' got2]);
% ---- contradictions ----------------------------------------------------------
b = struct('model', 'guzzo_noxa', 'run_dir', 'x');
c1 = err(@() sa_settings('resolve', setf(b, 'mode', 'production'), Mn));
c1b = err(@() sa_settings('resolve', setf(setf(b, 'mode', 'production'), 'budget', struct('max_calls', 10)), Mn));
okp = err(@() sa_settings('resolve', setf(setf(b, 'mode', 'production'), 'budget', ...
    struct('max_calls', 10, 'max_wall_s', 60)), Mn));
c2 = err(@() sa_settings('resolve', setf(b, 'pilot', struct('sample_per_axis', 5)), Mn));
c2ok = err(@() sa_settings('resolve', setf(setf(b, 'mode', 'pilot'), 'pilot', struct('sample_per_axis', 5)), Mn));
c3 = err(@() sa_settings('resolve', setf(setf(b, 'targets', struct('policy', 'all_seedable')), ...
    'continuation', struct('box', [])), Mn));
c3p = err(@() sa_settings('resolve', setf(setf(setf(b, 'mode', 'pilot'), 'targets', struct('policy', 'all_seedable')), ...
    'continuation', struct('region', 'box')), Mn));
c3ok = err(@() sa_settings('resolve', setf(b, 'targets', struct('policy', 'all_seedable')), Mn));   % smoke box [1 1]
cid = 'SA:SETTINGS:E_CONTRADICTION';
rows = sa_check(rows, T, ['contradictions refused: production without user budgets (none / one), pilot.* outside pilot, ' ...
    'smoke or pilot with all_seedable over the whole lattice; the allowed neighbours resolve'], ...
    strcmp(c1, cid) && strcmp(c1b, cid) && isempty(okp) && strcmp(c2, cid) && isempty(c2ok) && ...
    strcmp(c3, cid) && strcmp(c3p, cid) && isempty(c3ok), strjoin({c1, c1b, okp, c2, c2ok, c3, c3p, c3ok}, ' / '));
u1 = err(@() sa_settings('resolve', setf(b, 'colour', 1), Mn));
u2 = err(@() sa_settings('resolve', setf(b, 'continuation', struct('regoin', 'box')), Mn));
u3 = err(@() sa_settings('resolve', setf(b, 'mode', 'fast'), Mn));
u4 = err(@() sa_settings('resolve', setf(b, 'continuation', struct('region', 'cross', 'box', [1 1])), Mn));
rows = sa_check(rows, T, 'unknown leaf (top and nested) and unknown mode refused; region cross with a box refused', ...
    strcmp(u1, 'SA:SETTINGS:E_UNKNOWN') && strcmp(u2, 'SA:SETTINGS:E_UNKNOWN') && ...
    strcmp(u3, 'SA:SETTINGS:E_MODE') && strcmp(u4, cid), strjoin({u1, u2, u3, u4}, ' / '));
% ---- run folder files, mode change on resume ----------------------------------
addpath(E.stub);
r = sa_m3_run(fullfile(d, 'smoke'), [1 1]);
Z = load(fullfile(sa_layout('state', fullfile(d, 'smoke')), 'settings_user.mat'));
s0 = sa_m3_stub('settings', fullfile(d, 'smoke'), [1 1]);
txt = fileread(fullfile(d, 'smoke', 'settings_effective.txt'));
rows = sa_check(rows, T, 'settings_user.mat is the struct exactly as passed; settings_effective.txt lists every leaf with its source', ...
    isequaln(Z.settings, s0) && contains(txt, 'continuation.box') && ~isempty(regexp(txt, ...
    'budget\.max_calls\s+5000\s+user', 'once')) && ~isempty(regexp(txt, 'rescue\.ring\s+true\s+user', 'once')) ...
    && ~isempty(regexp(txt, 'budget\.max_wall_s\s+3600\s+mode', 'once')) ...
    && ~isempty(regexp(txt, 'mode\s+''smoke''\s+package', 'once')), sprintf('calls %d', r.calls));
n0 = numel(dir(fullfile(sa_layout('state', fullfile(d, 'smoke')), 'attempts', 'A*.mat')));
g = err(@() sa_m3_run(fullfile(d, 'smoke'), [1 1], 'settings', struct('resume', true, 'mode', 'pilot', 'continuation', struct('region', 'box'))));
n1 = numel(dir(fullfile(sa_layout('state', fullfile(d, 'smoke')), 'attempts', 'A*.mat')));
rows = sa_check(rows, T, 'resume with another mode refused (SA:RESUME:E_SETTINGS), 0 calls', ...
    strcmp(g, 'SA:RESUME:E_SETTINGS') && n0 == n1, g);
% ---- template and from_run ------------------------------------------------------
t1 = sa_settings('template', 'guzzo_noxa', 'production');
t1.run_dir = 'x';
e1 = err(@() sa_settings('resolve', t1, Mn));
t1.budget = struct('max_calls', 50, 'max_wall_s', 600);
e1b = err(@() sa_settings('resolve', t1, Mn));
t2 = sa_settings('template', 'original_polarity', 'pilot'); t2.run_dir = 'x';
e2 = err(@() sa_settings('resolve', t2, Mf));
rows = sa_check(rows, T, 'template: production refused until the user sets budgets; pilot resolves as is', ...
    strcmp(e1, cid) && isempty(e1b) && isempty(e2), strjoin({e1, e1b, e2}, ' / '));
pd = fullfile(d, 'pilotrun'); mkdir(pd);
Upilot = struct('model', 'guzzo_noxa', 'mode', 'pilot', 'run_dir', pd, 'lattice', struct('step', 0.1), ...
    'pilot', struct('sample_per_axis', 5), 'routes', struct('rescue', false), ...
    'budget', struct('max_calls', 400), 'identity', struct('candidate_mesh', 'dense'));
sa_write('mat', fullfile(pd, 'settings_user.mat'), struct('settings', Upilot));
[Up, dropped] = sa_settings('from_run', pd, 'production');
Up.run_dir = 'y'; Up.budget = struct('max_calls', 900, 'max_wall_s', 3600);
[Sp, Pp] = sa_settings('resolve', Up, Mn);
rows = sa_check(rows, T, ['from_run(pilot, production): user leaves kept (step, candidate_mesh), mode-layer leaves, ' ...
    'pilot.*, run_dir and mode dropped; resolves as production'], strcmp(Sp.mode, 'production') && ...
    Sp.lattice.step == 0.1 && strcmp(Pp.lattice.step, 'user') && strcmp(Sp.identity.candidate_mesh, 'dense') && ...
    ~isfield(Up, 'pilot') && ~isfield(Up, 'routes') && strcmp(Sp.targets.policy, 'per_seed') && ...   % production default is per_seed
    isempty(Sp.continuation.box) && all(ismember({'budget.max_calls', 'pilot.sample_per_axis', ...
    'routes.rescue', 'mode', 'run_dir'}, dropped)), strjoin(dropped, ','));
rmpath(E.stub);
end

function id = err(f)
id = '';
try
    f();
catch ME
    id = ME.identifier;
end
end

function s = setf(s, f, v)
s.(f) = v;
end
