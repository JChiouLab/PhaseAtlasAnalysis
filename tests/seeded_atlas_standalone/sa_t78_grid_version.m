function rows = sa_t78_grid_version(E)
%SA_T78_GRID_VERSION  Milestone 7c items 4 and 5: the version label 0.8; the
%   wrapper's points_per_axis (header grid line, refusal of a second place,
%   a folder on another grid refused E_SETTINGS, never mixed); an older
%   version's folder still loads, plots and traces, is refused on resume, and
%   its trace is recomputed once. Zero solves (stub model, stub corrector,
%   stub integrator; the shipped script is read as text).
T = 'T78 grid and version';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d0 = sa_m3_testdir(E, 't78');
v = sa_settings('version');
% P3 amended: 0.10 (model_options.main_variable and the seam rule of the family label); was 0.9
% P4 amended: 0.11 (mode names as the family label, uniform-state flatness); was 0.10
% M8 amended: 0.12 (seed_modes, mode-n libraries, length scaling); was 0.11
% M9 amended: 1.0.0-rc1 (the release self-check of the traced bvp4c, provenance.arch); was 0.12
% M10 amended: 1.0.0-rc2 (import, uncapped defaults); was 1.0.0-rc1
rows = sa_check(rows, T, 'software version is seeded_atlas_standalone/1.0.0-rc2 (no milestone text)', ...
    strcmp(v, 'seeded_atlas_standalone/1.0.0-rc2'), v);

% ---- the shipped script --------------------------------------------------------------------------------------
w = fileread(fullfile(E.root, 'run_seeded_atlas.m'));
live = regexp(w, '(?m)^points_per_axis\s*=\s*51;', 'once');
passed = contains(w, '''points_per_axis'', points_per_axis');
old = regexp(w, '(?m)^%\s*settings_overrides\.lattice\.points_per_axis', 'once');
rows = sa_check(rows, T, ['run_seeded_atlas.m: a live points_per_axis = 51 line, passed to the wrapper; no commented ' ...
    'settings_overrides.lattice.points_per_axis line left; the way back to a fixed step is stated'], ...
    ~isempty(live) && passed && isempty(old) && contains(w, '[] = a fixed step'));
% M9 amended (E8): two added lines put the script's own folder on the path and call setup_atlas('quiet')
% after the package is added, before the wrapper
i_pkg = strfind(w, 'addpath(fullfile(here, ''src'', ''pipelines'', ''seeded_atlas_standalone''));');
i_e8 = regexp(w, '(?m)^addpath\(here\);\r?\nsetup_atlas\(''quiet''\);', 'once');
i_w = strfind(w, 'W = seeded_atlas_wrapper(cfg);');
rows = sa_check(rows, T, ['run_seeded_atlas.m (E8): "addpath(here);" and "setup_atlas(''quiet'');" on two consecutive ' ...
    'lines after the package addpath and before the wrapper call'], isscalar(i_pkg) && ~isempty(i_e8) && ...
    isscalar(i_w) && i_pkg < i_e8 && i_e8 < i_w && isfile(fullfile(E.root, 'setup_atlas.m')), ...
    sprintf('package addpath at %s, E8 at %s, wrapper at %s', mat2str(i_pkg), mat2str(i_e8), mat2str(i_w)));

% ---- wrapper grid line and precedence ------------------------------------------------------------------------
s0 = sa_m3_stub('settings', 'x', [1 1]);
ov = rmfield(s0, {'model', 'run_dir', 'budget'});
ovp = ov; ovp.lattice = rmfield(ovp.lattice, 'step');
base = struct('model', 'stub_two_pool', 'par', sa_m3_stub('params'), 'mode', 'smoke', ...
    'library_manifest', struct(), 'settings_overrides', ovp, 'budget', struct('max_calls', 5000), 'status_only', true);
c = base; c.run_dir = fullfile(d0, 'g5'); c.points_per_axis = 5;
out = evalc('W5 = seeded_atlas_wrapper(c);');
gl = regexp(out, ' grid [^\n]*', 'match', 'once');
c = base; c.run_dir = fullfile(d0, 'g0'); c.settings_overrides = ov;
out0 = evalc('W0 = seeded_atlas_wrapper(c);');
g0 = regexp(out0, ' grid [^\n]*', 'match', 'once');
rows = sa_check(rows, T, ['the wrapper''s header prints the grid: points_per_axis 5 with the derived per-axis steps and ' ...
    '5 x 5 = 25 nodes; with [] the fixed step'], contains(gl, 'points_per_axis 5') && contains(gl, '5 x 5 = 25 nodes') && ...
    isequal(W5.status.settings.lattice.points_per_axis, 5) && W5.lattice.n_nodes == 25 && ...
    contains(g0, 'fixed step 0.65') && isempty(W0.status.settings.lattice.points_per_axis), sprintf('%s || %s', gl, g0));
c = base; c.run_dir = fullfile(d0, 'g2'); c.points_per_axis = 5;
c.settings_overrides.lattice.points_per_axis = 7;
id = ''; try, seeded_atlas_wrapper(c); catch ME, id = ME.identifier; end
rows = sa_check(rows, T, 'the grid set both in the edit block and in settings_overrides is refused (E_FIELD)', ...
    strcmp(id, 'SA:WRAPPER:E_FIELD'), id);

% ---- a folder on the fixed step, then the count: refused, never mixed ------------------------------------------
st = sa_m3_stub('state');
dS = fullfile(d0, 'step');
c = base; c.status_only = false; c.run_dir = dS; c.settings_overrides = ov; c.make_plot = false;
c.corrector = sa_m3_stub('corrector', st);
evalc('seeded_atlas_wrapper(c);');
L0 = listing(dS);
n0 = sa_m3_stub('calls', st);
c2 = c; c2.settings_overrides = ovp; c2.points_per_axis = 5;
c2.status_only = true;
evalc('Wq = seeded_atlas_wrapper(c2);');
c2.status_only = false; id2 = '';
try, evalc('seeded_atlas_wrapper(c2);'); catch ME, id2 = ME.identifier; end
rows = sa_check(rows, T, ['a folder run on the fixed step, opened with points_per_axis 5: status refused ' ...
    'SA:RESUME:E_SETTINGS naming lattice.points_per_axis; the wrapper stops (E_REFUSED), zero calls, the folder ' ...
    'byte-unchanged'], strcmp(Wq.status.state, 'refused') && strcmp(Wq.status.id, 'SA:RESUME:E_SETTINGS') && ...
    contains(Wq.status.message, 'lattice.points_per_axis') && strcmp(id2, 'SA:WRAPPER:E_REFUSED') && ...
    sa_m3_stub('calls', st) == n0 && isequal(L0, listing(dS)), Wq.status.message);

% ---- an at entry through the wrapper (the edit block's traces) ---------------------------------------------------
Qs = seeded_atlas_plot(dS, struct('coordinates_only', true));
cw = containers.Map({'n'}, {0});
c3 = c; c3.corrector = @(varargin) error('T78:CALLED', 'corrector called on a completed run');
c3.integrator = @(theta, Y0, entry) stub_int(theta, Y0, entry, cw);
c3.traces = struct('name', 'P', 'root', [], 'at', Qs.coord(end, :));
evalc('Wt = seeded_atlas_wrapper(c3);');
rows = sa_check(rows, T, 'an at entry in the wrapper''s traces (root left empty) runs and resolves to the root at that position', ...
    Wt.traces(1).ran && strcmp(Wt.traces(1).at_resolution.root_id, Qs.root_id{end}) && cw('n') == 1, ...
    Wt.traces(1).reason);

% ---- an older version's folder ------------------------------------------------------------------------------------
dO = fullfile(d0, 'older');
sa_m3_run(dO, [1 1]);
R = seeded_atlas_roots(dO, struct('verify', false));
cnt = containers.Map({'n'}, {0});
o = struct('integrator', @(theta, Y0, entry) stub_int(theta, Y0, entry, cnt));
eA = struct('name', 'A', 'root', R(1).root_id);
seeded_atlas_trace(dO, eA, o);
% make it a folder of the previous version: the stored version text and the trace key's
D = sa_layout('state', dO);
Z = load(fullfile(D, 'settings_effective.mat'));
pv = 'seeded_atlas_standalone/0.7 (milestone 7b)';
Z.settings.software_version = pv; Z.provenance.software_version = pv;
sa_write('mat', fullfile(D, 'settings_effective.mat'), Z);
for f = dir(fullfile(D, 'attempts', '*.mat')).'
    A = load(fullfile(f.folder, f.name)); A.provenance.software_version = pv;
    sa_write('mat', fullfile(f.folder, f.name), A);
end
tf = fullfile(dO, 'traces', 'A', 'trace.mat');
Zt = load(tf); Zt.key.software_version = pv;
Zt.key.settings_sha256 = sa_sha256(fullfile(D, 'settings_effective.mat'));
Zt.key.attempt_sha256 = sa_sha256(fullfile(D, 'attempts', [Zt.key.attempt_id '.mat']));
sa_write('mat', tf, Zt);
s = sa_m3_stub('settings', dO, [1 1]); s.resume = true;
M = sa_model('stub_two_pool');
cur = sa_resume('current', M, sa_m3_stub('params'), struct());
Zs = sa_resume('status', dO, sa_settings('resolve', s, M), cur);
Rl = seeded_atlas_roots(dO, struct('verify', true));
pl = seeded_atlas_plot(dO);
n1 = cnt('n');
t1 = seeded_atlas_trace(dO, eA, o);
t2 = seeded_atlas_trace(dO, eA, o);
rows = sa_check(rows, T, ['an older version''s folder: status refused (software version), roots load verified, the ' ...
    'plot redraws, the trace is recomputed once (reason: software_version) and then reused'], ...
    strcmp(Zs.state, 'refused') && contains(Zs.message, 'software') && numel(Rl) == numel(R) && pl.redrawn && ...
    t1.ran && contains(t1.reason, 'software_version') && ~t2.ran && cnt('n') == n1 + 1, ...
    sprintf('%s | %s', Zs.id, t1.reason));
end

% =========================================================================

function L = listing(d)
f = dir(fullfile(d, '**', '*'));
f = f(~[f.isdir]);
L = sort(arrayfun(@(x) sprintf('%s|%d|%s', fullfile(x.folder, x.name), x.bytes, sa_sha256(fullfile(x.folder, x.name))), ...
    f, 'UniformOutput', false));
end

function fr = stub_int(theta, Y0, entry, cnt)
cnt('n') = cnt('n') + 1;
P = sa_post_model('stub_two_pool');
t = entry.t_out; F = numel(t);
x = linspace(0, 2 * theta.Lh, size(Y0, 2));
O = P.observables(theta, x, Y0);
fr = struct('t', t, 'coord', repmat(O.c, F, 1), 'masses', repmat(O.masses, F, 1), ...
    'umax', repmat(O.umax, F, 1), 'rhs_inf', nan(1, F), 'profile_frames', [1 F], 'stop', 'horizon_reached', ...
    'stop_detail', '', 'wall_s', 0, 'x', x);
fr.profiles = repmat(Y0, 1, 1, 2);
fr.mass_drift_rel = zeros(F, size(O.masses, 2));
fr.Y_final = Y0;
end
