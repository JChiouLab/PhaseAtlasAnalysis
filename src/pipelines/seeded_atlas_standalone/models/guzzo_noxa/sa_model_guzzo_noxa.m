function M = sa_model_guzzo_noxa()
%SA_MODEL_GUZZO_NOXA  Model module of the M. xanthus model (coupled problem).
%
%   M = SA_MODEL_GUZZO_NOXA() returns the module struct that
%   sa_model('guzzo_noxa') validates. It puts the shared model core
%   src/models/guzzo_noxa/ (gnx_model, gnx_reaction; read only) on the path
%   itself, relative to this file. It never puts the uncoupled pipeline's
%   guzzo_noxa folder on the path (the two packages do not depend on each
%   other): the two functions it needs from there are verbatim copies
%   (sag_hermite_min, sag_bvp_convergence).
%
%   Coupled state (sag_coupled24, 24 rows, x in [0, Lh], dimensional um):
%     1:8 left half Rm Rm' Rc Rc' Bm Bm' Bc Bc'; 9:16 right half;
%     17:20 sources s_RL s_BL s_RR s_BR (named sigma in the code); 21:24 integrals of
%     the half masses.
%   Conserved quantities R = Rm + Rc (grid axis 1, dn_R), B = Bm + Bc (axis 2, dn_B).

here = fileparts(mfilename('fullpath'));
src = fileparts(fileparts(fileparts(fileparts(here))));
core = fullfile(src, 'models', 'guzzo_noxa');
if ~contains([pathsep path pathsep], [pathsep core pathsep])
    addpath(core);
end

copies = copy_records(here);
M = struct();
M.id        = 'guzzo_noxa';
M.version   = module_version(here, core);
M.unclipped = true;                      % gnx_reaction and the ode subfunction apply no max(0,.); positivity is checked afterwards
M.n_states  = 24;
M.fields    = named_rows({'Rm_L','Rc_L','Bm_L','Bc_L','Rm_R','Rc_R','Bm_R','Bc_R'}, ...
                         [1 3 5 7 9 11 13 15]);
M.identity_fields = M.fields;
M.pools     = struct('name', {'R', 'B'});          % the two conserved quantities
M.library_types = {'HSS', 'HPSS', 'MPSS'};
M.domain    = @(theta) [0 theta.Lh];
M.normalize = @normalize;
M.pool_targets = @pool_targets;
M.library   = @sag_library;
M.pair_guess = @(lib, theta, tL, cL, tR, cR, mo) sag_library('pair_guess', lib, theta, tL, cL, tR, cR, mo);
M.ode       = @(x, Y, sol, theta) sag_coupled24('ode', x, Y, theta);
M.solve     = @sag_solve;
M.convergence = @sag_bvp_convergence;
M.diagnostics = @sag_checks;
M.shifts    = @shifts;
M.shape_label = @shape_label;
M.shape_label_profile = @shape_label_profile;
M.mode_families = mode_families();
M.seed_shape = @sag_seed_shape;          % seed shapes, type order, mirror rule
% The left/right half swap (halves, sources and integrals in blocks; the layout of
% sag_coupled24) and the permutation it induces on the 8 identity rows (L then R)
M.reflection = struct('state_perm', [9:16 1:8 19 20 17 18 23 24 21 22], 'identity_perm', [5:8 1:4]);
M.scale_counts = @scale_counts;          % point counts that grow with L
M.field_scale = @field_scale;
% the settings leaves field_scale reads besides the parameters (a focus run
% takes the source's field scale only when they equal the source's)
M.field_scale_settings = {'lattice.anchor', 'identity.field_scale_floor', 'model_options.hss_branch'};
M.comparison_mesh = @(theta, S) linspace(0, theta.Lh, S.identity.n_comparison_nodes);
M.required_checks = {'real_valued', 'positivity_exact', 'seam_concentration', 'seam_flux', ...
    'mass_accumulator', 'mass_dense', 'mass_dense_refinement', 'shift_antisymmetry', ...
    'velocity_consistency'};
M.model_options_defaults = option_defaults();
M.validate_model_options = @validate_options;
M.copies    = copies;
M.setting_defaults = setting_defaults();     % model layer of the settings (optional member)
M.version   = [M.version setting_defaults_tag(M.setting_defaults)];
% The 'solver:' component: SHA-256
% prefixes of the traced bvp4c copy, its five helpers, sa_deval, the toolbox
% bvp4c.m / deval.m and the forwarders; verifies the generated copy first
% (SA:SOLVER:E_TRACED / E_TOOLBOX at zero solves).
M.version   = [M.version sa_bvp4c_setup('component', here)];
end

% =========================================================================

function sd = setting_defaults()
% Model layer: the shared settings at the model's default values (step
% 0.15; RelTol 1e-8 / AbsTol 1e-10; refinement 1e-10 / 1e-12; guess mesh
% 'subsample' with N = 482 = 2 x 241; candidate mesh 'comparison').
% Applied after the package defaults and before the run-mode layer.
sd = struct();
sd.lattice = struct('step', 0.15);
sd.solver = struct('reltol', 1e-8, 'abstol', 1e-10);
sd.validity = struct('refine_reltol', 1e-10, 'refine_abstol', 1e-12);
sd.continuation = struct('guess_mesh', 'subsample', 'guess_nodes', 482);
sd.identity = struct('candidate_mesh', 'comparison');
sd.report = struct('replay_s_per_record', 0.05);   % report-only replay rate (seconds per stored record on a resume)
% Cold-only mesh cap NMax = mesh_reject = 10,000 (Inf turns it off) and the
% per-root continuation pass (off for this model; root_pass.enabled = true
% turns it on).
sd.cold = struct('mesh_cap', 10000);
sd.root_pass = struct('enabled', false);
% The reference length of the point-count rule (L = 6, the value of Table 1;
% scale_counts). Enters the model digest.
sd.scaling = struct('reference_length', 6);
end

function [S, SRC] = scale_counts(theta, S, SRC)
% Point counts that grow with L. The length leaf
% is the total length L (theta.L = P.L; Lh = L/2). At L <= reference_length
% (on the stored double) nothing runs and every leaf keeps its source. Above
% it, for every leaf whose source is package or model (user values win):
% node counts N(L) = 1 + ceil((N_c - 1) L / L_ref), mesh limits
% ceil(N_c L / L_ref) (Inf stays Inf), source 'derived'.
L = theta.L;
Lref = S.scaling.reference_length;
if ~(L > Lref), return; end
nodes = {'model_options', 'guess_nodes'; 'continuation', 'guess_nodes'; 'model_options', 'shape_nodes'};
mesh  = {'solver', 'nmax'; 'solver', 'mesh_reject'; 'cold', 'mesh_cap'};
for k = 1:size(nodes, 1)
    [S, SRC] = scale_leaf(S, SRC, nodes(k, :), @(N) 1 + ceil((N - 1) * L / Lref));
end
for k = 1:size(mesh, 1)
    [S, SRC] = scale_leaf(S, SRC, mesh(k, :), @(N) ceil(N * L / Lref));
end
end

function [S, SRC] = scale_leaf(S, SRC, path, rule)
src = SRC.(path{1}).(path{2});
if ~any(strcmp(src, {'package', 'model'})), return; end
N = S.(path{1}).(path{2});
if ~isfinite(N), return; end
S.(path{1}).(path{2}) = rule(N);
SRC.(path{1}).(path{2}) = 'derived';
end

function t = setting_defaults_tag(sd)
% model_version carries a digest of the model layer (the module files
% change with the member). report.* is excluded from the
% digest: a report-only leaf never changes a solve, so recalibrating it
% never makes runs unresumable.
if isfield(sd, 'report'), sd = rmfield(sd, 'report'); end
txt = leaf_text(sd, '');
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(txt), 'int8'));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
t = sprintf(' setting_defaults:%s', h(1:12));
end

function t = leaf_text(s, path)
t = '';
f = fieldnames(s);
for k = 1:numel(f)
    p = [path '.' f{k}];
    v = s.(f{k});
    if isstruct(v), t = [t leaf_text(v, p)]; %#ok<AGROW>
    elseif ischar(v), t = [t p '=' v ';']; %#ok<AGROW>
    else, t = [t p '=' mat2str(v, 17) ';']; %#ok<AGROW>
    end
end
end


% =========================================================================

function v = module_version(here, core)
L = dir(fullfile(here, 'sag_*.m'));
names = sort({L.name});
parts = cell(1, numel(names));
for k = 1:numel(names)
    h = sa_sha256(fullfile(here, names{k}));
    parts{k} = sprintf('%s:%s', names{k}, h(1:12));
end
for f = {'gnx_model.m', 'gnx_reaction.m'}
    h = sa_sha256(fullfile(core, f{1}));
    parts{end+1} = sprintf('%s:%s', f{1}, h(1:12)); %#ok<AGROW>
end
v = sprintf('guzzo_noxa_coupled/1 (%s)', strjoin(parts, ','));
end

function theta = normalize(p)
% gnx_model('parameters') (rejects unknown or missing fields: L, r_speed,
% dba_factor, rho_R, rho_B, res, dn_half_range) plus the half domain Lh.
theta = gnx_model('parameters', p);
theta.Lh = theta.L / 2;
end

function s = named_rows(names, rows)
s = struct('name', names, 'row', num2cell(rows));
end

function t = pool_targets(theta, d)
% d = [dn_R dn_B] (grid axis order); values [n_RL n_RR n_BL n_BR].
tgt = sag_coupled24('targets', theta, d(1), d(2));
t = struct();
t.values = [tgt.nRL tgt.nRR tgt.nBL tgt.nBR];
t.labels = {'n_RL', 'n_RR', 'n_BL', 'n_BR'};
t.solver_target = tgt;
end

function s = shifts(sol)
s = struct('names', {{'sigma_RL', 'sigma_BL', 'sigma_RR', 'sigma_BR'}}, ...
    'values', sol.y(17:20, 1).');
s.max = max(abs(s.values));
end

function [lab, info] = shape_label(sol, theta, mo, fs) %#ok<INUSL>
% The shape label of a solution: the mode name (sa_mode_name) of the main
% variable (model_options.main_variable, default Bm; a sum of fields is summed
% per half, in the given order) in each half, with the flatness reference
% ref = the main variable's value in the uniform (HSS) state at the anchor (fs:
% field_scale's struct, its anchor_hss_values: the library HSS interpolated
% to the anchor). shape_nodes points per half (x = 0 the pole,
% x = Lh the middle). It never affects validity or root identity. It sets the
% attempt's observed_family and, through the per-root pass, that pass's
% unit family string and the child's seed_family, so it is part of the run's
% settings; it never chooses a starting profile, a node, an order or a cap.
% info: sa_mode_name's info, plus .legacy, the older three-state label (HSS /
% HPSS / MPSS per half: flat when the relative contrast (max - min)/max <
% 0.2, a maximum at either end of the half -> HPSS).
xi = linspace(sol.x(1), sol.x(end), mo.shape_nodes);
v = sa_main_variable('parse', mo.main_variable, species());
n = numel(v.index);
Y = deval(sol, xi, [2*v.index - 1, 2*v.index + 7]);   % rows 1 3 5 7, then 9 11 13 15
[lab, info] = label_pair(sa_main_variable('sum', Y, 1:n), sa_main_variable('sum', Y, n + (1:n)), mo, ...
    sa_mode_name('reference', fs, v.index));
end

function [lab, info] = shape_label_profile(xi, P, theta, mo, fs) %#ok<INUSL>
% The same rule on a stored profile instead of a solver solution: P holds
% the identity rows (left Rm Rc Bm Bc, then right) on the mesh xi, e.g. a
% root's representative profile on the comparison mesh. Each field the
% main variable needs is interpolated linearly to shape_label's
% shape_nodes points, then summed in the same order. For post-processing
% only (a plot relabelled with another main variable); a run's labels come
% from shape_label.
xs = linspace(xi(1), xi(end), mo.shape_nodes);
v = sa_main_variable('parse', mo.main_variable, species());
Y = zeros(8, numel(xs));
for r = [v.index, v.index + 4]
    Y(r, :) = interp1(xi, P(r, :), xs, 'linear');
end
[lab, info] = label_pair(sa_main_variable('sum', Y, v.index), sa_main_variable('sum', Y, v.index + 4), mo, ...
    sa_mode_name('reference', fs, v.index));
end

function s = species()
s = {'Rm', 'Rc', 'Bm', 'Bc'};
end

function f = mode_families()
% The library half type that a half's mode name corresponds to (the per-root
% pass's family lookup): flat -> HSS (no-peak); a pole peak of mode 1 (left 1+,
% right 1-) -> HPSS (half-peak); a middle peak 2- -> MPSS (full-peak). Other names have none.
f = struct('left', {{'0', 'HSS'; '1+', 'HPSS'; '2-', 'MPSS'}}, ...
    'right', {{'0', 'HSS'; '1-', 'HPSS'; '2-', 'MPSS'}});
end

function [lab, info] = label_pair(mL, mR, mo, ref)
[lab, info] = sa_mode_name('label', mL, mR, ref, mo);
info.legacy = sprintf('%s_%s', legacy_shape(mL, mo), legacy_shape(mR, mo));
end

function s = legacy_shape(b, mo)
% One half's older label (tests only): HSS when
% the relative contrast (max - min)/max < 0.2, HPSS
% when the maximum lies within shape_edge_frac of either end, otherwise
% MPSS.
mx = max(b);
if ~(mx > 0) || (mx - min(b)) / mx < 0.2, s = 'HSS'; return; end
n = numel(b); e = max(1, round(mo.shape_edge_frac * n));
[~, k] = max(b);
if k <= e || k >= n - e + 1, s = 'HPSS'; else, s = 'MPSS'; end
end

function fs = field_scale(theta, lib, S)
% Per field max(|value| of the library HSS at the anchor cell, floor), one
% value per species, applied to both halves. hss_values: the same HSS
% values, raw; anchor_hss_values: the library HSS interpolated to the
% anchor itself (the uniform-state reference of the shape label, sa_mode_name
% 'reference').
lib = sag_library('profiles', lib, {'HSS'});
c = sag_library('cell', lib, theta, S.lattice.anchor);
if ~c.in_range
    error('SAG:FIELD_SCALE', 'The anchor lies outside the library''s native range.');
end
h = sag_library('half', lib, theta, 'HSS', c.cell, S.model_options);
if ~h.available
    error('SAG:FIELD_SCALE', 'No HSS half at the anchor cell (%d,%d): %s.', ...
        c.cell(1), c.cell(2), h.reason);
end
v4 = max(abs(h.Y4(:, 1)), S.identity.field_scale_floor);
fs = struct();
fs.values = [v4; v4];
fs.species_values = v4.';
fs.species = {'Rm', 'Rc', 'Bm', 'Bc'};
fs.source = sprintf(['library HSS half, hss_branch = %s, native cell (%d,%d) at anchor ' ...
    '[%s], model units, floor %g'], S.model_options.hss_branch, c.cell(1), c.cell(2), ...
    num2str(S.lattice.anchor), S.identity.field_scale_floor);
fs.hss_values = h.Y4(:, 1).';
% the uniform-state reference of the shape label: the library HSS interpolated
% bilinearly to the anchor itself (the native grid has no node there in
% general); the nearest cell's values when a surrounding cell has none
val = @(i, j) hss_at(lib, theta, [i j], S.model_options);
[v, why] = sa_mode_name('bilinear', lib.axis_row, lib.axis_col, S.lattice.anchor, val);
if isempty(v)
    fs.anchor_hss_values = fs.hss_values;
    fs.anchor_source = sprintf('nearest native cell (%d,%d): bilinear interpolation impossible, %s', ...
        c.cell(1), c.cell(2), why);
else
    fs.anchor_hss_values = v.';
    fs.anchor_source = sprintf('library HSS interpolated bilinearly to the anchor [%s]', num2str(S.lattice.anchor));
end
end

function y = hss_at(lib, theta, cell, mo)
h = sag_library('half', lib, theta, 'HSS', cell, mo);
if h.available, y = h.Y4(:, 1); else, y = []; end
end

function d = option_defaults()
% Model-owned settings: default values of the numerical tolerances of the
% validity checks.
d = struct();
d.hss_branch      = 'single_root';   % fail closed unless the HSS cell has one root
d.guess_nodes     = 241;             % mesh points of the initial guess profile
d.seam_tol        = 1e-9;            % seam_concentration, seam_flux
d.acc_tol         = 1e-12;           % mass_accumulator (uncoupled acc_end_rel)
d.mass_dense_tol  = 1e-7;            % mass_dense (uncoupled mass_rel), exact Hermite integral
d.mass_refine_tol = 1e-12;           % mass_dense_refinement (bisected)
d.antisym_abs     = 1e-9;            % shift_antisymmetry gate = max(abs, rel*s_max)
d.antisym_rel     = 1e-6;
d.vel_abs         = 1e-9;            % velocity_consistency gate = max(abs, rel*max|J|)
d.vel_rel         = 1e-6;
d.main_variable   = 'Bm';            % the label's and dn2's variable: a '+'-sum of Rm Rc Bm Bc
d.shape_flat_tol  = 0.2;             % label: a half is 0 (flat) when (max - min) / (uniform value) is below;
                                     %   at 0.2 the depleted halves of the f = 1 runs, whose whole
                                     %   amplitude is a few per cent of the uniform Bm, are flat
d.shape_turn_frac = 0.2;             % label: a turn counts when the profile moves back by this x (max - min)
d.shape_edge_frac = 0.02;            % label: the seam window, as a fraction of a half
d.shape_nodes     = 401;
end

function mo = validate_options(mo)
id = 'SA:SETTINGS:E_FIELD';
v = sa_main_variable('parse', mo.main_variable, species());
mo.main_variable = v.text;               % canonical form (no whitespace)
if ~strcmp(mo.hss_branch, 'single_root')
    error(id, 'model_options.hss_branch must be ''single_root'' (the only rule implemented).');
end
for f = {'guess_nodes', 'shape_nodes'}
    v = mo.(f{1});
    if ~(isa(v, 'double') && isscalar(v) && v >= 3 && v == fix(v))
        error(id, 'model_options.%s must be an integer >= 3.', f{1});
    end
end
for f = {'seam_tol', 'acc_tol', 'mass_dense_tol', 'mass_refine_tol', 'antisym_abs', ...
         'antisym_rel', 'vel_abs', 'vel_rel', 'shape_flat_tol', 'shape_edge_frac'}
    v = mo.(f{1});
    if ~(isa(v, 'double') && isscalar(v) && isfinite(v) && v > 0)
        error(id, 'model_options.%s must be a finite positive scalar.', f{1});
    end
end
v = mo.shape_turn_frac;
if ~(isa(v, 'double') && isscalar(v) && v > 0 && v < 1)
    error(id, 'model_options.shape_turn_frac must be a scalar in (0, 1).');
end
end

function c = copy_records(here)
% Provenance of the verbatim copies (only the function-name line differs);
% sources and their SHA-256 are checked by the portable tests.
names = {'sag_hermite_min.m', 'sag_bvp_convergence.m'};
c = struct('file', {}, 'copy_sha256', {});
for k = 1:numel(names)
    c(end+1) = struct('file', names{k}, ...
        'copy_sha256', sa_sha256(fullfile(here, names{k}))); %#ok<AGROW>
end
end
