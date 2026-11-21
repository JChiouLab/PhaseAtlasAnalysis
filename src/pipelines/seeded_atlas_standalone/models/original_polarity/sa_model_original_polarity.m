function M = sa_model_original_polarity()
%SA_MODEL_ORIGINAL_POLARITY  Model module of the S. pombe model (five fields u q v w p).
%
%   M = SA_MODEL_ORIGINAL_POLARITY() returns the module struct that
%   sa_model('original_polarity') validates. Everything model-specific the
%   shared pipeline needs comes from here; shared code
%   never names a field, conserved quantity, half type or model constant.
%
%   Coupled state (sap_forced28, 28 rows, x in [0, 1/2] in units of L):
%     rows 1:10   left half  [u u' q q' v v' w w' p p']
%     rows 11:20  right half, same order
%     rows 21:24  sources s_GL s_PL s_GR s_PR (constant states)
%     rows 25:28  integrals of the half masses I_GL I_PL I_GR I_PR (accumulators)
%   Conserved quantities: PAK = p + q (grid axis 1), GEF = u + q + v + w (grid axis 2), so
%   the grid index [i j] is (d_P index, d_G index).

here = fileparts(mfilename('fullpath'));
copies = copy_records(here);

M = struct();
M.id        = 'original_polarity';
M.version   = sprintf('original_polarity/1 (corrector copies %s)', ...
    strjoin(cellfun(@(h) h(1:12), {copies.copy_sha256}, 'UniformOutput', false), ','));
M.unclipped = true;                      % reaction has no max(.,0); checked by a source scan test
M.n_states  = 28;
M.fields    = named_rows({'u_L','q_L','v_L','w_L','p_L','u_R','q_R','v_R','w_R','p_R'}, ...
                         [1 3 5 7 9 11 13 15 17 19]);
M.identity_fields = M.fields;            % the 10 concentrations
M.pools     = struct('name', {'P', 'G'});
M.library_types = {'HSS', 'HPSS', 'MPSS'};
M.domain    = @(theta) [0 0.5];
M.normalize = @sap_normalize;
M.pool_targets = @pool_targets;
M.library   = @sap_library;
M.pair_guess = @(lib, theta, tL, cL, tR, cR, mo) sap_library('pair_guess', lib, theta, tL, cL, tR, cR, mo);
M.ode       = @ode_matrix;
M.solve     = @sap_solve;
M.convergence = @sap_bvp_convergence;
M.diagnostics = @sap_checks;
M.shifts    = @shifts;
M.shape_label = @shape_label;
M.shape_label_profile = @shape_label_profile;
M.mode_families = mode_families();
M.seed_shape = @sap_seed_shape;          % seed shapes, type order, mirror rule
M.scale_counts = @scale_counts;          % point counts that grow with L
% The left/right half swap (sap_forced28 'swap', a pure block
% swap; the permutation lives here) and the
% permutation it induces on the 10 identity rows (L then R)
M.reflection = struct('state_perm', [11:20 1:10 23 24 21 22 27 28 25 26], 'identity_perm', [6:10 1:5]);
M.field_scale = @field_scale;
% the settings leaves field_scale reads besides the parameters (a focus run
% takes the source's field scale only when they equal the source's)
M.field_scale_settings = {'lattice.anchor', 'identity.field_scale_floor', 'model_options.hss_branch'};
M.comparison_mesh = @(theta, S) linspace(0, 0.5, S.identity.n_comparison_nodes);
M.required_checks = {'seam_concentration','seam_flux','mass_accumulator', ...
    'mass_dense','mass_dense_refinement','shift_antisymmetry','velocity_consistency'};
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
% 0.02; RelTol 1e-6 / AbsTol 1e-8; refinement 1e-8 / 1e-10; guess mesh
% 'subsample' with N = 401 (802 as an option for harder cases); candidate mesh
% 'dense'). Applied after the package defaults and before the run-mode layer.
sd = struct();
sd.lattice = struct('step', 0.02);
sd.solver = struct('reltol', 1e-6, 'abstol', 1e-8);
sd.validity = struct('refine_reltol', 1e-8, 'refine_abstol', 1e-10);
sd.continuation = struct('guess_mesh', 'subsample', 'guess_nodes', 401);
sd.identity = struct('candidate_mesh', 'dense');
sd.report = struct('replay_s_per_record', 0.35);   % report-only replay rate (seconds per stored record on a resume)
% Cold-only mesh cap NMax = mesh_reject = 10,000 (Inf turns it off) and the
% per-root continuation pass (on for this model).
sd.cold = struct('mesh_cap', 10000);
sd.root_pass = struct('enabled', true);
% The reference length of the point-count rule (the largest
% tested L; scale_counts). Enters the model digest.
sd.scaling = struct('reference_length', 14);
end

function [S, SRC] = scale_counts(theta, S, SRC)
% Point counts that grow with L. The length leaf
% is the total length L (theta.L = par.L). At L <= reference_length (on the
% stored double) nothing runs and every leaf keeps its source. Above it,
% for every leaf whose source is package or model (user values win):
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

function s = named_rows(names, rows)
s = struct('name', names, 'row', num2cell(rows));
end

function t = pool_targets(theta, d)
% d = [d_P d_G] (grid axis order). Half targets of the coupled solve,
% (M_c +/- d_c)/2 with M_G = 1 and M_P = rho.
tgt = sap_forced28('targets', theta, d(1), d(2));
t = struct();
t.values = [tgt.n_PL tgt.n_PR tgt.n_GL tgt.n_GR];
t.labels = {'n_PL','n_PR','n_GL','n_GR'};
t.solver_target = tgt;
end

function F = ode_matrix(x, Y, sol, theta) %#ok<INUSL>
% The right-hand side of the coupled BVP (the residual check needs F), for all
% columns at once. Every entry is computed by the same floating-point
% operations, in the same order, as sap_forced28('ode', x(k), Y(:, k), theta)
% on column k (elementwise only: no reduction, no reordering), so F is
% bit-identical to the column-by-column evaluation (checked by a portable
% test on stored records).
th = theta;
mu = sap_model('mu', th);
F = zeros(size(Y));
sh = [21 22; 23 24];          % [sG sP] per half
for hh = 1:2
    o  = (hh-1)*10;
    sG = Y(sh(hh,1), :); sP = Y(sh(hh,2), :);
    u = Y(o+1, :); q = Y(o+3, :); v = Y(o+5, :); w = Y(o+7, :); p = Y(o+9, :);
    Ru = (th.alpha*u.^2 + 1).*v - th.beta*u - th.gamma*u.*p;   % UNCLIPPED
    Rq =  th.gamma*u.*p - th.delta*q;
    Rv = -(th.alpha*u.^2 + 1).*v + th.beta*u + th.epsilon*w;
    Rw =  th.delta*q - th.epsilon*w;
    Rp = -th.gamma*u.*p + th.delta*q;
    F(o+1, :)  = Y(o+2, :);   F(o+2, :)  = -(Ru + sG)/mu(1);
    F(o+3, :)  = Y(o+4, :);   F(o+4, :)  = -Rq/mu(2);
    F(o+5, :)  = Y(o+6, :);   F(o+6, :)  = -Rv/mu(3);
    F(o+7, :)  = Y(o+8, :);   F(o+8, :)  = -Rw/mu(4);
    F(o+9, :)  = Y(o+10, :);  F(o+10, :) = -(Rp + sP)/mu(5);
end
F(21:24, :) = 0;
F(25, :) = Y(1, :) + Y(3, :) + Y(5, :) + Y(7, :);      % I_GL'
F(26, :) = Y(9, :) + Y(3, :);                          % I_PL'
F(27, :) = Y(11, :) + Y(13, :) + Y(15, :) + Y(17, :);  % I_GR'
F(28, :) = Y(19, :) + Y(13, :);                        % I_PR'
end

function s = shifts(sol)
s = struct('names', {{'s_GL','s_PL','s_GR','s_PR'}}, 'values', sol.y(21:24, 1).');
s.max = max(abs(s.values));
end

function [lab, info] = shape_label(sol, theta, mo, fs) %#ok<INUSL>
% The shape label: the mode name (sa_mode_name) of the main
% variable (model_options.main_variable, default u; a sum of fields is summed
% per half, in the given order) in each half, with the flatness reference
% ref = the main variable's value in the uniform (HSS) state at the anchor (fs:
% field_scale's struct, its anchor_hss_values: the library HSS interpolated
% to the anchor). shape_nodes points per half (x = 0 the pole,
% x = 1/2 the middle). It never affects validity or root identity. It sets the
% attempt's observed_family and, through the per-root pass, that pass's
% unit family string and the child's
% seed_family, so it is part of the run's settings; it never chooses a
% starting profile, a node, an order or a cap. info: sa_mode_name's info, plus .legacy, the
% older three-state label (HSS / HPSS / MPSS per half: flat when
% max - min < 0.25, a maximum at either end of the half -> HPSS).
xi = linspace(sol.x(1), sol.x(end), mo.shape_nodes);
Y  = deval(sol, xi);
v = sa_main_variable('parse', mo.main_variable, species());
[lab, info] = label_pair(sa_main_variable('sum', Y, 2*v.index - 1), ...
    sa_main_variable('sum', Y, 2*v.index + 9), mo, sa_mode_name('reference', fs, v.index));   % rows 1 3 .. 9 / 11 13 .. 19
end

function [lab, info] = shape_label_profile(xi, P, theta, mo, fs) %#ok<INUSL>
% The same rule on a stored profile instead of a solver solution: P holds
% the identity rows (left u q v w p, then right) on the mesh xi, e.g. a
% root's representative profile on the comparison mesh. Each field the
% main variable needs is interpolated linearly to shape_label's
% shape_nodes points, then summed in the same order. For post-processing
% only (a plot relabelled with another main variable); a run's labels come
% from shape_label.
xs = linspace(xi(1), xi(end), mo.shape_nodes);
v = sa_main_variable('parse', mo.main_variable, species());
Y = zeros(10, numel(xs));
for r = [v.index, v.index + 5]
    Y(r, :) = interp1(xi, P(r, :), xs, 'linear');
end
[lab, info] = label_pair(sa_main_variable('sum', Y, v.index), sa_main_variable('sum', Y, v.index + 5), mo, ...
    sa_mode_name('reference', fs, v.index));
end

function s = species()
s = {'u', 'q', 'v', 'w', 'p'};
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

function s = legacy_shape(u, mo)
% One half's older label (tests only): HSS when
% max - min < 0.25 (an absolute amplitude), HPSS when the
% maximum lies within shape_edge_frac of either end, otherwise MPSS.
if max(u) - min(u) < 0.25, s = 'HSS'; return; end
n = numel(u); e = max(1, round(mo.shape_edge_frac * n));
[~, k] = max(u);
if k <= e || k >= n - e + 1, s = 'HPSS'; else, s = 'MPSS'; end
end

function fs = field_scale(theta, lib, S)
% Per field max(|value| of the HSS root chosen by
% the hss_branch rule at the anchor, floor). One value per species, applied
% to both halves. hss_values: the same HSS values, raw; anchor_hss_values:
% the library HSS interpolated to the anchor itself (the uniform-state
% reference of the shape label, sa_mode_name 'reference').
lib = sap_library('profiles', lib, {'HSS'});
c = sap_library('cell', lib, theta, S.lattice.anchor);
if ~c.in_range
    error('SAP:FIELD_SCALE', 'The anchor lies outside the library''s native range.');
end
h = sap_library('half', lib, theta, 'HSS', c.cell, S.model_options);
if ~h.available
    error('SAP:FIELD_SCALE', ['No HSS half exists at the anchor cell (%d,%d): the ' ...
        'field-scale rule has no root to use.'], c.cell(1), c.cell(2));
end
v5 = max(abs(h.Y5(:, 1)), S.identity.field_scale_floor);
fs = struct();
fs.values  = [v5; v5];
fs.species_values = v5.';
fs.species = {'u','q','v','w','p'};
fs.source  = sprintf(['HSS half, hss_branch = %s, native cell (%d,%d) at anchor ' ...
    '[%s], normalized units, floor %g'], S.model_options.hss_branch, ...
    c.cell(1), c.cell(2), num2str(S.lattice.anchor), S.identity.field_scale_floor);
fs.hss_values = h.Y5(:, 1).';
% the uniform-state reference of the shape label: the library HSS interpolated
% bilinearly to the anchor itself (the native grid has no node there in
% general); the nearest cell's values when a surrounding cell has none
val = @(i, j) hss_at(lib, theta, [i j], S.model_options);
[v, why] = sa_mode_name('bilinear', lib.axis_row, lib.axis_col, S.lattice.anchor * theta.M_BG, val);
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
h = sap_library('half', lib, theta, 'HSS', cell, mo);
if h.available, y = h.Y5(:, 1); else, y = []; end
end

function d = option_defaults()
% Model-owned settings (documented defaults). The validity tolerances are the
% values of the coupled solve, which receives the same numbers.
d = struct();
d.hss_branch        = 'representative';  % HSS seed = the library's uniform state
d.guess_nodes       = 401;     % cold-guess mesh in x
d.seam_tol          = 1e-6;    % middle concentration / middle flux
d.g6a_tol           = 1e-8;    % accumulator mass
d.g6b_tol           = 1e-7;    % dense mass
d.g6b_refine_tol    = 1e-9;    % dense mass refinement 20001 vs 40001
d.antisym_abs       = 1e-9;    % gate = max(abs, rel * s_max)
d.antisym_rel       = 1e-6;
d.vel_abs           = 1e-9;    % gate = max(abs, rel * |velocity|)
d.vel_rel           = 1e-6;
d.main_variable     = 'u';     % the label's and dn2's variable: a '+'-sum of u q v w p
d.shape_flat_tol    = 0.5606;  % label: a half is 0 (flat) when (max - min) / (uniform value) is below;
                               %   0.5606 = 0.25 / 0.445953, where 0.445953 is u of the uniform state at the anchor (L = 7)
d.shape_turn_frac   = 0.2;     % label: a turn counts when the profile moves back by this x (max - min)
d.shape_edge_frac   = 0.02;    % label: the window at the middle, as a fraction of a half
d.shape_nodes       = 401;
end

function mo = validate_options(mo)
id = 'SA:SETTINGS:E_FIELD';
v = sa_main_variable('parse', mo.main_variable, species());
mo.main_variable = v.text;               % canonical form (no whitespace)
if ~strcmp(mo.hss_branch, 'representative')
    error(id, 'model_options.hss_branch must be ''representative'' (the only rule implemented).');
end
for f = {'guess_nodes','shape_nodes'}
    v = mo.(f{1});
    if ~(isa(v,'double') && isscalar(v) && v >= 3 && v == fix(v))
        error(id, 'model_options.%s must be an integer >= 3.', f{1});
    end
end
for f = {'seam_tol','g6a_tol','g6b_tol','g6b_refine_tol','antisym_abs','antisym_rel', ...
         'vel_abs','vel_rel','shape_flat_tol','shape_edge_frac'}
    v = mo.(f{1});
    if ~(isa(v,'double') && isscalar(v) && isfinite(v) && v > 0)
        error(id, 'model_options.%s must be a finite positive scalar.', f{1});
    end
end
v = mo.shape_turn_frac;
if ~(isa(v,'double') && isscalar(v) && v > 0 && v < 1)
    error(id, 'model_options.shape_turn_frac must be a scalar in (0, 1).');
end
end

function c = copy_records(here)
% Provenance of the verbatim copies: copy SHA-256 only;
% the source path and source SHA-256 are recorded and checked by the tests.
names = {'sap_p28_solve.m','sap_forced28.m','sap_bvp_convergence.m','sap_model.m'};
c = struct('file', {}, 'copy_sha256', {});
for k = 1:numel(names)
    c(end+1) = struct('file', names{k}, ...
        'copy_sha256', sa_sha256(fullfile(here, names{k}))); %#ok<AGROW>
end
end
