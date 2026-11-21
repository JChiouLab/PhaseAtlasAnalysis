function rows = sa_t68_lattice_axes(E)
%SA_T68_LATTICE_AXES  lattice.points_per_axis.
%   Zero solves (stub corrector on the stub model; the real modules are only
%   resolved and their admissibility predicate evaluated).
T = 'T68 per-axis lattice';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d0 = sa_m3_testdir(E, 't68');

% ---- widths by bisection: stub and both real modules ------------------------
Ms = sa_model('stub_two_pool');
ths = struct('R', 4, 'B', 4, 'Lh', 1.5);
S = res(Ms, struct('lattice', struct('points_per_axis', [7 9])));
S2 = sa_lattice('derive_step', Ms, ths, S);
Rs = 4 - 2 * 1e-4;                                     % (R - d)/2 >= margin
rows = sa_check(rows, T, 'stub: open interval, step_c = 2 R_c / n_c (R = 3.9998 both axes)', ...
    max(abs(S2.lattice.step - 2 * Rs ./ [7 9])) <= 1e-12 * Rs, mat2str(S2.lattice.step, 17));
Mf = sa_model('original_polarity');
thf = Mf.normalize(E.p14);
Sf = sa_lattice('derive_step', Mf, thf, res(Mf, struct('lattice', struct('points_per_axis', [25 99])), 'original_polarity'));
Rf = [thf.rho 1] - 2 * 1e-4;
rows = sa_check(rows, T, 'fission yeast L = 14: widths rho - 2 margin and 1 - 2 margin; steps 2R/n', ...
    max(abs(Sf.lattice.step - 2 * Rf ./ [25 99]) ./ (2 * Rf ./ [25 99])) <= 1e-12, ...
    sprintf('steps %s (expected %s)', mat2str(Sf.lattice.step, 8), mat2str(2 * Rf ./ [25 99], 8)));
Mn = sa_model('guzzo_noxa');
thn = Mn.normalize(sag_test_fixture('params', 1));
Sn = sa_lattice('derive_step', Mn, thn, res(Mn, struct('lattice', struct('points_per_axis', [27 27])), 'guzzo_noxa'));
Rn = 2 * (1 - 1e-4 / thn.Lh);
rows = sa_check(rows, T, 'noxa f = 1: width 2 (rho - margin / Lh); [27 27] gives step 0.148143', ...
    max(abs(Sn.lattice.step - 2 * Rn / 27)) <= 1e-12 && abs(Sn.lattice.step(1) - 0.148143) < 1e-6, ...
    mat2str(Sn.lattice.step, 10));
Lf = sa_lattice('build', Mf, thf, Sf);
rows = sa_check(rows, T, 'fission yeast [25 99]: 25 x 99 nodes, outermost node half a step inside', ...
    isequal(Lf.index_max - Lf.index_min + 1, [25 99]) && ...
    all(abs(max(Lf.coordinate) - (Rf - Sf.lattice.step / 2)) <= 1e-12), ...
    sprintf('%d x %d', Lf.index_max - Lf.index_min + 1));

% ---- refusals ------------------------------------------------------------------
e1 = errmsg(@() res(Ms, struct('lattice', struct('points_per_axis', 50))));
rows = sa_check(rows, T, 'an even count is refused, naming the two nearest odd counts (50 -> 49 / 51)', ...
    contains(e1, '50 is even') && contains(e1, 'use 49 or 51') && contains(e1, 'Unlike res, the half-cell library grid'), e1);
e2 = errmsg(@() res(Ms, struct('lattice', struct('points_per_axis', [1 5]))));
rows = sa_check(rows, T, 'n < 3 is refused', contains(e2, 'below 3'), e2);
[e3, i3] = errmsg(@() res(Ms, struct('lattice', struct('points_per_axis', [7 9], 'step', 0.02))));
rows = sa_check(rows, T, 'both set by the user is refused (by source): E_CONTRADICTION', ...
    strcmp(i3, 'SA:SETTINGS:E_CONTRADICTION') && contains(e3, 'not both'), e3);
[~, i4] = errmsg(@() res(Mf, struct('lattice', struct('points_per_axis', [25 99])), 'original_polarity'));
rows = sa_check(rows, T, 'a model-layer step plus a user count is accepted', isempty(i4), i4);
[e5, i5] = errmsg(@() res(Ms, struct('lattice', struct('points_per_axis', [7 9], 'anchor', [0.65 0]))));
rows = sa_check(rows, T, 'a moved anchor is refused', strcmp(i5, 'SA:SETTINGS:E_CONTRADICTION') && ...
    contains(e5, 'anchor'), e5);
Ma = sa_model('stub_asym');
[e6, i6] = errmsg(@() sa_lattice('derive_step', Ma, ths, res(Ma, struct('lattice', ...
    struct('points_per_axis', [7 9])), 'stub_asym')));
rows = sa_check(rows, T, 'asymmetric pools are refused (SA:LATTICE:E_STEP)', strcmp(i6, 'SA:LATTICE:E_STEP'), e6);
da = fullfile(d0, 'asym_run');
sa_ = sa_m3_stub('settings', da, [1 1]); sa_.model = 'stub_asym';
sa_.lattice = struct('points_per_axis', [7 9]);
[~, i7] = errmsg(@() seeded_atlas_run(sa_m3_stub('params'), struct(), sa_, ...
    @(varargin) error('T68:CALLED', 'corrector called')));
rows = sa_check(rows, T, 'a refused derivation inside seeded_atlas_run creates no folder', ...
    strcmp(i7, 'SA:LATTICE:E_STEP') && exist(da, 'dir') ~= 7, i7);

% ---- cap rule: closed interval, the end node on the cap --------------------------------
Sc = sa_lattice('derive_step', Ms, ths, res(Ms, struct('lattice', struct('points_per_axis', [5 7], ...
    'max_abs_coordinate', [1 2]))));
Lc = sa_modes('cap_lattice', sa_lattice('build', Ms, ths, Sc), Sc);
cmax = max(abs(Lc.coordinate));
rows = sa_check(rows, T, 'cap: step 2 cap / (n - 1); the capped lattice has n nodes per axis; end nodes on the cap', ...
    max(abs(Sc.lattice.step - [2/4 4/6])) == 0 && isequal(Lc.index_max - Lc.index_min + 1, [5 7]) && ...
    all(abs(cmax - [1 2]) <= 1e-12 * [1 2]), sprintf('steps %s; max |d| %s', mat2str(Sc.lattice.step, 17), ...
    mat2str(cmax, 17)));

% ---- the walk count at a large odd n -------------------------------------------------------
S9 = sa_lattice('derive_step', Mf, thf, res(Mf, struct('lattice', struct('points_per_axis', [999 999])), ...
    'original_polarity'));
cnt = zeros(1, 2);
for c = 1:2
    kp = 0; while sa_lattice('admissible', Mf, thf, S9, unit(c) * S9.lattice.step(c) * (kp + 1)), kp = kp + 1; end
    km = 0; while sa_lattice('admissible', Mf, thf, S9, -unit(c) * S9.lattice.step(c) * (km + 1)), km = km + 1; end
    cnt(c) = kp + km + 1;
end
rows = sa_check(rows, T, 'the walk rule gives exactly 999 nodes per axis at n = 999 (fission yeast L = 14)', ...
    isequal(cnt, [999 999]), mat2str(cnt));

% ---- scalar 0.02 == row [0.02 0.02] bit for bit -----------------------------------------------
Sa = res(Mf, struct(), 'original_polarity');
Sb = res(Mf, struct('lattice', struct('step', [0.02 0.02])), 'original_polarity');
La = sa_lattice('build', Mf, thf, Sa); Lb = sa_lattice('build', Mf, thf, Sb);
same = isequal(La.coordinate, Lb.coordinate) && isequal(La.key, Lb.key) && isequal(La.index, Lb.index) && ...
    bits(La.coordinate, Lb.coordinate);
Sa.lattice.max_abs_coordinate = 0.06; Sb.lattice.max_abs_coordinate = 0.06;
Ca = sa_modes('cap_lattice', La, Sa); Cb = sa_modes('cap_lattice', Lb, Sb);
same_cap = isequal(Ca.key, Cb.key) && bits(Ca.coordinate, Cb.coordinate);
C = struct('in_range', true(La.n_nodes, 1));
Sa.pilot = struct('sample_per_axis', 7); Sb.pilot = Sa.pilot;
same_pilot = isequal(sa_modes('pilot_sample', La, C, Sa), sa_modes('pilot_sample', Lb, C, Sb));
rng(68); v = rand(200, 2) * 20 - 10;
fa = La.anchor + La.step * v(1, :); fb = Lb.anchor + sa_lattice('steps', Lb) .* v(1, :);
ok_front = bits(fa, fb);
for k = 1:size(v, 1)
    ok_front = ok_front && bits(Sa.lattice.anchor + Sa.lattice.step * v(k, :), ...
        Sb.lattice.anchor + Sb.lattice.step .* v(k, :));
end
tola = 1e-9 * La.step; tolb = 1e-9 * sa_lattice('steps', Lb);
ok_match = isequal(abs(La.coordinate - La.coordinate(5, :)) <= tola, abs(Lb.coordinate - Lb.coordinate(5, :)) <= tolb);
rows = sa_check(rows, T, ['row step [0.02 0.02] = scalar 0.02 bit for bit: coordinates, keys, cap, pilot sample, ' ...
    'frontier and edge coordinates (anchor + step .* idx), compare matching'], ...
    same && same_cap && same_pilot && ok_front && ok_match, sprintf('%d %d %d %d %d', same, same_cap, ...
    same_pilot, ok_front, ok_match));

% ---- a stub run with points_per_axis ------------------------------------------------------------
d = fullfile(d0, 'ppa_run');
s = sa_m3_stub('settings', d, [1 1]);
s.lattice = struct('points_per_axis', [7 9]);
st = sa_m3_stub('state');
r = seeded_atlas_run(sa_m3_stub('params'), struct(), s, sa_m3_stub('corrector', st));
Ez = load(fullfile(sa_layout('state', d), 'settings_effective.mat'));
Qf = load(fullfile(sa_layout('state', d), 'queue_final.mat'));
h = Ez.lattice.steps; anc = Ez.lattice.anchor;
A = r.attempts; okc = true; nsub = 0; nnode = 0;
for k = 1:numel(A)
    a = A(k);
    if strcmp(a.outcome, 'interrupted'), continue; end
    if a.fraction == 1 || strcmp(a.kind, 'cold') || strcmp(a.kind, 'refinement')
        okc = okc && bits(a.coordinate, anc + h .* a.index); nnode = nnode + 1;
    else
        e = Qf.Q.edges(a.edge_id);
        kp = Qf.Q.index(e.p, :); kq = Qf.Q.index(e.q, :);
        okc = okc && bits(a.coordinate, anc + h .* (kp + a.fraction * (kq - kp))); nsub = nsub + 1;
    end
end
rows = sa_check(rows, T, ['a stub run with points_per_axis [7 9] completes; node and substep coordinates equal ' ...
    'anchor + steps .* idx bit for bit'], strcmp(r.stop_reason, 'completed') && okc && ...
    isequal(Ez.lattice.index_max - Ez.lattice.index_min + 1, [7 9]), sprintf('%d node rows, %d substep rows, stop %s', ...
    nnode, nsub, r.stop_reason));
txt = fileread(fullfile(d, 'settings_effective.txt'));
ln = regexp(txt, '(?m)^lattice\.step\s+.*$', 'match', 'once');
rows = sa_check(rows, T, 'settings_effective.txt marks the step derived', contains(ln, 'derived'), strtrim(ln));
s2 = s; s2.resume = true; s2.lattice.points_per_axis = [7 11];
[~, i8] = errmsg(@() seeded_atlas_run(sa_m3_stub('params'), struct(), s2, sa_m3_stub('corrector', st)));
rows = sa_check(rows, T, 'resume refuses a changed count (E_SETTINGS)', strcmp(i8, 'SA:RESUME:E_SETTINGS'), i8);
U = sa_settings('from_run', d, 'smoke');
rows = sa_check(rows, T, 'from_run keeps points_per_axis', isequal(U.lattice.points_per_axis, [7 9]), ...
    mat2str(U.lattice.points_per_axis));
end

% =========================================================================

function S = res(M, extra, model)
if nargin < 3, model = 'stub_two_pool'; end
u = struct('model', model, 'run_dir', 'x');
f = fieldnames(extra);
for k = 1:numel(f), u.(f{k}) = extra.(f{k}); end
S = sa_settings('resolve', u, M);
end

function [m, id] = errmsg(fn)
m = ''; id = '';
try
    fn();
catch ME
    m = ME.message; id = ME.identifier;
end
end

function u = unit(c)
u = zeros(1, 2); u(c) = 1;
end

function tf = bits(a, b)
tf = isequal(size(a), size(b)) && isequal(typecast(double(a(:)), 'uint64'), typecast(double(b(:)), 'uint64'));
end
