function rows = sa_t71_traces(E)
%SA_T71_TRACES  Full-PDE traces, with a stub
%   integrator on stub runs (no time stepping), and the real post modules'
%   IC builder and stop rule on synthetic states.
T = 'T71 traces';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d0 = sa_m3_testdir(E, 't71');
d = fullfile(d0, 'stub');
sa_m3_run(d, [1 1]);
R = seeded_atlas_roots(d, struct('verify', false));
cnt = containers.Map({'n'}, {0});
drift = containers.Map({'v'}, {0});
integ = @(theta, Y0, entry) stub_integrate(theta, Y0, entry, cnt, drift);
o = struct('integrator', integ);
eA = struct('name', 'A', 'root', R(1).root_id);

% ---- schema, reuse, regeneration ----------------------------------------------------
a1 = seeded_atlas_trace(d, eA, o);
tf = fullfile(d, 'traces', 'A', 'trace.mat');
Z = load(tf);
need = {'t', 'coord', 'coord_labels', 'masses', 'mass_drift_rel', 'umax', 'rhs_inf', 'profiles', 'profile_frames', ...
    'x', 'Y0', 'Y_root', 'Y_final', 'stop', 'stop_detail', 'wall_s', 'entry', 'key', 'root_id', 'plotted_coord', ...
    'frame0_gap', 'provenance'};
rows = sa_check(rows, T, 'trace.mat schema; profile frames include the first and the last; summary and figure written', ...
    a1.ran && all(isfield(Z, need)) && Z.profile_frames(1) == 1 && Z.profile_frames(end) == numel(Z.t) && ...
    size(Z.profiles, 3) == numel(Z.profile_frames) && exist(fullfile(d, 'traces', 'A', 'trace_summary.txt'), 'file') == 2 && ...
    exist(fullfile(d, 'traces', 'A', 'roots_3d_trace.fig'), 'file') == 2, strjoin(setdiff(need, fieldnames(Z)), ','));
n0 = cnt('n');
a2 = seeded_atlas_trace(d, eA, o);
rows = sa_check(rows, T, 'an unchanged entry is skipped (zero integrator calls)', ~a2.ran && cnt('n') == n0, a2.reason);
eA2 = eA; eA2.t_out = 0:2:30;
a3 = seeded_atlas_trace(d, eA2, o);
rows = sa_check(rows, T, 'a changed t_out regenerates and keeps trace_previous.mat', a3.ran && cnt('n') == n0 + 1 && ...
    exist(fullfile(d, 'traces', 'A', 'trace_previous.mat'), 'file') == 2 && contains(a3.reason, 'entry'), a3.reason);
eB = struct('name', 'B', 'root', R(end).root_id);
n1 = cnt('n');
b1 = seeded_atlas_trace(d, eB, o);
a4 = seeded_atlas_trace(d, eA2, o);
rows = sa_check(rows, T, 'a new entry runs alone', b1.ran && ~a4.ran && cnt('n') == n1 + 1);

% ---- comment-only vs code edits (the key's code part) -----------------------------------------------
sd = fullfile(d0, 'digest'); mkdir(fullfile(sd, 'c')); mkdir(fullfile(sd, 'k'));
src = fileread(fullfile(E.pkg, 'seeded_atlas_trace.m'));
wr(fullfile(sd, 'c', 'seeded_atlas_trace.m'), strrep(src, 'if nargin < 3, opts = struct(); end', ...
    sprintf('%% a planted comment\nif nargin < 3, opts = struct(); end   %% planted trailing comment')));
wr(fullfile(sd, 'k', 'seeded_atlas_trace.m'), strrep(src, 'if nargin < 3, opts = struct(); end', ...
    'if nargin < 3, opts = struct(); end; opts = opts;'));
h0 = sa_code_digest('files', {fullfile(E.pkg, 'seeded_atlas_trace.m')});
hc = sa_code_digest('files', {fullfile(sd, 'c', 'seeded_atlas_trace.m')});
hk = sa_code_digest('files', {fullfile(sd, 'k', 'seeded_atlas_trace.m')});
rows = sa_check(rows, T, 'the key''s code digest ignores a planted comment edit and detects a planted code edit', ...
    strcmp(h0, hc) && ~strcmp(h0, hk), sprintf('%s / %s / %s', h0(1:12), hc(1:12), hk(1:12)));
Zb = load(fullfile(d, 'traces', 'B', 'trace.mat'));
Zb.provenance.files(1).sha256 = repmat('0', 1, 64);           % raw SHA-256: provenance only
sa_write('mat', fullfile(d, 'traces', 'B', 'trace.mat'), Zb);
n2 = cnt('n');
b2 = seeded_atlas_trace(d, eB, o);
rows = sa_check(rows, T, 'a changed raw file SHA-256 (provenance) does not regenerate', ~b2.ran && cnt('n') == n2, b2.reason);
Zb.key.code = 'planted-code-digest';
sa_write('mat', fullfile(d, 'traces', 'B', 'trace.mat'), Zb);
b3 = seeded_atlas_trace(d, eB, o);
rows = sa_check(rows, T, 'a changed code digest regenerates, keeping trace_previous.mat', b3.ran && cnt('n') == n2 + 1 && ...
    exist(fullfile(d, 'traces', 'B', 'trace_previous.mat'), 'file') == 2 && contains(b3.reason, 'code'), b3.reason);

% ---- figure-only redraws -----------------------------------------------------------------------------
kf = fullfile(d, 'traces', 'B', 'roots_3d_trace_key.mat');
K0 = load(kf); K1 = K0; K1.key.code = 'planted-plot-digest';
sa_write('mat', kf, K1);
n3 = cnt('n');
seeded_atlas_trace(d, eB, o);
K2 = load(kf);
rows = sa_check(rows, T, 'a plot-code digest change redraws only the figure (zero integration)', cnt('n') == n3 && ...
    isequal(K2.key, K0.key));
db = fullfile(d0, 'budget');
sa_m3_run(db, [1 1], 'settings', struct('budget', struct('max_calls', 6)));
Rb = seeded_atlas_roots(db, struct('verify', false));
eC = struct('name', 'C', 'root', Rb(1).root_id);
seeded_atlas_trace(db, eC, o);
kc = fullfile(db, 'traces', 'C', 'roots_3d_trace_key.mat');
C0 = load(kc);
sa_m3_run(db, [1 1], 'settings', struct('budget', struct('max_calls', 5000), 'resume', true));
n4 = cnt('n');
c2 = seeded_atlas_trace(db, eC, o);
C1 = load(kc);
Rb2 = seeded_atlas_roots(db, struct('verify', false));
rows = sa_check(rows, T, 'a resume that adds roots redraws only the figure (zero integration, the root set changed)', ...
    ~c2.ran && cnt('n') == n4 && numel(Rb2) > numel(Rb) && ~isequal(C0.key.roots, C1.key.roots), ...
    sprintf('roots %d -> %d', numel(Rb), numel(Rb2)));

% ---- the IC comes from the representative attempt's sol ----------------------------------------------
P = sa_post_model('stub_two_pool');
Ef = load(fullfile(sa_layout('state', d), 'settings_effective.mat'), 'theta');
r = R(end);
S = load(fullfile(sa_layout('state', d), 'attempts', [r.representative_attempt_id '.mat']), 'sol');
icr = P.trace_ic(Ef.theta, S.sol, Zb.entry);
Zb2 = load(fullfile(d, 'traces', 'B', 'trace.mat'));
Fp = P.full_profile(Ef.theta, r.comparison_mesh, r.representative_profile(P.profile_rows(1, :), :), ...
    r.representative_profile(P.profile_rows(2, :), :), Zb2.entry.nx);
dm = fullfile(d0, 'planted_sol'); copyfile(d, dm);
af = fullfile(dm, 'internal', 'attempts', [r.representative_attempt_id '.mat']);
Za = load(af); Za.sol.y(1, :) = 1.01 * Za.sol.y(1, :); save(af, '-struct', 'Za');
pm = seeded_atlas_trace(dm, eB, o);
Zm = load(fullfile(dm, 'traces', 'B', 'trace.mat'));
icm = P.trace_ic(Ef.theta, Za.sol, Zm.entry);
rows = sa_check(rows, T, ['the IC is the representative attempt''s sol on its own mesh: equal to trace_ic(sol) bit for ' ...
    'bit; a planted change of that sol (comparison profile untouched) moves the IC and regenerates the trace'], ...
    isequal(Zb2.Y_root, icr.Y_root) && pm.ran && isequal(Zm.Y_root, icm.Y_root) && ~isequal(Zm.Y_root, Zb2.Y_root), ...
    sprintf('|sol route - profile route| %.3g; planted reason: %s', max(abs(Zb2.Y_root(:) - Fp.Y(:))), pm.reason));

% ---- stops per model ------------------------------------------------------------------------------------
Pf = sa_post_model('original_polarity');
[k1, ~] = Pf.stop_rule([1 + 2e-6, 0.2], [1 0.2], 1);
[k2, ~] = Pf.stop_rule([1 + 5e-7, 0.2], [1 0.2], 1e-9);
[k3, ~] = Pf.stop_rule([1, 0.2 * (1 + 1.1e-6)], [1 0.2], 1);
[k4, ~] = Pf.stop_rule([1 + 5e-7, 0.2], [1 0.2], 1);
rows = sa_check(rows, T, 'fission yeast: mass_drift_exceeded above 1e-6 (either pool), stationary at max|dY/dt| <= 1e-8', ...
    strcmp(k1, 'mass_drift_exceeded') && strcmp(k2, 'stationary') && strcmp(k3, 'mass_drift_exceeded') && isempty(k4));
Pn = sa_post_model('guzzo_noxa');
rows = sa_check(rows, T, 'noxa: no drift or stationarity stop (stops horizon_reached / pdepe_incomplete), WARNING above 1e-4', ...
    isequal(Pn.stops, {'horizon_reached', 'pdepe_incomplete'}) && Pn.drift_warn == 1e-4 && isinf(Pf.drift_warn));
drift('v') = 2e-4;
seeded_atlas_trace(d, struct('name', 'W', 'root', R(1).root_id), o);
drift('v') = 5e-5;
seeded_atlas_trace(d, struct('name', 'W2', 'root', R(1).root_id), o);
drift('v') = 0;
s1 = fileread(fullfile(d, 'traces', 'W', 'trace_summary.txt'));
s2 = fileread(fullfile(d, 'traces', 'W2', 'trace_summary.txt'));
rows = sa_check(rows, T, 'trace_summary.txt warns when the drift exceeds the post module''s warning level', ...
    contains(s1, 'WARNING') && ~contains(s2, 'WARNING'));
rows = sa_check(rows, T, 'the frame-0 gap is reported (trace.mat and trace_summary.txt)', ...
    all(isfinite(Z.frame0_gap)) && contains(s1, 'frame-0 gap'), mat2str(Z.frame0_gap, 3));

% ---- cos perturbation: pools unchanged at each node, positivity kept --------------------------------------
[okf, df] = cos_check(Pf, 'fy');
[okn, dn] = cos_check(Pn, 'noxa');
rows = sa_check(rows, T, 'cos perturbation (a = 0.9): every pool unchanged at each node (<= 4 ulp), positivity kept', ...
    okf && okn, sprintf('%s; %s', df, dn));

% ---- flat folder via opts.out; atomic write --------------------------------------------------------------------
fl = fullfile(d0, 'flat');
mkdir(fl);
copyfile(fullfile(sa_layout('state', d), '*'), fl);
for f = sa_layout('human_files'), if exist(fullfile(d, f{1}), 'file') == 2, copyfile(fullfile(d, f{1}), fl); end, end
id = '';
try, seeded_atlas_trace(fl, eA, o); catch ME, id = ME.identifier; end
fo = seeded_atlas_trace(fl, eA, struct('integrator', integ, 'out', fullfile(d0, 'flat_out')));
rows = sa_check(rows, T, 'a flat folder: refused without opts.out (E_LEGACY), traced into opts.out', ...
    strcmp(id, 'SA:TRACE:E_LEGACY') && fo.ran && exist(fullfile(d0, 'flat_out', 'traces', 'A', 'trace.mat'), 'file') == 2 && ...
    strcmp(sa_layout('layout', fl), 'flat') && exist(fullfile(fl, 'traces'), 'dir') ~= 7, id);
bad = @(theta, Y0, entry) error('SA:TEST:KILL', 'planted kill before the save');
try, seeded_atlas_trace(d, struct('name', 'K', 'root', R(1).root_id), struct('integrator', bad)); catch, end
rows = sa_check(rows, T, 'a planted kill before the save leaves no trace.mat', ...
    exist(fullfile(d, 'traces', 'K', 'trace.mat'), 'file') ~= 2);
end

% =========================================================================

function fr = stub_integrate(theta, Y0, entry, cnt, drift)
cnt('n') = cnt('n') + 1;
P = sa_post_model('stub_two_pool');
t = entry.t_out; F = numel(t);
x = linspace(0, 2 * theta.Lh, size(Y0, 2));
O = P.observables(theta, x, Y0);
fr = struct('t', t, 'coord', repmat(O.c, F, 1) + (0:F-1).' * [1e-3 0 0], ...
    'masses', repmat(O.masses, F, 1) .* (1 + (0:F-1).' / max(F - 1, 1) * drift('v')), ...
    'umax', repmat(O.umax, F, 1), 'rhs_inf', nan(1, F), 'profile_frames', [], 'stop', 'horizon_reached', ...
    'stop_detail', '', 'wall_s', 0, 'x', x);
pf = unique([1:entry.profile_every:F, F]);
fr.profile_frames = pf;
fr.profiles = repmat(Y0, 1, 1, numel(pf));
fr.mass_drift_rel = abs(fr.masses - fr.masses(1, :)) ./ fr.masses(1, :);
fr.Y_final = Y0;
end

function [ok, txt] = cos_check(P, which)
switch which
    case 'fy'
        x = linspace(0, 0.5, 61); y = zeros(28, numel(x));
        base = [1.2 0.3 0.8 0.5 0.4];
        for k = 1:5
            y(P.sol_rows(1, k), :) = base(k) * (1 + 0.3 * cos(2 * pi * x));
            y(P.sol_rows(2, k), :) = base(k) * (1 - 0.2 * cos(2 * pi * x));
        end
        th = struct('L', 1);
        pools = {[1 2 3 4], [2 5]};
    case 'noxa'
        th = struct('L', 6); x = linspace(0, 3, 61); y = zeros(24, numel(x));
        base = [0.4 0.6 0.7 0.3];
        for k = 1:4
            y(P.sol_rows(1, k), :) = base(k) * (1 + 0.3 * cos(pi * x / 3));
            y(P.sol_rows(2, k), :) = base(k) * (1 - 0.2 * cos(pi * x / 3));
        end
        pools = {[1 2], [3 4]};
end
sol = struct('x', x, 'y', y);
e = struct('nx', 101, 'perturbation', struct('kind', 'cos', 'amplitude', 0.9, 'mode', 1));
ic = P.trace_ic(th, sol, e);
ok = all(ic.Y(:) > 0) && ~isequal(ic.Y, ic.Y_root);
worst = 0;
for k = 1:numel(pools)
    a = sum(ic.Y(pools{k}, :), 1); b = sum(ic.Y_root(pools{k}, :), 1);
    worst = max(worst, max(abs(a - b) ./ eps(b)));
end
ok = ok && worst <= 4;
txt = sprintf('%s: worst pool change %g ulp, min %.3g', which, worst, min(ic.Y(:)));
end

function wr(f, t)
fid = fopen(f, 'w'); fwrite(fid, t, 'char'); fclose(fid);
end
