function out = seeded_atlas_trace(run_dir, entry, opts)
%SEEDED_ATLAS_TRACE  A full-PDE trace started from one root of a run folder.
%
%   out = SEEDED_ATLAS_TRACE(run_dir, entry)
%   out = SEEDED_ATLAS_TRACE(run_dir, entry, opts)
%
%   entry (one struct; [] / missing fields take the model's defaults)
%     name           folder traces/<name>/ (letters, digits, _ and -)
%     root           a root_id of roots_index.csv, e.g. 'k:-1,4/R000001'
%     at             instead of root: a position [x y z] in the plotted
%                    coordinates of plots/roots_3d.fig (fission yeast dG, dP,
%                    dn2; Guzzo_noXA dn_R, dn_B, dn2 of the main variable), for
%                    example pasted from the figure's data tip. The trace
%                    starts from the root closest on the figure: each axis
%                    difference divided by that axis's plotted range (max -
%                    min over the run's roots, the coordinates the plot
%                    draws; an axis whose extent on the figure is below 1 %
%                    of the largest axis's extent counts as that 1 %), then
%                    the Euclidean norm; ties go to the first
%                    root in roots_index.csv order. The chosen root, its
%                    coordinates and the distance (scaled, and raw per axis)
%                    are printed and stored (trace.mat at_resolution,
%                    trace_summary.txt). Give root or at, not both.
%     perturbation   struct('kind', 'none') or struct('kind', 'cos',
%                    'amplitude', a, 'mode', m), 0 < a < 1, m a positive
%                    integer: for each mass-conserving field pair (a, b) of
%                    the model, a is raised and b lowered pointwise by
%                    a min(Y_a, Y_b) cos(m pi eta) (every pool unchanged at
%                    each node, positivity kept)
%     t_out          output times
%     nx             uniform full-domain mesh (odd, a node on the seam)
%     profile_every  full profiles every Nth frame (first and last always)
%     reltol, abstol integrator tolerances
%   opts
%     out            folder that receives traces/<name>/; REQUIRED for a
%                    flat (software 0.5) folder, which is never written into
%     integrator     test hook: a function with the signature of the post
%                    module's integrate(theta, Y0, entry)
%     plot           true (default): draw roots_3d_trace.fig
%     main_variable  '' (default: the run's own, as seeded_atlas_plot) or a
%                    '+'-sum of the model's fields: the trace's third
%                    coordinate (dn2), the "at" resolution and the trace
%                    figure use the plot relabelled with it
%                    (seeded_atlas_plot opts.main_variable); it enters
%                    the trace key
%
%   Writes traces/<name>/trace_summary.txt, trace.mat (atomically, after the
%   integration) and roots_3d_trace.fig (seeded_atlas_plot with the trace
%   overlaid). The initial condition is the root's representative attempt
%   sol, assembled on its own bvp4c mesh (post trace_ic); the integrator is
%   the post module's integrate; both are pure functions (no folder, no
%   write).
%
%   Reuse. trace.mat stores a key: the software version; the root id (for
%   an at entry, the root it resolved to on this call), its
%   representative attempt id and that attempt file's SHA-256; the SHA-256
%   of internal/settings_effective.mat; the resolved entry; and the
%   comment-free token digest (sa_code_digest) of this file and of the
%   model's post/ files. An equal key and a complete file: skipped, zero
%   integration. A different key: regenerated, the old file kept as
%   trace_previous.mat, the reason printed (for an at entry whose nearest
%   root changed, for example after a resume added roots: the old and the
%   new root). A comment edit never changes the
%   key; a code edit does. The raw SHA-256 values of those files are stored
%   as provenance only. The figure is redrawn when its own key changed
%   (seeded_atlas_plot: the trace key, the root set, the plot code).
%
%   out: name, dir, file, ran (true when integrated), reason, stop, frames,
%   wall_s, frame0_gap.

if nargin < 3, opts = struct(); end
here = fileparts(mfilename('fullpath'));
addpath_once(fullfile(here, 'post'));
D = sa_layout('state', run_dir);
lay = sa_layout('layout', run_dir);
if strcmp(lay, 'none'), error('SA:TRACE:E_NO_RUN', '%s holds no run.', run_dir); end
base = getf(opts, 'out', '');
if isempty(base)
    if strcmp(lay, 'flat')
        error('SA:TRACE:E_LEGACY', ['%s is a flat (software 0.5 or earlier) run folder: it is evidence and ' ...
            'is never written into; pass opts.out (another folder).'], run_dir);
    end
    base = run_dir;
end
E = load(fullfile(D, 'settings_effective.mat'), 'settings', 'theta');
mvreq = getf(opts, 'main_variable', '');
mv = mvreq;                                      % the coordinates the plot draws (seeded_atlas_plot)
if isempty(mv) && isfield(E.settings, 'model_options') && isfield(E.settings.model_options, 'main_variable')
    mv = E.settings.model_options.main_variable;
end
P = sa_post_model(E.settings.model, mv);
entry = resolve_entry(entry, P);
tdir = fullfile(base, 'traces', entry.name);
tfile = fullfile(tdir, 'trace.mat');

% ---- the root and the key -----------------------------------------------------
R = seeded_atlas_roots(run_dir, struct('verify', false));
at_res = struct();
if isempty(entry.at)
    k = find(strcmp({R.root_id}, entry.root), 1);
    if isempty(k)
        error('SA:TRACE:E_ROOT', 'trace %s: root %s is not a root of %s (see roots_index.csv).', ...
            entry.name, entry.root, run_dir);
    end
else
    [k, at_res] = nearest_root(run_dir, R, entry, mvreq);
end
r = R(k);
af = fullfile(D, 'attempts', [r.representative_attempt_id '.mat']);
code_files = [{fullfile(here, 'seeded_atlas_trace.m')}, P.files];
key = struct('software_version', sa_settings('version'), 'root_id', r.root_id, ...
    'attempt_id', r.representative_attempt_id, 'attempt_sha256', sa_sha256(af), ...
    'settings_sha256', sa_sha256(fullfile(D, 'settings_effective.mat')), 'entry', entry, ...
    'code', sa_code_digest('files', code_files), 'main_variable', '');
if isfield(P, 'main_variable'), key.main_variable = P.main_variable; end
out = struct('name', entry.name, 'dir', tdir, 'file', tfile, 'ran', false, 'reason', '', 'stop', '', ...
    'frames', NaN, 'wall_s', 0, 'frame0_gap', [NaN NaN NaN], 'at_resolution', []);
[need, why] = decide(tfile, key);
if need && ~isempty(entry.at)
    why = at_reason(tfile, key, why);
end
out.reason = why;
out.at_resolution = at_res;
if ~need
    Z = load(tfile, 'stop', 't');
    out.stop = Z.stop; out.frames = numel(Z.t);
    fprintf('trace %s: reused (%s)\n', entry.name, why);
else
    fprintf('trace %s: integrating (%s)\n', entry.name, why);
    Zs = load(af, 'sol');
    ic = P.trace_ic(E.theta, Zs.sol, entry);
    integ = getf(opts, 'integrator', P.integrate);
    fr = integ(E.theta, ic.Y, entry);
    % the plotted coordinates of the root (the plot's route: stored profile)
    Yp = r.representative_profile;
    Fp = P.full_profile(E.theta, r.comparison_mesh, Yp(P.profile_rows(1, :), :), Yp(P.profile_rows(2, :), :), ...
        P.plot_nx);
    Op = P.observables(E.theta, Fp.x, Fp.Y);
    gap = abs(fr.coord(1, :) - Op.c);
    T = struct();
    T.t = fr.t; T.coord = fr.coord; T.coord_labels = {P.coord_names};
    T.masses = fr.masses; T.mass_names = {P.mass_names}; T.mass_drift_rel = fr.mass_drift_rel;
    T.umax = fr.umax; T.umax_names = {P.umax_names}; T.rhs_inf = fr.rhs_inf;
    T.profiles = fr.profiles; T.profile_frames = fr.profile_frames; T.x = fr.x;
    T.Y0 = ic.Y; T.Y_root = ic.Y_root; T.Y_final = fr.Y_final;
    T.stop = fr.stop; T.stop_detail = fr.stop_detail; T.wall_s = fr.wall_s;
    T.entry = entry; T.key = key; T.root_id = r.root_id; T.plotted_coord = Op.c; T.frame0_gap = gap;
    T.at_resolution = at_res;
    T.frame0_gap_note = ['|frame 0 (the attempt sol on its own mesh, nx = entry.nx) - the plotted ' ...
        'coordinates (the stored comparison profile, plot nx)|; with a perturbation it includes it'];
    T.provenance = struct('software_version', sa_settings('version'), 'matlab_version', version, ...
        'files', {cellfun(@(f) struct('file', f, 'sha256', sa_sha256(f)), code_files)}, ...
        'attempt_file', af, 'run_dir', run_dir, 'written_at', datestr(now, 'yyyy-mm-ddTHH:MM:SS')); %#ok<DATST,TNOW1>
    if exist(tdir, 'dir') ~= 7, mkdir(tdir); end
    sa_write('text', fullfile(tdir, 'trace_summary.txt'), summary_text(T, P, run_dir));
    if exist(tfile, 'file') == 2
        movefile(tfile, fullfile(tdir, 'trace_previous.mat'), 'f');
    end
    sa_write('mat', tfile, T);
    out.ran = true; out.stop = T.stop; out.frames = numel(T.t); out.wall_s = T.wall_s; out.frame0_gap = gap;
    fprintf('trace %s: %d frames, stop %s, %.1f s\n', entry.name, numel(T.t), T.stop, T.wall_s);
end
if getf(opts, 'plot', true)
    po = struct('trace', tfile, 'out', tdir, 'file', 'roots_3d_trace', 'main_variable', mvreq);
    seeded_atlas_plot(run_dir, po);
end
end

% =========================================================================

function e = resolve_entry(e, P)
if ~(isstruct(e) && isscalar(e)), error('SA:TRACE:E_ENTRY', 'a trace entry must be a scalar struct.'); end
if ~isfield(e, 'name') || ~ischar(e.name) || isempty(e.name)
    error('SA:TRACE:E_ENTRY', 'a trace entry needs name (text).');
end
has_root = isfield(e, 'root') && ~isempty(e.root);
has_at = isfield(e, 'at') && ~isempty(e.at);
if has_root == has_at
    error('SA:TRACE:E_ENTRY', 'trace %s: give root (a root_id) or at ([x y z]), exactly one of them.', e.name);
end
if has_root && ~ischar(e.root)
    error('SA:TRACE:E_ENTRY', 'trace %s: root must be a root_id (text).', e.name);
end
if has_at
    if ~(isnumeric(e.at) && isreal(e.at) && numel(e.at) == 3 && all(isfinite(e.at(:))))
        error('SA:TRACE:E_ENTRY', 'trace %s: at must be three finite numbers [x y z].', e.name);
    end
    e.at = reshape(double(e.at), 1, 3);
    e.root = '';
else
    e.at = [];
end
if isempty(regexp(e.name, '^[A-Za-z0-9_-]+$', 'once'))
    error('SA:TRACE:E_ENTRY', 'trace name %s: use letters, digits, _ and - only.', e.name);
end
d = P.entry_defaults();
known = [{'name', 'root', 'at'}, fieldnames(d).'];
extra = setdiff(fieldnames(e), known);
if ~isempty(extra)
    error('SA:TRACE:E_ENTRY', 'unknown trace entry field(s): %s.', strjoin(extra(:).', ', '));
end
for f = fieldnames(d).'
    if ~isfield(e, f{1}) || isempty(e.(f{1})), e.(f{1}) = d.(f{1}); end
end
e = orderfields(e, known);
p = e.perturbation;
if ~(isstruct(p) && isfield(p, 'kind') && any(strcmp(p.kind, {'none', 'cos'})))
    error('SA:TRACE:E_ENTRY', 'perturbation.kind must be none or cos.');
end
if strcmp(p.kind, 'cos')
    if ~(isfield(p, 'amplitude') && isscalar(p.amplitude) && p.amplitude > 0 && p.amplitude < 1)
        error('SA:TRACE:E_ENTRY', 'perturbation.amplitude must lie in (0, 1).');
    end
    if ~isfield(p, 'mode') || isempty(p.mode), p.mode = 1; end
    if ~(isscalar(p.mode) && p.mode >= 1 && p.mode == fix(p.mode))
        error('SA:TRACE:E_ENTRY', 'perturbation.mode must be a positive integer.');
    end
    e.perturbation = struct('kind', 'cos', 'amplitude', p.amplitude, 'mode', p.mode);
else
    e.perturbation = struct('kind', 'none');
end
if ~(isscalar(e.nx) && e.nx >= 5 && mod(e.nx, 2) == 1)
    error('SA:TRACE:E_ENTRY', 'nx must be an odd integer >= 5 (a node on the seam).');
end
if ~(isnumeric(e.t_out) && isvector(e.t_out) && numel(e.t_out) >= 2 && all(diff(e.t_out) > 0))
    error('SA:TRACE:E_ENTRY', 't_out must be an increasing vector of at least two times.');
end
e.t_out = reshape(double(e.t_out), 1, []);
if ~(isscalar(e.profile_every) && e.profile_every >= 1 && e.profile_every == fix(e.profile_every))
    error('SA:TRACE:E_ENTRY', 'profile_every must be a positive integer.');
end
end

function [k, A] = nearest_root(run_dir, R, entry, mvreq)
% The root closest to entry.at on the figure: plotted coordinates of every
% root (seeded_atlas_plot's own route), each axis divided by its plotted
% range, Euclidean norm; the first of ties.
Q = seeded_atlas_plot(run_dir, struct('coordinates_only', true, 'main_variable', mvreq));
if ~isequal(Q.root_id, {R.root_id})
    error('SA:TRACE:E_ROOT', 'trace %s: the plotted root set differs from the loaded one.', entry.name);
end
ok = all(isfinite(Q.coord), 2);
if ~any(ok)
    error('SA:TRACE:E_ROOT', 'trace %s: no root of %s has finite plotted coordinates.', entry.name, run_dir);
end
raw = Q.coord - entry.at;
dist = sqrt(sum((raw ./ Q.range) .^ 2, 2));
dist(~ok) = Inf;
[ds, o] = sort(dist);
k = o(1);
A = struct('at', entry.at, 'root_id', Q.root_id{k}, 'family', Q.family{k}, 'coord', Q.coord(k, :), ...
    'scaled_distance', ds(1), 'raw_difference', raw(k, :), 'axis_range', Q.range, ...
    'coord_names', {Q.coord_names}, 'next_root_id', '', 'next_scaled_distance', NaN, 'n_roots', numel(R));
if numel(o) >= 2 && isfinite(ds(2))
    A.next_root_id = Q.root_id{o(2)}; A.next_scaled_distance = ds(2);
end
fprintf(['trace %s: at %s -> nearest root %s (%s) at %s; distance %.3g scaled (axis ranges %s), raw per ' ...
    'axis %s (%s)%s\n'], entry.name, mat2str(entry.at, 9), A.root_id, A.family, mat2str(A.coord, 9), ...
    A.scaled_distance, mat2str(A.axis_range, 4), mat2str(A.raw_difference, 3), strjoin(Q.coord_names, ', '), ...
    next_text(A));
end

function t = next_text(A)
t = '';
if ~isempty(A.next_root_id)
    t = sprintf('; next nearest %s at %.3g', A.next_root_id, A.next_scaled_distance);
end
end

function why = at_reason(tfile, key, why)
% A clearer reason when an at entry now resolves to another root
try
    Z = load(tfile, 'key');
    if isfield(Z.key, 'entry') && isfield(Z.key.entry, 'at') && isequal(Z.key.entry, key.entry) && ...
            ~strcmp(Z.key.root_id, key.root_id)
        why = sprintf('the nearest root to at %s changed from %s to %s (the root set changed); %s', ...
            mat2str(key.entry.at, 9), Z.key.root_id, key.root_id, why);
    end
catch
end
end

function [need, why] = decide(tfile, key)
if exist(tfile, 'file') ~= 2
    need = true; why = 'no trace file'; return
end
try
    Z = load(tfile, 'key', 'stop', 't');
catch
    need = true; why = 'trace file unreadable'; return
end
if ~isfield(Z, 'key') || ~isfield(Z, 't') || ~isfield(Z, 'stop')
    need = true; why = 'trace file incomplete'; return
end
if isequal(Z.key, key)
    need = false; why = 'key unchanged'; return
end
need = true;
parts = fieldnames(key);
d = parts(cellfun(@(f) ~isfield(Z.key, f) || ~isequal(Z.key.(f), key.(f)), parts));
why = sprintf('key changed: %s', strjoin(d(:).', ', '));
end

function t = summary_text(T, P, run_dir)
L = {};
L{end+1} = sprintf('trace %s  (%s)', T.entry.name, T.provenance.software_version);
L{end+1} = sprintf('run folder %s', run_dir);
L{end+1} = sprintf('root %s (representative attempt %s)', T.root_id, T.key.attempt_id);
if ~isempty(fieldnames(T.at_resolution))
    A = T.at_resolution;
    L{end+1} = sprintf(['started from at %s: nearest root on the figure, coordinates %s (%s); distance %.3g ' ...
        'scaled by the axis ranges %s, raw per axis %s%s'], mat2str(A.at, 9), mat2str(A.coord, 9), ...
        strjoin(A.coord_names, ', '), A.scaled_distance, mat2str(A.axis_range, 4), mat2str(A.raw_difference, 3), ...
        next_text(A));
end
p = T.entry.perturbation;
if strcmp(p.kind, 'cos')
    L{end+1} = sprintf('perturbation cos: amplitude %g, mode %d', p.amplitude, p.mode);
else
    L{end+1} = 'perturbation none';
end
L{end+1} = sprintf('mesh nx %d; t_out %g .. %g (%d times); tolerances RelTol %g AbsTol %g', T.entry.nx, ...
    T.entry.t_out(1), T.entry.t_out(end), numel(T.entry.t_out), T.entry.reltol, T.entry.abstol);
L{end+1} = sprintf('stop %s%s', T.stop, tern(isempty(T.stop_detail), '', [' (' T.stop_detail ')']));
L{end+1} = sprintf('frames %d, t_end %.6g; profiles stored at %d frames; wall %.1f s', numel(T.t), T.t(end), ...
    numel(T.profile_frames), T.wall_s);
L{end+1} = sprintf('coordinates (%s): frame 0 %s, last %s', strjoin(P.coord_names, ', '), ...
    mat2str(T.coord(1, :), 6), mat2str(T.coord(end, :), 6));
dr = max(T.mass_drift_rel, [], 1);
L{end+1} = sprintf('pool masses (%s): frame 0 %s; largest relative drift %s', strjoin(P.mass_names, ', '), ...
    mat2str(T.masses(1, :), 10), mat2str(dr, 3));
if any(dr > P.drift_warn)
    L{end+1} = sprintf('WARNING: the relative mass drift exceeds %g (%s): %s', P.drift_warn, mat2str(dr, 3), ...
        P.drift_warn_text);
end
L{end+1} = sprintf('frame-0 gap to the plotted coordinates %s (%s)', mat2str(T.frame0_gap, 3), T.frame0_gap_note);
L{end+1} = 'A trace end state is not shown to be an attractor; no stability is claimed.';
t = [strjoin(L, newline) newline];
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function addpath_once(p)
if ~contains([pathsep path pathsep], [pathsep p pathsep]), addpath(p); end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
