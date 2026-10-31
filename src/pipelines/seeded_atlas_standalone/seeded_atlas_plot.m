function out = seeded_atlas_plot(run_dir, opts)
%SEEDED_ATLAS_PLOT  The root plots of a run folder: a rotatable 3-D .fig and 2-D sheets.
%
%   out = SEEDED_ATLAS_PLOT(run_dir)
%   out = SEEDED_ATLAS_PLOT(run_dir, opts)
%
%   Writes <run_dir>/plots/roots_3d.fig and roots_3d.csv (or into opts.out), and
%   the redraw key roots_3d_key.mat beside them (a trace figure: the .fig and
%   its key only); and the 2-D sheets plots/sheets_2d/<family>.fig (one per
%   observed family; the file name is sa_mode_name 'file_tag' of the
%   family: '+' -> p, '-' -> m, '/' -> _, so 1+/0 is 1p_0.fig and 2-/1- is
%   2m_1m.fig) and plots/sheets_2d/all_families.fig (every family
%   overlaid, with a legend), with their key sheets_2d_key.mat.
%   Zero solves, no deval: roots come from seeded_atlas_roots (decoded
%   shards, read one at a time; opts.verify = true re-derives every verdict
%   first), loaded once for the 3-D plot and the sheets. Each root's
%   coordinates are computed from its own stored representative profile
%   (identity rows on the comparison mesh): the two halves are assembled
%   into the full domain (post full_profile, opts.nx points, default the
%   model's plot_nx) and passed to the model's observables (plot order:
%   fission yeast dG, dP, dn2; Guzzo_noXA dn_R, dn_B, dn2; dn2 =
%   int cos(2 pi eta) m d(eta) of the main variable m). The profile-computed
%   lattice coordinates are compared with the node's lattice coordinate
%   (the gap, in the CSV and UserData).
%
%   Main variable. By default the figures show the run's own labels (the
%   observed families as stored) and dn2 of the run's main variable
%   (settings model_options.main_variable). A folder written before
%   software 0.10 records none: its labels were computed on the post
%   module's legacy_main_variable (fission yeast u, Guzzo_noXA Bm) and dn2
%   is computed on the current default (u / Bm); a line is printed when
%   that differs from the dn2 drawn before 0.10 (fission yeast u + q).
%   A folder written before software 0.11 stores HSS / HPSS / MPSS family
%   names (a printed line says so); a relabel shows it with mode names.
%   opts.main_variable relabels for plotting: every root's family is
%   recomputed by the model's shape_label_profile (the current label rule,
%   the mode names of sa_mode_name, with that main variable and the run's
%   other model_options where they still exist, else the current defaults;
%   the flatness reference is the folder's stored anchor-interpolated HSS,
%   or for a folder written before software 0.11 the same value recomputed
%   from the folder's library files, or failing that the stored HSS of the
%   nearest native cell, with a printed note; never a solve) on its stored
%   representative profile, and dn2 uses
%   that variable. The run's records
%   are not touched; the files get the suffix _mv-<variable> (roots_3d_mv-u+q.fig,
%   plots/sheets_2d_mv-u+q/) and the titles say "relabelled for plotting; run
%   labels as stored". The relabel reads the stored profile on the
%   comparison mesh (interpolated linearly to the label's points), not the
%   attempt's solver solution the run labelled with, so it can differ from a
%   run-time label at a half whose decision sits at a threshold.
%
%   One scatter per observed family, filled markers half transparent (face
%   and edge alpha 0.5) so overlapping sheets show, colours fixed per family
%   name; legend 'FAMILY (n)'. 3-D: view(50, 18) (= matplotlib elevation 18,
%   azimuth -40, in the display frame); orthographic, no grid, no box;
%   fission yeast uses the Fig 3b stretch dG : 3.6 dP : 1.4 dn2. 2-D sheets:
%   the mass-imbalance plane (the 3-D plot's x and y, the same stretch). The
%   titles give the run folder, the software version, the root count, the
%   main variable and PARTIAL RUN (<stop reason>) unless the run completed.
%   Saved with savefig only (no raster or PDF export); each figure is built
%   invisible with a CreateFcn that shows it on openfig.
%
%   opts
%     verify   false (default) | true    seeded_atlas_roots verification
%     out      output folder; REQUIRED for a flat (software 0.5) folder, which
%              is evidence and is never written into
%     nx       full-domain mesh of the coordinates (odd)
%     main_variable  '' (default: the run's own) | a '+'-sum of the model's
%              fields: relabel for plotting (above)
%     sheets   true (default without a trace) | false: the 2-D sheets
%     trace    a trace.mat file: overlay the trace (line, start circle, end
%              square); default file name roots_3d_trace; no sheets
%     file     output file name without extension (default roots_3d)
%     force    redraw even when the stored keys are unchanged (default false)
%     coordinates_only  true: compute and return the plotted coordinates
%              only (out.coord, out.root_id, out.family, out.range); no
%              figure, no file, no key (seeded_atlas_trace's "at" entries)
%   out: fig, csv, redrawn (false when the 3-D key was unchanged), n_roots,
%   gap_max (1 x 2), label, main_variable, relabelled, sheets (the sheet
%   files), sheets_redrawn.
%
%   Data tips: clicking a root in any figure shows two rows, "at" (the
%   root's plotted 3-D coordinates as [x y z], ready to paste into a trace
%   entry's at field) and "root" (its root_id).
%
%   The figure keys (plots/<file>_key.mat, sheets_2d/sheets_2d_key.mat): the
%   root set (SHA-256 of every shard), the plot code digest
%   (seeded_atlas_plot.m and the model's post/ files, comment-free tokens;
%   relabelled: also the model module), the main variable, the relabel
%   options, the trace key (when overlaid) and the options; a figure set is
%   redrawn when any part changed or one of its files is missing.

if nargin < 2, opts = struct(); end
verify = getf(opts, 'verify', false);
D = sa_layout('state', run_dir);
lay = sa_layout('layout', run_dir);
if strcmp(lay, 'none')
    error('SA:PLOT:E_NO_RUN', '%s holds no run.', run_dir);
end
od = getf(opts, 'out', '');
if isempty(od) && ~getf(opts, 'coordinates_only', false)   % coordinates_only writes nothing
    if strcmp(lay, 'flat')
        error('SA:PLOT:E_LEGACY', ['%s is a flat (software 0.5 or earlier) run folder: it is evidence and ' ...
            'is never written into; pass opts.out (another folder).'], run_dir);
    end
    od = fullfile(run_dir, 'plots');
end
trace_file = getf(opts, 'trace', '');
here = fileparts(mfilename('fullpath'));
addpath_once(fullfile(here, 'post'));
E = load(fullfile(D, 'settings_effective.mat'), 'settings', 'lattice', 'theta', 'field_scale', ...
    'library_manifest', 'parameters');
V = main_variable(E, getf(opts, 'main_variable', ''));
P = V.P;
fname = getf(opts, 'file', '');
if isempty(fname)
    if isempty(trace_file), fname = 'roots_3d'; else, fname = 'roots_3d_trace'; end
end
sheets_dir = fullfile(od, 'sheets_2d');
if V.relabel
    fname = [fname '_' V.tag];
    sheets_dir = fullfile(od, ['sheets_2d_' V.tag]);
end
nx = getf(opts, 'nx', P.plot_nx);
if getf(opts, 'coordinates_only', false)
    [R, C, ~, ~, fam] = coordinates(run_dir, verify, E, P, nx, V);
    out = struct('coord', C, 'root_id', {{R.root_id}}, 'family', {fam}, 'range', axis_range(C, P), ...
        'coord_names', {P.coord_names}, 'nx', nx, 'main_variable', V.text, 'relabelled', V.relabel);
    return
end
if ~isempty(V.note), fprintf('%s\n', V.note); end

% ---- keys --------------------------------------------------------------------
sh = sa_store_io('shard_files', D);              % skips sa_write temporary files
names = sort({sh.name});
shas = cellfun(@(n) sa_sha256(fullfile(D, 'roots', n)), names, 'UniformOutput', false);
key = struct('roots', {[names; shas]}, ...
    'code', sa_code_digest('files', [{fullfile(here, 'seeded_atlas_plot.m')}, P.files, V.files]), ...
    'verify', logical(verify), 'nx', nx, 'trace', struct(), 'software_version', sa_settings('version'), ...
    'main_variable', V.text, 'relabel', V.relabel, 'relabel_options', V.mo);
if ~isempty(trace_file)
    Zt = load(trace_file, 'key');
    key.trace = Zt.key;
end
figf = fullfile(od, [fname '.fig']);
keyf = fullfile(od, [fname '_key.mat']);
csvf = fullfile(od, [fname '.csv']);
skeyf = fullfile(sheets_dir, 'sheets_2d_key.mat');
want_sheets = getf(opts, 'sheets', isempty(trace_file));
force = getf(opts, 'force', false);
out = struct('fig', figf, 'csv', csvf, 'redrawn', false, 'n_roots', NaN, 'gap_max', [NaN NaN], ...
    'label', '', 'main_variable', V.text, 'relabelled', V.relabel, 'sheets', {{}}, 'sheets_redrawn', false);
need3 = force || ~key_equal(keyf, key, {figf});
needS = false;
if want_sheets
    needS = force || ~key_equal(skeyf, key, {}) || ~sheet_files_exist(skeyf, sheets_dir);
end
if ~need3
    fprintf('plot %s: unchanged (key equal), not redrawn\n', figf);
end
if want_sheets && ~needS
    out.sheets = sheet_paths(skeyf, sheets_dir);
    fprintf('sheets %s: unchanged (key equal), not redrawn\n', sheets_dir);
end
if ~need3 && ~needS, return; end

% ---- roots and coordinates (one load for the 3-D plot and the sheets) -----------
[R, C, lat, gap, fam, fam_stored] = coordinates(run_dir, verify, E, P, nx, V);
A = load(fullfile(D, 'attempts.mat'), 'stop_reason');
n = numel(R);
if verify
    label = 'verified shard load (verdicts re-derived)';
else
    label = 'unverified shard load (decoded, not re-derived)';
end
stop = A.stop_reason;
tt = sprintf('%s | %s | %d roots', run_dir, sa_settings('version'), n);
if ~isempty(V.text), tt = sprintf('%s | main variable %s', tt, V.text); end
if V.relabel
    tt = sprintf('%s: relabelled for plotting; run labels as stored (%d of %d roots changed family)', tt, ...
        sum(~strcmp(fam, fam_stored)), n);
end
if ~strcmp(stop, 'completed'), tt = sprintf('%s | PARTIAL RUN (%s)', tt, stop); end
gap_max = max(gap, [], 1, 'omitnan');
if isempty(gap_max), gap_max = [NaN NaN]; end
ud = struct('label', label, 'run_dir', run_dir, 'software_version', sa_settings('version'), ...
    'n_roots', n, 'stop_reason', stop, 'coord_names', {P.coord_names}, 'nx', nx, ...
    'gap_max', gap_max, 'gap_note', ['|profile-computed lattice coordinates - node coordinate| ' ...
    '(stored comparison-mesh profile, pchip to nx)'], 'key', key, 'trace', struct(), ...
    'main_variable', V.text, 'relabelled', V.relabel, 'main_variable_note', V.note);
if exist(od, 'dir') ~= 7, mkdir(od); end

if need3
    % ---- 3-D figure -----------------------------------------------------------------
    [fig, ax, h, leg] = family_figure(C, R, fam, P, 3);
    trace_info = struct();
    if ~isempty(trace_file)
        Zt = load(trace_file, 'coord', 'stop', 'entry', 'root_id');
        tc = Zt.coord;
        h(end+1) = plot3(ax, tc(:, 1), tc(:, 2), tc(:, 3), '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 1, ...
            'DisplayName', sprintf('trace %s', Zt.entry.name), 'Tag', 'trace');
        plot3(ax, tc(1, 1), tc(1, 2), tc(1, 3), 'o', 'Color', [0.1 0.1 0.1], 'MarkerSize', 7, ...
            'LineWidth', 1.2, 'Tag', 'trace_start', 'HandleVisibility', 'off');
        plot3(ax, tc(end, 1), tc(end, 2), tc(end, 3), 's', 'Color', [0.1 0.1 0.1], 'MarkerSize', 7, ...
            'LineWidth', 1.2, 'Tag', 'trace_end', 'HandleVisibility', 'off');
        leg{end+1} = sprintf('trace %s', Zt.entry.name);
        trace_info = struct('file', trace_file, 'root_id', Zt.root_id, 'stop', Zt.stop, 'n_frames', size(tc, 1));
        Cs = [C; tc];
    else
        Cs = C;
    end
    if ~isempty(h), legend(ax, h, leg, 'Location', 'northeastoutside', 'Interpreter', 'none'); end
    xlabel(ax, P.coord_labels{1}); ylabel(ax, P.coord_labels{2}); zlabel(ax, P.coord_labels{3});
    grid(ax, 'off'); box(ax, 'off');
    if size(Cs, 1) >= 1 && all(isfinite(Cs(:)))
        set(ax, 'DataAspectRatio', P.aspect(Cs));
    end
    camproj(ax, 'orthographic');
    view(ax, 50, 18);
    title(ax, tt, 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 9);
    ud.trace = trace_info;
    fig.UserData = ud;
    savefig(fig, figf);
    close(fig);
    sa_write('mat', keyf, struct('key', key));

    % ---- CSV (the plotted data; the root plot only, not a trace figure) --------------------
    if ~isempty(trace_file), csvf = ''; out.csv = ''; end
    hdr = sprintf('root_id,key,observed_family,%s,%s,%s,lattice_%s,lattice_%s,gap_%s,gap_%s', P.coord_names{:}, ...
        P.coord_names{1}, P.coord_names{2}, P.coord_names{1}, P.coord_names{2});
    if V.relabel, hdr = [hdr ',stored_family']; end
    lines = {sprintf('# %s; %s', label, sa_settings('version')), hdr};
    if ~isempty(V.text)
        lines{1} = sprintf('%s; main variable %s', lines{1}, V.text);
        if V.relabel, lines{1} = [lines{1} ' (relabelled for plotting; run labels as stored: stored_family)']; end
    end
    for k = 1:n
        s = sprintf('"%s","%s",%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.6g,%.6g', R(k).root_id, R(k).key, ...
            fam{k}, C(k, :), lat(k, :), gap(k, :));
        if V.relabel, s = [s ',' fam_stored{k}]; end %#ok<AGROW>
        lines{end+1} = s; %#ok<AGROW>
    end
    if ~isempty(csvf), sa_write('text', csvf, [strjoin(lines, newline) newline]); end
    out.redrawn = true;
    fprintf('plot %s: %d roots (%s), largest lattice-coordinate gap %s\n', figf, n, label, mat2str(gap_max, 3));
end

if needS
    out.sheets = draw_sheets(sheets_dir, skeyf, key, C, R, fam, P, tt, ud);
    out.sheets_redrawn = true;
    fprintf('sheets %s: %d files\n', sheets_dir, numel(out.sheets));
end
out.n_roots = n; out.gap_max = gap_max; out.label = label;
out.coord = C; out.family = fam; out.root_id = {R.root_id};
end

% =========================================================================

function V = main_variable(E, req)
% The main variable the figures use, the post struct built with it and,
% for a relabel, the model module and model_options that recompute labels.
model = E.settings.model;
stored = '';
if isfield(E.settings, 'model_options') && isfield(E.settings.model_options, 'main_variable')
    stored = E.settings.model_options.main_variable;
end
V = struct('relabel', ~isempty(req), 'P', [], 'text', '', 'tag', '', 'note', '', 'M', [], 'mo', struct(), ...
    'files', {{}}, 'fs', []);
if V.relabel
    V.P = sa_post_model(model, req);
    V.M = sa_model(model);
    if ~isfield(V.M, 'shape_label_profile')
        error('SA:PLOT:E_RELABEL', 'Model %s has no shape_label_profile: its roots cannot be relabelled.', model);
    end
    mo = sa_settings('defaults', V.M);
    mo = mo.model_options;
    for f = fieldnames(E.settings.model_options).'        % the run's own model_options where stored
        if isfield(mo, f{1}), mo.(f{1}) = E.settings.model_options.(f{1}); end
    end
    mo.main_variable = V.P.main_variable;
    V.mo = V.M.validate_model_options(mo);
    V.fs = anchor_field_scale(E, V.M, V.mo);
    V.files = {which(['sa_model_' model])};
else
    V.P = sa_post_model(model, stored);
end
if isfield(V.P, 'main_variable')
    V.text = V.P.main_variable;
    V.tag = sa_main_variable('file_tag', V.text);
    if ~V.relabel && isempty(stored)
        V.note = sprintf(['%s records no main variable (written before software 0.10): the family labels are ' ...
            'as stored (computed on %s with the pre-0.10 rule); the third coordinate is dn2 of %s, the current ' ...
            'default'], E.settings.model, V.P.legacy_main_variable, V.text);
        if ~strcmp(V.P.legacy_dn2_variable, V.text)
            V.note = sprintf('%s (figures drawn before 0.10 used %s)', V.note, V.P.legacy_dn2_variable);
        end
        V.note = [V.note '.'];
    end
    if ~V.relabel && ~isfield(E.settings.model_options, 'shape_flat_tol')
        V.note = strtrim(sprintf(['%s The family labels are the names written by software before ' ...
            '0.11; for mode names (1+/0, ...) relabel: opts.main_variable = ''%s''.'], V.note, V.text));
    end
end
end

function fs = anchor_field_scale(E, M, mo)
% The field scale a relabel uses: the folder's own when it carries the
% anchor-interpolated HSS (software 0.11 on); else recomputed from the
% folder's library and settings (files only, no solve); else the stored
% nearest-cell values, with a printed note.
fs = E.field_scale;
if isfield(fs, 'anchor_hss_values'), return; end
try
    lib = M.library('load', E.library_manifest, E.theta, E.parameters);
    S = sa_settings('legacy_fill', E.settings);
    S.model_options = mo;
    f2 = M.field_scale(E.theta, lib, S);
    fs.anchor_hss_values = f2.anchor_hss_values;
    fs.anchor_source = [f2.anchor_source ' (recomputed for the relabel from the folder''s library)'];
    fprintf('relabel: uniform-state reference %s.\n', fs.anchor_source);
catch ME
    fprintf(['relabel: the folder''s library could not be read (%s); the uniform-state reference is the ' ...
        'stored uniform-state values of the native cell nearest the anchor.\n'], ME.message);
end
end

function [R, C, lat, gap, fam, fam_stored] = coordinates(run_dir, verify, E, P, nx, V)
% The plotted coordinates of every root: its stored representative profile
% assembled on the full domain (nx points) and passed to the observables;
% the families as stored, or recomputed from the same profile (relabel)
R = seeded_atlas_roots(run_dir, struct('verify', verify));
n = numel(R);
C = nan(n, 3); gap = nan(n, 2); lat = nan(n, 2);
fam_stored = {R.observed_family};
fam_stored(cellfun(@isempty, fam_stored)) = {'unlabelled'};
fam = fam_stored;
for k = 1:n
    Y = R(k).representative_profile;
    F = P.full_profile(E.theta, R(k).comparison_mesh, Y(P.profile_rows(1, :), :), Y(P.profile_rows(2, :), :), nx);
    O = P.observables(E.theta, F.x, F.Y);
    C(k, :) = O.c;
    lat(k, :) = R(k).coordinate(P.lattice_axes);
    gap(k, :) = abs(O.c(1:2) - lat(k, :));
    if V.relabel
        fam{k} = V.M.shape_label_profile(R(k).comparison_mesh, Y, E.theta, V.mo, V.fs);
    end
end
end

function [fig, ax, h, leg] = family_figure(C, R, fam, P, dim)
% One half-transparent filled scatter per family (dim 3: scatter3 of the
% three coordinates; dim 2: scatter of the first two), in family order
fig = figure('Visible', 'off', 'Color', 'w', 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
ax = axes(fig);
hold(ax, 'on');
order = P.family_order(unique(fam));
h = gobjects(0); leg = {};
for f = 1:numel(order)
    m = strcmp(fam, order{f});
    if ~any(m), continue; end
    nm = sprintf('%s (%d)', order{f}, sum(m));
    if dim == 3
        h(end+1) = scatter3(ax, C(m, 1), C(m, 2), C(m, 3), 12, P.family_colour(order{f}), 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', nm, 'Tag', ['family:' order{f}]); %#ok<AGROW>
    else
        h(end+1) = scatter(ax, C(m, 1), C(m, 2), 12, P.family_colour(order{f}), 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, 'DisplayName', nm, 'Tag', ['family:' order{f}]); %#ok<AGROW>
    end
    h(end).DataTipTemplate.DataTipRows = [dataTipTextRow('at', at_text(C(m, :))), ...
        dataTipTextRow('root', string({R(m).root_id}))];
    leg{end+1} = nm; %#ok<AGROW>
end
end

function files = draw_sheets(sd, skeyf, key, C, R, fam, P, tt, ud)
% The 2-D sheets: one figure per observed family and one with all of them,
% in the mass-imbalance plane (the 3-D plot's x and y, its stretch). Sheet
% files of an earlier drawing that are no longer drawn are removed. The key
% file stores base names only, resolved against sd, so a copied run folder
% never reads or deletes the original folder's sheets.
if exist(sd, 'dir') ~= 7, mkdir(sd); end
order = P.family_order(unique(fam));
names = [order, {'all_families'}];
base = cellfun(@(f) [sa_mode_name('file_tag', f) '.fig'], names, 'UniformOutput', false);
files = cellfun(@(f) fullfile(sd, f), base, 'UniformOutput', false);
if exist(skeyf, 'file') == 2
    try
        Zs = load(skeyf, 'files');
        old = setdiff(sheet_names(Zs.files), base);
        for k = 1:numel(old)
            f = fullfile(sd, old{k});
            if exist(f, 'file') == 2, delete(f); end
        end
    catch
    end
end
ok = all(isfinite(C(:, 1:2)), 2);
a = [1 1 1];
if any(ok), a = P.aspect(C(ok, :)); end
for k = 1:numel(names)
    if strcmp(names{k}, 'all_families'), m = true(size(fam)); what = 'all families';
    else, m = strcmp(fam, names{k}); what = names{k};
    end
    [fig, ax, h, leg] = family_figure(C(m, :), R(m), fam(m), P, 2);
    if ~isempty(h), legend(ax, h, leg, 'Location', 'northeastoutside', 'Interpreter', 'none'); end
    xlabel(ax, P.coord_labels{1}); ylabel(ax, P.coord_labels{2});
    grid(ax, 'off'); box(ax, 'on');
    set(ax, 'DataAspectRatio', [a(1) a(2) 1]);
    title(ax, sprintf('%s | sheet: %s', tt, what), 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 9);
    ud.sheet = what;
    fig.UserData = ud;
    savefig(fig, files{k});
    close(fig);
end
sa_write('mat', skeyf, struct('key', key, 'files', {base}));
end

function same = key_equal(keyf, key, need)
% the stored key equals key and every file of need exists
same = exist(keyf, 'file') == 2 && all(cellfun(@(f) exist(f, 'file') == 2, need));
if ~same, return; end
try
    Zk = load(keyf, 'key');
    same = isequal(Zk.key, key);
catch
    same = false;
end
end

function ok = sheet_files_exist(skeyf, sd)
ok = false;
try
    f = sheet_paths(skeyf, sd);
    ok = ~isempty(f) && all(cellfun(@(p) exist(p, 'file') == 2, f));
catch
end
end

function f = sheet_paths(skeyf, sd)
% the stored sheet names resolved against this folder's sheets_dir
Zs = load(skeyf, 'files');
f = cellfun(@(n) fullfile(sd, n), sheet_names(Zs.files), 'UniformOutput', false);
end

function n = sheet_names(stored)
% bare sheet file names; a key from before this repair stored full paths, of
% which only the name is kept; anything that is not a plain name.fig is dropped
[~, b, e] = cellfun(@fileparts, stored, 'UniformOutput', false);
n = strcat(b, e);
n = n(~cellfun(@isempty, regexp(n, '^[A-Za-z0-9_]+\.fig$', 'once')));
end

function r = axis_range(C, P)
% The plotted range of each axis: max - min over the finite root
% coordinates. An axis whose extent on the figure (range divided by the
% plot's DataAspectRatio component) is below 1 % of the largest axis's
% extent is given that 1 % instead, so a nearly flat axis (few roots, or
% roots on one line) does not dominate the distance. All ranges 1 when no
% axis has any spread.
ok = all(isfinite(C), 2);
r = ones(1, 3);
if ~any(ok), return; end
r = max(C(ok, :), [], 1) - min(C(ok, :), [], 1);
a = P.aspect(C(ok, :));
e = r ./ a;
if ~(max(e) > 0), r = ones(1, 3); return; end
r = max(r, 0.01 * max(e) * a);
end

function s = at_text(C)
% "[x y z]" per root, nine significant digits: pasted into a trace entry's
% at field it picks this root
s = strings(size(C, 1), 1);
for k = 1:size(C, 1), s(k) = sprintf('[%.9g %.9g %.9g]', C(k, :)); end
end

function addpath_once(p)
if ~contains([pathsep path pathsep], [pathsep p pathsep]), addpath(p); end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
