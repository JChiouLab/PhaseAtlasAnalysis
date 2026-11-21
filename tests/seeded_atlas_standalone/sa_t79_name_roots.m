function rows = sa_t79_name_roots(E)
%SA_T79_NAME_ROOTS  The default run-folder name shows
%   the grid the run uses (n<count> or step<step>, after the parameters, before
%   the hash; the hash over the parameters only, unchanged), the wrapper takes
%   it from the resolved settings (points_per_axis, or a grid set in
%   settings_overrides); and the 'final reports' / 'pilot report' progress
%   lines give the roots held then. Zero solves: sa_run_name is pure text;
%   the wrapper runs status_only on the stub model (nothing written); the
%   runs use the stub corrector on the stub model.
T = 'T79 run name and roots count';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d0 = sa_m3_testdir(E, 't79');

% ---- the name rule (text only) ---------------------------------------------------------------------------
par = struct('a', 2.38/4, 'b', 0.185, 'c', 151.6/2, 'd', 0.552, 'e', 0.828, 'f', 1, ...
    'Dm', 0.05, 'Dc', 5, 'L', 7, 'M_BG', 2*7, 'M_P', 2*0.246*7);
pre = 'original_polarity_pilot_a0p595_b0p185_c75p8_d0p552_e0p828_f1_Dm0p05_Dc5_L7_MBG14_MP3p444';
owner = [pre '_92260295'];                                   % the folder name before the grid was shown
s = struct('model', 'original_polarity', 'mode', 'pilot');
lat = @(n, h) struct('points_per_axis', n, 'step', h);
n51 = sa_run_name(par, setf(s, 'lattice', lat(51, 0.02)));
n2 = sa_run_name(par, setf(s, 'lattice', lat([25 99], 0.02)));
h1 = sa_run_name(par, setf(s, 'lattice', lat([], 0.02)));
h2 = sa_run_name(par, setf(s, 'lattice', lat([], [0.01 0.05])));
none = sa_run_name(par, s);
rows = sa_check(rows, T, ['the grid follows the parameters, before the unchanged hash: n51, n25x99 (per-axis counts), ' ...
    'step0p02 (fixed step), step0p01x0p05 (per-axis steps); the count wins over the step when both are present; ' ...
    'the shipped L = 7 name differs from the pre-grid folder name'], ...
    strcmp(n51, [pre '_n51_92260295']) && strcmp(n2, [pre '_n25x99_92260295']) && ...
    strcmp(h1, [pre '_step0p02_92260295']) && strcmp(h2, [pre '_step0p01x0p05_92260295']) && ...
    strcmp(none, owner) && ~strcmp(n51, owner) && ~strcmp(h1, owner), sprintf('%s | %s', n51, h1));
pn = struct('L', 6, 'r_speed', 0.35, 'dba_factor', 1, 'rho_R', 1, 'rho_B', 1, 'dn_half_range', 1.8, 'res', 25);
sn = struct('model', 'guzzo_noxa', 'mode', 'pilot');
g51 = sa_run_name(pn, setf(sn, 'lattice', lat(51, 0.15)));
g0 = sa_run_name(pn, setf(sn, 'lattice', lat([], 0.15)));
gp = 'guzzo_noxa_pilot_L6_rspeed0p35_dbafactor1_rhoR1_rhoB1_dnhalfrange1p8_res25';
rows = sa_check(rows, T, 'Guzzo_noXA: the same rule (n51; the fixed step step0p15), hash unchanged', ...
    strcmp(g51, [gp '_n51_6d40c23e']) && strcmp(g0, [gp '_step0p15_6d40c23e']), sprintf('%s | %s', g51, g0));
src = {};
for f = dir(fullfile(E.pkg, '*.m')).'
    if strcmp(f.name, 'sa_run_name.m'), continue; end
    if contains(fileread(fullfile(f.folder, f.name)), 'sa_run_name('), src{end+1} = f.name; end %#ok<AGROW>
end
rows = sa_check(rows, T, ['the name enters no solve, row or key: the only caller is the wrapper, where it gives the ' ...
    'default run_dir (a resume-free leaf)'], isequal(src, {'seeded_atlas_wrapper.m'}) && ...
    any(strcmp(sa_settings('resume_free'), 'run_dir')), strjoin(src, ', '));

% ---- the wrapper: the name from the grid the run uses (status_only, nothing written) ---------------------------
s0 = sa_m3_stub('settings', 'x', [1 1]);
ov = rmfield(s0, {'model', 'run_dir', 'budget'});
ovp = ov; ovp.lattice = rmfield(ovp.lattice, 'step');
base = struct('model', 'stub_two_pool', 'par', sa_m3_stub('params'), 'mode', 'smoke', ...
    'library_manifest', struct(), 'settings_overrides', ovp, 'budget', struct('max_calls', 5000), 'status_only', true);
c = base; c.points_per_axis = 5;
evalc('Wa = seeded_atlas_wrapper(c);');
c = base; c.settings_overrides.lattice.points_per_axis = 7;
evalc('Wb = seeded_atlas_wrapper(c);');
c = base; c.settings_overrides = ov;                         % [] and the stub's step 0.65
evalc('Wc = seeded_atlas_wrapper(c);');
c = base; c.settings_overrides = ov; c.settings_overrides.lattice.step = [0.65 0.5];
evalc('Wd = seeded_atlas_wrapper(c);');
nm = cellfun(@(W) last_part(W.run_dir), {Wa, Wb, Wc, Wd}, 'UniformOutput', false);
hs = cellfun(@(x) x(end-7:end), nm, 'UniformOutput', false);
parent = fullfile(E.root, 'data', 'standalone_runs');
same = isequal(cellfun(@(W) fileparts(W.run_dir), {Wa, Wb, Wc, Wd}, 'UniformOutput', false), repmat({parent}, 1, 4));
rows = sa_check(rows, T, ['wrapper default run_dir (status_only): points_per_axis 5 -> _n5_; a count set in ' ...
    'settings_overrides.lattice.points_per_axis 7 -> _n7_; points_per_axis [] -> _step0p65_ (the step the run ' ...
    'resolves), per-axis steps -> _step0p65x0p5_; one hash; each name''s grid is the status''s resolved grid; ' ...
    'status new and no folder created'], same && endsWith(nm{1}, ['_n5_' hs{1}]) && ...
    endsWith(nm{2}, ['_n7_' hs{1}]) && endsWith(nm{3}, ['_step0p65_' hs{1}]) && ...
    endsWith(nm{4}, ['_step0p65x0p5_' hs{1}]) && numel(unique(hs)) == 1 && ...
    isequal(Wa.status.settings.lattice.points_per_axis, 5) && isequal(Wb.status.settings.lattice.points_per_axis, 7) && ...
    isempty(Wc.status.settings.lattice.points_per_axis) && ...
    all(cellfun(@(W) strcmp(W.status.state, 'new') && exist(W.run_dir, 'dir') ~= 7, {Wa, Wb, Wc, Wd})), ...
    strjoin(nm, ' | '));
c = base; c.points_per_axis = 5; c.run_dir = fullfile(d0, 'explicit');
evalc('We = seeded_atlas_wrapper(c);');
rows = sa_check(rows, T, 'an explicit run_dir is used as given (no grid added)', strcmp(We.run_dir, c.run_dir), We.run_dir);

% ---- the roots count on the last stage lines --------------------------------------------------------------------
% Stage lines only (interval Inf) and a call budget that stops the run inside a stage, so roots are found after the
% last stage line: the count that line held (what the final lines used to repeat) is below the result's n_roots.
sC = sa_m3_stub('settings', fullfile(d0, 'C'), [1 1]); sC.report = struct('progress_interval_s', Inf);
sC.budget.max_calls = 6;
stC = sa_m3_stub('state');
oC = evalc('rC = seeded_atlas_run(sa_m3_stub(''params''), struct(), sC, sa_m3_stub(''corrector'', stC));');
[fr, prev, was] = final_lines(oC, 'final reports');
rows = sa_check(rows, T, ['the final reports line gives the run''s calls and roots (the result''s n_roots); the stage ' ...
    'line before it held a smaller count (the one the line used to repeat)'], rC.n_roots > 0 && ...
    contains(fr, sprintf('calls %d of', rC.calls)) && contains(fr, sprintf('roots %d |', rC.n_roots)) && ...
    ~isnan(was) && was < rC.n_roots, sprintf('stop %s, %d roots; final: %s; before: %s', rC.stop_reason, ...
    rC.n_roots, fr, prev));
sP = sa_m3_stub('settings', fullfile(d0, 'P'), []); sP.report = struct('progress_interval_s', Inf);
sP.mode = 'pilot'; sP.targets = struct('policy', 'pilot_sample', 'list', zeros(0, 2));   % as in sa_t47_pilot_report
sP.continuation.region = 'cross'; sP.budget.max_calls = 30;
stP = sa_m3_stub('state');
oP = evalc('rP = seeded_atlas_run(sa_m3_stub(''params''), struct(), sP, sa_m3_stub(''corrector'', stP));');
[fp, prevP, wasP] = final_lines(oP, 'final reports');
pr = final_lines(oP, 'pilot report');
rows = sa_check(rows, T, ['pilot mode: the final reports and pilot report lines give the result''s calls and roots; ' ...
    'the stage line before them held a smaller count'], rP.n_roots > 0 && ...
    contains(pr, sprintf('calls %d of', rP.calls)) && contains(pr, sprintf('roots %d |', rP.n_roots)) && ...
    contains(fp, sprintf('roots %d |', rP.n_roots)) && ~isnan(wasP) && wasP < rP.n_roots, ...
    sprintf('stop %s, %d roots; %s; before: %s', rP.stop_reason, rP.n_roots, pr, prevP));
end

% =========================================================================

function s = setf(s, f, v)
s.(f) = v;
end

function [line, prev, was] = final_lines(out, name)
% the stage line called name, the stage line before it, and that line's roots count (NaN: none)
st = regexp(out, '\[progress\] stage:[^\n]*', 'match');
k = find(contains(st, ['stage: ' name]), 1);
line = ''; prev = ''; was = NaN;
if isempty(k), return; end
line = st{k};
if k > 1
    prev = st{k - 1};
    t = regexp(prev, 'roots (\d+) \|', 'tokens', 'once');
    if ~isempty(t), was = str2double(t{1}); end
end
end

function n = last_part(p)
[~, n, e] = fileparts(p);
n = [n e];
end
