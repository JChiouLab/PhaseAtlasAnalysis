function rows = sa_t46_mode_equivalence(E)
%SA_T46_MODE_EQUIVALENCE  Milestone-4 T46 (packet v0.3 section 2a): with every
%   mode-layer leaf set explicitly, the mode changes only the reports, never
%   the calls: smoke, pilot and production stub runs are strictly equivalent
%   (sa_m3_equiv: units, sol, validity, opts, attempt rows, shards) and only
%   the pilot run writes pilot_report.txt. Every run writes roots_index.csv.
T = 'T46 mode equivalence';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't46');
fail = @(dd, b, ex) sa_m3_wall(dd, ex);
o = struct('targets', struct('policy', 'explicit', 'list', [0 0]), ...
    'routes', struct('cold', true, 'continuation', true, 'rescue', true), ...
    'continuation', struct('region', 'box', 'box', [3 1]), 'rescue', struct('ring', true), ...
    'budget', struct('max_calls', 5000, 'max_wall_s', 3600, 'wall_s', Inf));
R = struct();
for m = {'smoke', 'pilot', 'production'}
    oo = o; oo.mode = m{1};
    R.(m{1}) = sa_m3_run(fullfile(d, m{1}), [3 1], 'fail', fail, 'settings', oo);
end
D1 = sa_m3_equiv(fullfile(d, 'smoke'), fullfile(d, 'pilot'));
D2 = sa_m3_equiv(fullfile(d, 'smoke'), fullfile(d, 'production'));
rows = sa_check(rows, T, 'explicit leaves: smoke == pilot == production call for call (units, sol, validity, opts, rows, shards)', ...
    D1.pass && D2.pass && R.smoke.calls == R.pilot.calls && R.smoke.calls == R.production.calls && ...
    R.smoke.calls > 20, sprintf('calls %d/%d/%d; %s', R.smoke.calls, R.pilot.calls, R.production.calls, ...
    strjoin([D1.diffs, D2.diffs], ' | ')));
ex = @(m, f) exist(fullfile(d, m, f), 'file') == 2;
rows = sa_check(rows, T, 'reports differ only: pilot_report.txt in pilot only; summary.txt and roots_index.csv in all', ...
    ex('pilot', 'pilot_report.txt') && ~ex('smoke', 'pilot_report.txt') && ~ex('production', 'pilot_report.txt') && ...
    ex('smoke', 'summary.txt') && ex('production', 'roots_index.csv') && ex('pilot', 'roots_index.csv'));
ri = strsplit(strtrim(fileread(fullfile(d, 'smoke', 'roots_index.csv'))), newline);
rows = sa_check(rows, T, 'roots_index.csv: one row per root', numel(ri) - 1 == R.smoke.n_roots, ...
    sprintf('%d rows, %d roots', numel(ri) - 1, R.smoke.n_roots));
rmpath(E.stub);
end
