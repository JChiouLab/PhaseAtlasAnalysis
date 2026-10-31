function varargout = sa_layout(action, varargin)
%SA_LAYOUT  Run-folder layout: human-facing files on top, state in internal/.
%
%   D = SA_LAYOUT('state', run_dir)
%       The state folder of a run folder: <run_dir>/internal when it exists
%       or when the folder is new (missing, empty or holding only the human
%       files); <run_dir> itself for a flat folder written before the
%       internal/ layout (state files on top).
%   t = SA_LAYOUT('layout', run_dir)   'v2' | 'flat' | 'none'
%   SA_LAYOUT('create', run_dir)
%       Makes internal/attempts, internal/roots, plots/ and traces/ in an
%       existing empty run folder and writes README.txt.
%   t = SA_LAYOUT('readme')            the README.txt text
%   f = SA_LAYOUT('human_files')       the human-facing file names on top
%
%   Layout:
%     <run_dir>/  README.txt summary.txt roots_index.csv attempts.csv
%                 settings_effective.txt [pilot_report.txt]
%       plots/    roots_3d.fig, roots_3d.csv
%       traces/   <name>/trace.mat, roots_3d_trace.fig, trace_summary.txt
%       compare/  (only after seeded_atlas_compare)
%       internal/ settings_effective.mat settings_user.mat journal.txt
%                 attempts/ roots/ attempts.mat queue.mat queue_final.mat
%                 process_time.mat resume_<n>.mat identity_pairs.csv
%                 pilot/*.csv cutset_* compaction.txt/.mat
%   internal/ holds exactly what a flat folder held minus the human files,
%   so every state reader and writer works on the state folder unchanged.

switch action
    case 'state',       varargout{1} = state_dir(varargin{1});
    case 'layout',      varargout{1} = layout_of(varargin{1});
    case 'create',      create(varargin{1});
    case 'readme',      varargout{1} = readme_text();
    case 'human_files', varargout{1} = human_files();
    otherwise, error('SA:LAYOUT', 'Unknown sa_layout action %s.', action);
end
end

% =========================================================================

function f = human_files()
f = {'README.txt', 'summary.txt', 'roots_index.csv', 'attempts.csv', 'settings_effective.txt', ...
    'pilot_report.txt'};
end

function t = layout_of(d)
if exist(fullfile(d, 'internal'), 'dir') == 7
    t = 'v2';
elseif exist(fullfile(d, 'settings_effective.mat'), 'file') == 2 || ...
        exist(fullfile(d, 'settings_user.mat'), 'file') == 2 || ...
        exist(fullfile(d, 'journal.txt'), 'file') == 2 || exist(fullfile(d, 'attempts'), 'dir') == 7
    t = 'flat';
else
    t = 'none';
end
end

function D = state_dir(d)
if strcmp(layout_of(d), 'flat')
    D = d;
else
    D = fullfile(d, 'internal');
end
end

function create(d)
for s = {fullfile('internal', 'attempts'), fullfile('internal', 'roots'), 'plots', 'traces'}
    [ok, msg] = mkdir(fullfile(d, s{1}));
    if ~ok, error('SA:LAYOUT', 'Cannot create %s: %s', fullfile(d, s{1}), msg); end
end
sa_write('text', fullfile(d, 'README.txt'), readme_text());
end

function t = readme_text()
L = {
sprintf('Seeded coupled atlas run folder (%s)', sa_settings('version'))
''
'Read first'
'  summary.txt             what the run did: stop reason, calls, roots, and why calls failed'
'  roots_index.csv         one row per root: root_id, lattice key, coordinate, routes, observed family'
'  attempts.csv            one row per solver call; the last column (reason) says in plain words'
'                          why a call failed'
'  settings_effective.txt  every setting: its value and where it came from (package, model, mode,'
'                          user, derived)'
'  pilot_report.txt        pilot mode only: cost, yield and the projection for a production run'
'  plots/roots_3d.fig      the 3-D root plot (openfig or double-click; rotate with the figure tool)'
'  plots/roots_3d.csv      the plotted coordinates (roots_3d_key.mat: when to redraw)'
'  traces/<name>/          full-PDE traces started from roots: trace.mat, roots_3d_trace.fig,'
'                          trace_summary.txt (roots_3d_trace_key.mat: when to redraw)'
'  compare/                only after seeded_atlas_compare'
''
'Evidence'
'  Load roots for analysis with seeded_atlas_roots(run_dir): it re-derives every verdict from'
'  the stored solver calls. The plot and plots/roots_3d.csv use an unverified shard load (decoded,'
'  not re-derived): they are a display, not evidence.'
''
'internal/'
'  The pipeline''s own state: settings_effective.mat, settings_user.mat, journal.txt,'
'  attempts/A######.mat (one file per solver call), roots/ (one root store per lattice node),'
'  attempts.mat, queue.mat, queue_final.mat, process_time.mat, resume_<n>.mat,'
'  identity_pairs.csv, pilot/*.csv. Do not edit, move or delete these files.'
'  The file column of attempts.csv is relative to internal/.'
'  Every call''s file is kept, so the run can be resumed and every root re-verified.'
'  After a completed run, seeded_atlas_compact(run_dir) can remove the files of calls that'
'  found no root (listed with their SHA-256 in internal/compaction.txt); those attempts.csv'
'  rows then point to removed files, and the run can no longer be resumed.'
};
t = [strjoin(L, newline) newline];
end
