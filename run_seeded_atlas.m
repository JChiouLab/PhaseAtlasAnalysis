% run_seeded_atlas.m
% Seeded coupled atlas, from the uncoupled seed library to coupled roots,
% plots and full-PDE traces. Edit the block below and run this file.
%
% What it does, in order (paths are printed before any work):
%   1. the uncoupled library: reused when its three files exist and hold
%      these parameters, otherwise built (Stages 1-3 of
%      run_atlas_surface_pipeline.m only; a build opens a parallel pool)
%   2. the coupled run: new, resumed, or reused when it is complete; a run
%      folder whose settings or parameters differ is refused (choose a new
%      run_dir), never overwritten
%   3. the 3-D root plot, plots/roots_3d.fig, and the 2-D sheets,
%      plots/sheets_2d/<family>.fig and all_families.fig (zero solves; the
%      root files are read one at a time, so a production run plots in
%      about the memory of one root file plus the plotted roots)
%   4. the traces listed below; rerunning reuses every unchanged trace and
%      integrates only new or changed entries
% Rerunning on a completed run solves and replays nothing (steps 3-4 only).
% If the package code changed since the run started, the run folder is
% refused instead; plot and trace it directly (USER_GUIDE section 7):
%   seeded_atlas_plot(run_dir)  and  seeded_atlas_trace(run_dir, entry)
% On the screen: a [progress] line when each stage starts and every 5 minutes.
% bvp4c's mesh-limit warnings are recorded for every call but not printed.
% Guide: docs/studies/standalone_seeded_atlas/USER_GUIDE.md

%% ---- edit here --------------------------------------------------------------
model    = 'original_polarity';            % 'original_polarity' | 'guzzo_noxa'
par      = struct('a', 2.38/4, 'b', 0.185, 'c', 151.6/2, 'd', 0.552, 'e', 0.828, 'f', 1, ...
                  'Dm', 0.05, 'Dc', 5, 'L', 7, 'M_BG', 2*7, 'M_P', 2*0.246*7);  % the uncoupled script's expressions
% guzzo_noxa example: par = struct('L', 6, 'r_speed', 0.35, 'dba_factor', 1, 'rho_R', 1, 'rho_B', 1, ...
%                                  'dn_half_range', 1.8); res = 25;
res      = 50;                             % uncoupled seed library grid (noxa: also passed in par)
points_per_axis = 51;                      % coupled grid points per pool axis (odd), each over its own pool's range
                                           %   [] = a fixed step instead: 0.02 (noxa 0.15), or
                                           %   settings_overrides.lattice.step = 0.01 below
library_dir = '';                          % '' = data/<uncoupled run name>/outputs; built if missing
% library_confirmation = '';               % noxa only: a pre-seeding confirmation file
run_dir  = '';                             % '' = data/standalone_runs/<model>_<mode>_<parameters>_<grid>_<hash>
                                           %   e.g. original_polarity_pilot_a0p595_..._L7_MBG14_MP3p444_n51_92260295
                                           %   (grid: n51 = points_per_axis, step0p02 = a fixed step)
run_mode = 'pilot';                        % 'smoke' | 'pilot' | 'production'
from_pilot = '';                           % production: a pilot run folder whose settings to keep
budget   = struct('max_calls', 500, 'max_wall_s', 7200);  % ALWAYS the limit, in every mode (production too):
                                           %   the run stops at max_calls calls in total or after max_wall_s seconds
                                           %   of this process, whichever comes first. Running the file again resumes
                                           %   (after a max_calls stop, raise max_calls first). Production has no
                                           %   budget of its own, so raise these for a full atlas.
settings_overrides = struct();             % further settings; see the list below
traces   = struct('name', {}, 'root', {}); % trace entries; see below
make_plot = true;
compact_after = false;                     % true: after a completed run, remove the files of calls
                                           % that found no root (no resume afterwards)

% ---- useful overrides (USER_GUIDE section 6 explains each) -------------------
% settings_overrides.model_options.main_variable = 'u+q';   % the variable the family label (a mode name per
%                                                            %   half, e.g. 1+/0) and the plots' third axis dn2 look at: a
%                                                            %   '+'-sum of the model's fields (fission yeast
%                                                            %   u q v w p, default 'u'; Guzzo_noXA Rm Rc Bm Bc,
%                                                            %   default 'Bm'); a new run. To redraw a finished
%                                                            %   run with another one (its records unchanged):
%                                                            %   seeded_atlas_plot(run_dir, struct('main_variable', 'u+q'))
% settings_overrides.model_options.shape_flat_tol = 0.5606;  % a half is flat (0) when max - min of the main
%                                                            %   variable, divided by its value in the uniform
%                                                            %   state, is below this (fission yeast 0.5606,
%                                                            %   Guzzo_noXA 0.2)
% settings_overrides.model_options.shape_turn_frac = 0.2;    % a turn of the profile counts when it moves back
%                                                            %   by this fraction of the half's max - min
% (the grid is points_per_axis above; one count per pool, e.g. [25 99], is also allowed there)
% settings_overrides.lattice.step = 0.01;                    % a fixed step: only with points_per_axis = []
% settings_overrides.lattice.max_abs_coordinate = [0.1 0.3]; % keep only |coordinate| <= cap on each axis
% settings_overrides.solver.nmax = 40000;                    % more mesh points per call (also raise mesh_reject)
% settings_overrides.solver.mesh_reject = 40000;
% settings_overrides.continuation.guess_nodes = 802;         % continuation guess mesh (the pilot lever)
% settings_overrides.cold.mesh_cap = Inf;                    % no cold mesh cap (slower cold calls)
% settings_overrides.root_pass.enabled = false;              % no per-root continuation pass
% settings_overrides.rescue.rungs = {'diagonal', 'secant', 'deep'};  % the full rescue ladder
% settings_overrides.rescue.ring = true;                     % rescue also outside the library's range
% settings_overrides.targets.policy = 'all_seedable';        % a cold call at every seedable (family, node) pair
%                                                            %   (production default 'per_seed': one per library seed)
% settings_overrides.targets.thin = struct('families', {{'HSS_HSS'}}, 'stride', 2);  % fewer cold HSS_HSS targets
%                                                            %   (all_seedable only)
% settings_overrides.budget.wall_s = 3600;                   % stop this process cleanly after 3600 s (resume later)
% settings_overrides.solver.reltol = 1e-7;                   % solver tolerances (a new run)
% settings_overrides.solver.abstol = 1e-9;
% settings_overrides.report.progress_interval_s = 60;        % progress status lines every 60 s (default 300)
% settings_overrides.report.files_interval_s = 1800;         % rewrite attempts.csv, summary.txt, roots_index.csv
%                                                            %   every 30 min during the run (default 600 s)
% settings_overrides.report.solver_warnings = true;          % print bvp4c's warnings again (debugging)

% ---- traces (one struct per trace; [] fields take the model's defaults) -------
% Start from a root id, or from a position on plots/roots_3d.fig: click a root in
% the figure, copy the data tip's "at" row, and paste it; the root closest on the
% figure (each axis scaled by its plotted range) is used and printed.
% traces(1).name = 'C1';                                     % folder traces/C1/
% traces(1).root = 'k:1,4/R000001';                          % a root_id from roots_index.csv, or
% traces(1).at = [0.05 -0.12 0.31];                          % a position [x y z] on the figure (not with root)
% traces(1).perturbation = struct('kind', 'cos', 'amplitude', 0.1, 'mode', 1);
% traces(1).t_out = [0:0.05:50, 51:1:500];                   % output times
% traces(1).nx = 1001;                                       % full-domain mesh (odd)

%% ---- run (no edits below) ------------------------------------------------------
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src', 'pipelines', 'seeded_atlas_standalone'));
addpath(here);
setup_atlas('quiet');                     % the traced bvp4c copy: made or remade when missing or stale
cfg = struct('model', model, 'par', par, 'res', res, 'library_dir', library_dir, 'run_dir', run_dir, ...
    'mode', run_mode, 'points_per_axis', points_per_axis, 'from_pilot', from_pilot, 'budget', budget, ...
    'settings_overrides', settings_overrides, ...
    'traces', traces, 'make_plot', make_plot, 'compact_after', compact_after);
if exist('library_confirmation', 'var'), cfg.library_confirmation = library_confirmation; end
W = seeded_atlas_wrapper(cfg);
