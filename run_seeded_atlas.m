% run_seeded_atlas.m
% Phase atlas of a polarity model: from the half-cell steady states (the "uncoupled
% library") to the coupled quasi-steady states ("roots"), plots and full-model
% traces. Edit the block below and run this file.
%
% Terms: a root is one quasi-steady state of the coupled two-half problem at a grid
% point of the plane. A family label names the sheet it lies on, as the shape of the
% left half / right half: 0 = no-peak, 1+ (left) or 1- (right) = half-peak,
% 2- = full-peak (so 1+/0 = half-peak/no-peak, 2-/1- = full-peak/half-peak).
%
% What it does, in order (paths are printed before any work):
%   1. the uncoupled library: reused when its three files exist and hold
%      these parameters, otherwise built (Stages 1-3 of
%      run_atlas_surface_pipeline.m only; a build opens a parallel pool)
%   2. the coupled run: new, resumed, or reused when it is complete; a run
%      folder whose settings or parameters differ is refused (choose a new
%      run_dir), never overwritten
%   3. the 3-D root plot (the 3-D atlas, as in Fig. 3B), plots/roots_3d.fig, and the
%      2-D sheets, plots/sheets_2d/<family>.fig and all_families.fig (zero solves; the
%      root files are read one at a time, so a production run plots in
%      about the memory of one root file plus the plotted roots)
%   4. the traces (full-model simulations) listed below; rerunning reuses every
%      unchanged trace and integrates only new or changed entries
% Rerunning on a completed run solves and replays nothing (steps 3-4 only).
% If the package code changed since the run started, the run folder is
% refused instead; plot and trace it directly (docs/USER_GUIDE.md section 7):
%   seeded_atlas_plot(run_dir)  and  seeded_atlas_trace(run_dir, entry)
% On the screen: a [progress] line when each stage starts and every 5 minutes.
% bvp4c's mesh-limit warnings are recorded for every call but not printed.
% Guide: docs/USER_GUIDE.md

%% ---- edit here --------------------------------------------------------------
model    = 'original_polarity';            % 'original_polarity' (S. pombe, Eq. 5) | 'guzzo_noxa' (M. xanthus, Eqs. 1-4)
par      = struct('a', 2.38/4, 'b', 0.185, 'c', 151.6/2, 'd', 0.552, 'e', 0.828, 'f', 1, ...
                  'Dm', 0.05, 'Dc', 5, 'L', 7, 'M_BG', 2*7, 'M_P', 2*0.246*7);  % S. pombe, Table 2 (unrounded; see below)
% The values above are the S. pombe model of Fig. 3 (L = 7). Fields a b c d e f = kgg, dg,
% kpg, dgp, k*, kg of Table 2 (Table 2 rounds them); M_BG = MG = 2L, M_P = MP = 2*rho*L with
% rho = 0.246 (Table 2 rounds rho to 0.25). Ashbya (Fig. 5/S5): L = 30 and rho = 0.25/3 (Table 2).
% M. xanthus (Table 1): model = 'guzzo_noxa' and, for example,
% par = struct('L', 6, 'r_speed', 0.35, 'dba_factor', 1, 'rho_R', 1, 'rho_B', 1, ...
%                'dn_half_range', 1.8); res = 25;
% r_speed = 0.35 scales kR, kRB, dR to Table 1's kr, krb, dr; dba_factor is the signal factor f
% that multiplies dba (Fig. 2E varies it); dn_half_range is the half range of the dn_R, dn_B
% axes (density differences dn = Delta n/(L/2), so the paper's fraction Delta n/M = dn/2 for rho = 1).
res      = 50;                             % uncoupled library grid (M. xanthus: also passed in par)
max_mode = 2;                              % higher-mode half-cell libraries built beside the base files: 2 = none;
                                           %   e.g. 12 builds branch_MODE3..12_latest.mat (stops after the first
                                           %   empty mode; existing files reused), which seed_modes 'all' uses
points_per_axis = 51;                      % coupled grid points per axis (odd), each over its own range of Delta n
                                           %   [] = a fixed step instead: 0.02 (M. xanthus 0.15), or
                                           %   settings_overrides.lattice.step = 0.01 below
library_dir = '';                          % '' = data/<uncoupled run name>/outputs; built if missing
% library_confirmation = '';               % M. xanthus only: a confirmation file with seedable masks
run_dir  = '';                             % '' = data/standalone_runs/<model>_<mode>_<parameters>_<grid>_<hash>
                                           %   e.g. original_polarity_pilot_a0p595_..._L7_MBG14_MP3p444_n51_92260295
                                           %   (grid: n51 = points_per_axis, step0p02 = a fixed step)
run_mode = 'pilot';                        % 'smoke' | 'pilot' | 'production' | 'focus'
                                           %   (smoke: only the 3 x 3 grid points around the centre; pilot: those plus a
                                           %   sample over the library's range, stepped along two lines through the centre;
                                           %   production: the whole grid, i.e. the full atlas: raise budget below)
                                           %   focus: chosen roots of an earlier run (import_from, required) are the
                                           %   only starting profiles, the search fenced to their families; budget required;
                                           %   no library is built or read (res, max_mode, library_dir and
                                           %   seed_modes are ignored: the field scale and the family lookup
                                           %   come from import_from) (docs/USER_GUIDE.md section 3, "Focus mode")
from_pilot = '';                           % production: a pilot run folder whose settings to keep
import_from = '';                          % '' = none; a finished run folder whose roots seed this run
                                           %   (each re-solved once at the nearest node of this grid; same
                                           %   model and parameters). Import only (no cold calls: no solves started
                                           %   from the library):  settings_overrides.routes.cold = false below
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

% ---- useful overrides (docs/USER_GUIDE.md section 6 explains each) -----------
% settings_overrides.model_options.main_variable = 'u+q';   % the variable the family label (a shape per
%                                                            %   half, e.g. 1+/0) and the plots' third axis dn2 look at: a
%                                                            %   '+'-sum of the model's fields (S. pombe
%                                                            %   u q v w p, default 'u', as in the paper; M. xanthus Rm Rc Bm Bc,
%                                                            %   default 'Bm'); a new run. To redraw a finished
%                                                            %   run with another one (its records unchanged):
%                                                            %   seeded_atlas_plot(run_dir, struct('main_variable', 'u+q'))
% settings_overrides.model_options.shape_flat_tol = 0.5606;  % a half is flat (0) when max - min of the main
%                                                            %   variable, divided by its value in the uniform
%                                                            %   state, is below this (S. pombe 0.5606,
%                                                            %   M. xanthus 0.2)
% settings_overrides.model_options.shape_turn_frac = 0.2;    % a turn of the profile counts when it moves back
%                                                            %   by this fraction of the half's max - min
% (the grid is points_per_axis above; one count per axis, e.g. [25 99], is also allowed there)
% settings_overrides.lattice.step = 0.01;                    % a fixed step: only with points_per_axis = []
% settings_overrides.lattice.max_abs_coordinate = [0.1 0.3]; % keep only |coordinate| <= cap on each axis
% settings_overrides.solver.nmax = 40000;                    % more mesh points per call (also raise mesh_reject)
% settings_overrides.solver.mesh_reject = 40000;
% settings_overrides.continuation.guess_nodes = 802;         % continuation guess mesh
% settings_overrides.cold.mesh_cap = Inf;                    % no mesh cap on solves from the library (slower)
% settings_overrides.root_pass.enabled = false;              % no per-root pass (stepping to neighbouring points)
% settings_overrides.rescue.rungs = {'diagonal', 'secant', 'deep'};  % the full rescue ladder (retries of failed points)
% settings_overrides.rescue.ring = true;                     % rescue also outside the library's range
% settings_overrides.targets.policy = 'all_seedable';        % a solve from the library at every seedable (family, node) pair
%                                                            %   (production default 'per_seed': one per library profile)
% settings_overrides.targets.thin = struct('families', {{'HSS_HSS'}}, 'stride', 2);  % fewer cold HSS_HSS targets
%                                                            %   (all_seedable only)
% settings_overrides.budget.wall_s = 3600;                   % stop this process cleanly after 3600 s (resume later)
% settings_overrides.solver.reltol = 1e-7;                   % solver tolerances (a new run)
% settings_overrides.solver.abstol = 1e-9;
% settings_overrides.report.progress_interval_s = 60;        % progress status lines every 60 s (default 300)
% settings_overrides.report.files_interval_s = 1800;         % rewrite attempts.csv, summary.txt, roots_index.csv
%                                                            %   every 30 min during the run (default 600 s)
% settings_overrides.report.solver_warnings = true;          % print bvp4c's warnings again (debugging)
% settings_overrides.routes.cold = false;                    % with import_from: import only (no cold calls)
% settings_overrides.import.mirror = false;                  % no mirror images of the imported roots
% settings_overrides.import.source_edges = 'same_or_stronger';  % same grid, solver settings and model
%                                                            %   only: skip edges the source refused ('all': also
%                                                            %   its weaker per-root pass refusals; default 'none';
%                                                            %   focus default 'auto')
% settings_overrides.import.root_ids = {'k:1,4/R000001'};    % import only these source roots (or the path of a
%                                                            %   CSV with a root_id column, e.g. rows of the
%                                                            %   source's roots_index.csv)
% settings_overrides.import.families = {'2-/1-'};            % and/or the source roots with these labels
% settings_overrides.continuation.families = {'2-/1-', '1+/2-'};  % the family fence: roots with other labels are
%                                                            %   kept but never stepped from ({} = none; focus
%                                                            %   default 'selected': the chosen labels + mirrors)
% settings_overrides.import.source_stopped = true;           % import_from was killed and will not write again

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
    'mode', run_mode, 'points_per_axis', points_per_axis, 'from_pilot', from_pilot, 'import_from', import_from, ...
    'max_mode', max_mode, ...
    'budget', budget, ...
    'settings_overrides', settings_overrides, ...
    'traces', traces, 'make_plot', make_plot, 'compact_after', compact_after);
if exist('library_confirmation', 'var'), cfg.library_confirmation = library_confirmation; end
W = seeded_atlas_wrapper(cfg);
