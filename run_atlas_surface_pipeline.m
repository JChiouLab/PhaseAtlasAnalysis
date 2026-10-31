% run_atlas_surface_pipeline.m
% Driver for the HSS/HPSS/MPSS atlas workflow. Runs Stages 1-12 in order:
%   1-3  build the HPSS/HSS/MPSS branch libraries
%   4-7  visualize the static atlas (branch surfaces, per-sheet flux plots,
%        combined 3-D surfaces + nullclines)
%   8-9  run one atlas-seeded PDE trajectory and classify it against all 9 sheets
%   10-12 plot the classified trajectory (phase-plane pages, classification
%        summary) and render it to a movie
%
% Each stage reuses its existing output file unless forced to rebuild (Stages
% 1-3,5) or auto-detects a settings change (Stages 8-9); Stages 4,6,7,10-12
% always run (display only). See docs/atlas_pipeline.md for the full reference.
%
% Related tools NOT run here: infer_atlas_transitions/plot_atlas_transitions
% (transition events), run_atlas_pde_survey (batched trajectories), and
% debug_atlas_point (single-cell inspection) -- call those directly.

%% ---- model parameters ----------------------------------------------------
% Reaction rates a..f, slow/fast diffusion, domain length, two conserved-pool
% masses, and grid resolution. See CLAUDE.md for parameter sets.
% Resolve bundled solvers and results relative to this file rather than the
% MATLAB current folder.  Distribute this file together with
% src/solvers/atlas_surface_pipeline_dependencies/.
pipeline_dir = fileparts(mfilename('fullpath'));
solver_dir = fullfile(pipeline_dir, 'src', 'solvers', 'atlas_surface_pipeline_dependencies');
if ~exist(solver_dir, 'dir')
    error('run_atlas_surface_pipeline:MissingSolverBundle', ...
        'Missing bundled solvers folder: %s', solver_dir);
end
addpath(solver_dir);
model = 'original_polarity';            % 'original_polarity' | 'guzzo_noxa'
if strcmp(model, 'guzzo_noxa')
    addpath(fullfile(pipeline_dir, 'src', 'models', 'guzzo_noxa'));   % shared (R4)
    addpath(fullfile(solver_dir, 'models', 'guzzo_noxa'));
    gnx_run_libraries(gnx_model('parameters', struct('L', 6, 'r_speed', 0.35, ...
        'dba_factor', 1, 'rho_R', 1, 'rho_B', 1, 'res', 25, 'dn_half_range', 1.8)), ...
        pipeline_dir);
    return                               % stages 2b and 4-12 are original-only
end

p.a = 2.38/4; p.b = 0.185; p.c = 151.6/2;
p.d = 0.552;  p.e = 0.828; p.f = 1.0;
p.Dm = 0.05; p.Dc = 5.0;
p.L = 7.94;
p.M_BG = 2*p.L;
p.M_P  = 2*0.246*p.L;    % M_P scan sweet spot: static HPSS present at centre (nP=3, amp~10, well
                 % inside the branch -> best atlas coverage) AND sustained boundary-to-boundary
                 % shuttle. M_P>=8 is past the HPSS fold (nP~4.0) -> no static branch.
p.res = 50;

% Solver options -- only consulted when a sweep/branch stage actually rebuilds.
p.seeding_mode   = 'auto';           % 'auto' (headless) or 'interactive'
p.use_parallel   = true;
p.skip_stability = true;             % true = skip the HPSS eigenvalue stability check
p.mpss_strategy  = 'continuation';   % MPSS construction strategy
p.phase0_construct = true;           % Phase 0 Tier-0: seed HPSS by constructed edge-peak BVP
                                     % (no PDE relaxation) -- reaches the dynamically-unstable
                                     % static half-peak this oscillatory regime lives on.

% Oscillatory-regime settings. This parameter set drives sustained PDE
% oscillation on the half-domain (measured: L=7 oscillates with period ~65;
% the quarter/small-L domains do not). Without these the PDE-relaxation
% fallback never settles and burns its whole budget per cell (hours at res=50).
p.osc_fallback   = true;             % oscillation-aware fallback: bounded relax + snapshot/mean bvp4c seeds
% Short CHUNK length is the key bound: it makes the oscillation-detected
% early-stop (3 chunks) fast (~4s at L=7) while spanning ~30 oscillation
% periods so oscillation is detectable. The generous TOTAL budgets don't slow
% non-oscillatory relaxations (they converge and early-stop) -- they only cap
% a non-converging relaxation.
p.pde_T_chunk    = 2000;             % ~30 oscillation periods per chunk
p.pde_max_chunks = 12;               % caps the fallback at ceil(12/2)=6 chunks (~24s/cell at L=7);
                                     % non-oscillatory seed/cold-start relaxations converge and
                                     % early-stop well before this, so they are unaffected.
p.pde_T_end      = 1e5;              % seed / cold-start budget
p.pde_fallback_T = 1e4;              % Phase-2 / quarter continuation fallback budget

%% ---- rebuild flags --------------------------------------------------------
% Each stage reuses its existing output unless forced true (or the file is
% missing). Set true when parameters change or to rebuild from scratch.
force_hpss   = false;    % Stage 1: reuse saved M_P=6 HPSS (just rebuilt w/ construct-seed)  [restore to true for a full rebuild]
force_hss    = false;    % Stage 2: reuse saved M_P=6 HSS  [restore to true for a full rebuild]
force_mpss   = false;     % Stage 3: branch_MPSS_latest.mat
force_sheets = false;     % Stage 5: atlas_sheets_latest.mat
% force_trace unconditionally rebuilds the trace + classification. Rarely
% needed: Stages 8-9 auto-detect staleness, so editing any trace_* setting
% below already triggers the rebuild. Use only to force a rerun with no
% settings change (e.g. after editing run_atlas_pde_trace.m itself).
force_trace  = true;

% ---- Stage 8: PDE trajectory (part of the Stage 8/9 staleness key) ---------
% Seed sheet + grid cell for the PDE initial condition, and a TWO-PHASE linear
% time grid: step trace_dt_fine over [0, trace_t_fine_end] to resolve the fast
% early transient, then the coarser trace_dt_coarse out to trace_T_end.
trace_sheet       = 'HPSS_HSS';
trace_cell        = [40 40];   % [i_GEF, j_PAK] of the left half domain
trace_T_end       = 2000;
trace_t_fine_end  = 2000;      % end of the fine-dt regime
trace_dt_fine     = 1;
trace_dt_coarse   = 50;
trace_N_mesh      = 200;       % spatial mesh points over [0, L]

% Classification thresholds (Stage 9, part of the staleness key). Provisional
% residual flags, not calibrated decisions (see classify_pde_to_atlas_sheets.m).
trace_ambiguity_gap = 0.10;
trace_off_atlas_tol = 0.15;

% Editing anything above triggers exactly the stage(s) that depend on it.

% ---- Stage 10 (phase-plane pages), display only ---------------------------
% How many sheets to draw a trace page for (not part of any staleness key):
%   []  -> one page per sheet the trajectory dwelled on.
%   N   -> the N best-fitting sheets by min residual, even non-matching ones.
%   Inf / 9 -> all nine sheets.
trace_n_top_sheets = 9;

% ---- Stage 12 (movie), display only ---------------------------------------
render_movie           = true;   % false skips Stage 12 (slowest plotting stage)
movie_fps              = 15;
movie_frame_stride     = 10;      % render every Nth stored PDE frame (raise for a fine trace)
movie_resolution_scale = 2;      % raster supersampling: 1 = 900x750 base, 3 = tripled
movie_quality          = 90;

% ---- figure display -------------------------------------------------------
% true = figure windows pop up on screen; false = fully headless. Either way
% every plotting stage still exports its figure(s) to this run's plots folder.
show_figures = false;

% Convert the boolean to the 'on'/'off' string MATLAB's 'Visible' expects.
% (Lookup-table form, not if/else, to avoid a checkcode dead-branch warning.)
vis_map = {'off', 'on'};
fig_visible = vis_map{show_figures + 1};

%% ---- per-parameter-set results folder ------------------------------------
% Each model parameter set (the block at the top) gets its own folder
% data/<params>/ with outputs/ and plots/ subfolders, so different parameter
% sets never overwrite each other. Re-running the same parameters reuses the
% same folder. Every stage below reads/writes inside it.
run_name = atlas_run_name(p);
run_dir  = fullfile(pipeline_dir, 'data', run_name);
out_dir  = fullfile(run_dir, 'outputs');
plot_dir = fullfile(run_dir, 'plots');
if ~exist(out_dir, 'dir'),  mkdir(out_dir);  end
if ~exist(plot_dir, 'dir'), mkdir(plot_dir); end
p.out_dir = out_dir;   % consumed by the sweep + branch build stages

hpss_sweep_path  = fullfile(out_dir, 'bvp_sweep_full_results.mat');
hss_sweep_path   = fullfile(out_dir, 'hss_sweep_multibranch_results.mat');
branch_hpss_path = fullfile(out_dir, 'branch_HPSS_latest.mat');
branch_hss_path  = fullfile(out_dir, 'branch_HSS_latest.mat');
branch_mpss_path = fullfile(out_dir, 'branch_MPSS_latest.mat');
atlas_sheets_path = fullfile(out_dir, 'atlas_sheets_latest.mat');

% ---- initial-condition-specific subfolder --------------------------------
% Every Stage 1-7/2b output above is a property of the PARAMETER SET and lives
% in the shared run_dir. Everything downstream of the PDE trajectory (Stages
% 8-12: the trace .mat, its classification .mat, the trace-overlaid phase-plane
% pages, the classification-summary PDFs, and the movie) depends on the
% trajectory's INITIAL CONDITION (trace_sheet + trace_cell). Those go into a
% per-IC subfolder data/<params>/<sheet>_<i>_<j>/ so that different ICs at the
% SAME parameters never overwrite each other (e.g. HSS_HSS_25_26). Change the
% trace_sheet/trace_cell settings at the top and you get a fresh IC folder,
% while the expensive shared branch libraries/atlas are reused as-is.
ic_name = sprintf('%s_%d_%d', trace_sheet, trace_cell(1), trace_cell(2));
ic_dir  = fullfile(run_dir, ic_name);
if ~exist(ic_dir, 'dir'), mkdir(ic_dir); end

% Per-stage plotting options. 'visible' obeys show_figures; the branch/atlas/
% save paths all point into this run's folder. Two display-only toggles
% (default off):
%   interpolate_missing - NaN-fill grid cells from neighbours for smoother
%                         surfaces (display only; never alters the data).
%   show_markers        - mark cells that exist but are not provisionally
%                         stable (uninformative until stability is computed).
viz_opts = struct('interpolate_missing', false, 'show_markers', false, 'visible', fig_visible, ...
    'hss_path', branch_hss_path, 'hpss_path', branch_hpss_path, 'mpss_path', branch_mpss_path, ...
    'out_dir', out_dir, 'save_dir', plot_dir);                                                 % Stage 4
plot_opts = struct('show_markers', false, 'visible', fig_visible, ...
    'atlas_path', atlas_sheets_path, 'save_dir', plot_dir);                                     % Stage 6
combined_opts = struct('visible', fig_visible, ...
    'hss_path', branch_hss_path, 'hpss_path', branch_hpss_path, 'mpss_path', branch_mpss_path, ...
    'atlas_path', atlas_sheets_path, ...
    'save_path', fullfile(plot_dir, 'Atlas_Surfaces_Combined.pdf'));                            % Stage 7

fprintf('\n============================================================\n');
fprintf(' Atlas surface pipeline (run_atlas_surface_pipeline.m)\n');
fprintf(' results folder:        %s\n', run_dir);
fprintf(' IC-specific subfolder: %s   (Stages 8-12)\n', ic_dir);
fprintf('============================================================\n');

%% ---- Stage 1: HPSS branch library ---------------------------------------
if force_hpss || ~exist(branch_hpss_path, 'file')
    fprintf('\n[HPSS] building (force=%d, exists=%d)...\n', force_hpss, exist(branch_hpss_path,'file')>0);
    run_hybrid_sweeper_pdefb(p);                       % -> bvp_sweep_full_results.mat
    build_branch_library(hpss_sweep_path, hss_sweep_path, out_dir);
else
    fprintf('\n[HPSS] using existing %s\n', branch_hpss_path);
end

%% ---- Stage 2: HSS branch library ----------------------------------------
if force_hss || ~exist(branch_hss_path, 'file')
    fprintf('\n[HSS] building (force=%d, exists=%d)...\n', force_hss, exist(branch_hss_path,'file')>0);
    sweep_hss_multibranch_new(p);                      % -> hss_sweep_multibranch_results.mat
    build_branch_library(hpss_sweep_path, hss_sweep_path, out_dir);
else
    fprintf('[HSS]  using existing %s\n', branch_hss_path);
end

%% ---- Stage 2b: HSS dispersion-relation stability atlas -------------------
% Linear stability of the per-cell homogeneous steady state by dispersion
% relation: WHERE the HSS is Turing-unstable across (Dn_PAK, Dn_GEF), and the
% fastest-growing ONSET wavelength expressed as spatial periods on the
% half-domain (0.5 = one half-peak = HPSS morphology). Onset only, NOT the
% coarsened attractor peak count. Display-only, cheap; reads branch_HSS.
fprintf('\n=== HSS dispersion-relation stability atlas ===\n');
if exist(branch_hss_path, 'file')
    plot_hss_dispersion_atlas(struct('hss_path', branch_hss_path, ...
        'save_dir', plot_dir, 'domain', 'half', 'visible', fig_visible));
else
    fprintf('  No HSS branch file present; nothing to analyse (Stage 2 must run first).\n');
end

%% ---- Stage 3: MPSS branch library ---------------------------------------
if force_mpss || ~exist(branch_mpss_path, 'file')
    fprintf('\n[MPSS] building via continuation (force=%d, exists=%d)...\n', force_mpss, exist(branch_mpss_path,'file')>0);
    build_mpss_branch(p, out_dir);                     % -> branch_MPSS_latest.mat
else
    fprintf('[MPSS] using existing %s\n', branch_mpss_path);
end

%% ---- provenance + grid status -------------------------------------------
fprintf('\n=== loaded branch files ===\n');
print_provenance('HSS ', branch_hss_path);
print_provenance('HPSS', branch_hpss_path);
print_provenance('MPSS', branch_mpss_path);

%% ---- Stage 4: surface visualization -------------------------------------
fprintf('\n=== visualizing atlas surfaces ===\n');
if exist(branch_hss_path,'file') || exist(branch_hpss_path,'file') || exist(branch_mpss_path,'file')
    visualize_atlas_surfaces(viz_opts);
else
    fprintf('  No branch files present; nothing to visualize.\n');
end

%% ---- Stage 5: atlas sheet construction ----------------------------------
if force_sheets || ~exist(atlas_sheets_path, 'file')
    fprintf('\n[SHEETS] building atlas sheets (force=%d, exists=%d)...\n', ...
        force_sheets, exist(atlas_sheets_path,'file')>0);
    build_atlas_sheets(branch_hpss_path, branch_hss_path, branch_mpss_path, out_dir);
else
    fprintf('\n[SHEETS] using existing %s\n', atlas_sheets_path);
end
print_sheets_provenance(atlas_sheets_path);

%% ---- Stage 6: atlas sheet plotting ---------------------------------------
fprintf('\n=== plotting atlas sheets ===\n');
if exist(atlas_sheets_path, 'file')
    plot_atlas_sheets(plot_opts);
else
    fprintf('  No atlas sheets file present; nothing to plot (Stage 5 must run first).\n');
end

%% ---- Stage 7: combined 3-D surfaces + all-sheets nullclines ---------------
% One figure: all three branch flux surfaces overlaid, with every sheet's
% nullcline draped on top at its exact height.
fprintf('\n=== plotting combined atlas surfaces ===\n');
if exist(branch_hss_path,'file') && exist(branch_hpss_path,'file') && exist(branch_mpss_path,'file')
    plot_atlas_surfaces_combined(combined_opts);
else
    fprintf('  Not all branch files present; nothing to plot.\n');
end

%% ---- Stage 8: atlas-seeded PDE trajectory --------------------------------
% Assembles the two-phase time grid, packs the trajectory-defining settings +
% source branch/atlas ids into trace_cfg, and runs the PDE. The trace is
% regenerated only if force_trace, the file is missing, or trace_cfg differs
% from the run_config stored in the existing trace (decide_trace).
trace_path = fullfile(ic_dir, 'atlas_pde_trace_latest.mat');
trace_t_span = unique([0:trace_dt_fine:trace_t_fine_end, ...
                       trace_t_fine_end:trace_dt_coarse:trace_T_end, trace_T_end]);
trace_cfg = struct( ...
    'sheet',   trace_sheet, ...
    'cell',    trace_cell, ...
    't_span',  trace_t_span, ...
    'N_mesh',  trace_N_mesh, ...
    'src_ids', struct('atlas', file_id(atlas_sheets_path), ...
                      'HSS',   file_id(branch_hss_path), ...
                      'HPSS',  file_id(branch_hpss_path), ...
                      'MPSS',  file_id(branch_mpss_path)));

fprintf('\n=== Stage 8: PDE trajectory trace ===\n');
[need_trace, trace_reason] = decide_trace(trace_path, trace_cfg, force_trace);
if need_trace
    fprintf(['[TRACE] regenerating (%s): sheet=%s cell=[%d,%d] T_end=%g\n', ...
             '        t_span: dt=%g on [0,%g], dt=%g on [%g,%g] -> %d frames...\n'], ...
        trace_reason, trace_sheet, trace_cell(1), trace_cell(2), trace_T_end, ...
        trace_dt_fine, trace_t_fine_end, trace_dt_coarse, trace_t_fine_end, trace_T_end, ...
        numel(trace_t_span));
    run_atlas_pde_trace(struct('sheet_name', trace_sheet, 'cell', trace_cell, ...
        't_span', trace_t_span, 'N_mesh', trace_N_mesh, ...
        'run_config', trace_cfg, 'visible', fig_visible, 'out_dir', ic_dir, ...
        'atlas_path', atlas_sheets_path, 'hss_path', branch_hss_path, ...
        'hpss_path', branch_hpss_path, 'mpss_path', branch_mpss_path));
    trace_ran = true;
else
    fprintf('[TRACE] reusing existing %s (%s)\n', trace_path, trace_reason);
    trace_ran = false;
end

%% ---- Stage 9: PDE-to-sheet classification --------------------------------
% Reclassify only if the trace was regenerated, the classification is missing
% or was built from a different trace, or the thresholds changed (decide_class).
class_path = fullfile(ic_dir, 'pde_atlas_classification_latest.mat');
fprintf('\n=== Stage 9: classifying trace against all 9 sheets ===\n');
[need_class, class_reason] = decide_class(class_path, trace_path, trace_ran, ...
    trace_ambiguity_gap, trace_off_atlas_tol);
if need_class
    fprintf('[CLASSIFY] regenerating (%s)...\n', class_reason);
    classify_pde_to_atlas_sheets(struct('trace_path', trace_path, ...
        'ambiguity_gap', trace_ambiguity_gap, 'off_atlas_tol', trace_off_atlas_tol, ...
        'out_dir', ic_dir, 'atlas_path', atlas_sheets_path, 'hss_path', branch_hss_path, ...
        'hpss_path', branch_hpss_path, 'mpss_path', branch_mpss_path));
else
    fprintf('[CLASSIFY] reusing existing %s (%s)\n', class_path, class_reason);
end

%% ---- Stage 10: atlas phase-plane plots -----------------------------------
% One page per selected sheet (see trace_n_top_sheets): that sheet's flux
% field/nullclines with the full trajectory overlaid (solid where on this
% sheet, grey dotted elsewhere).
fprintf('\n=== Stage 10: plotting atlas phase plane(s) ===\n');
if exist(atlas_sheets_path, 'file')
    plot_atlas_phase_plane(struct('trace_path', trace_path, 'class_path', class_path, ...
        'n_top', trace_n_top_sheets, 'visible', fig_visible, ...
        'atlas_path', atlas_sheets_path, 'save_dir', ic_dir));
else
    fprintf('  No atlas sheets file present; nothing to plot (Stage 5 must run first).\n');
end

%% ---- Stage 11: PDE-to-atlas classification summary -----------------------
% The classifier's own diagnostic view of the Stage 9 table: best sheet per
% frame over time, and all-sheet residual vs time, plus a frame-vs-best-sheet
% profile overlay figure.
fprintf('\n=== Stage 11: plotting PDE-to-atlas classification summary ===\n');
if exist(class_path, 'file')
    plot_pde_to_atlas_classification(struct('class_path', class_path, ...
        'trace_path', trace_path, 'visible', fig_visible, 'save_dir', ic_dir, ...
        'atlas_path', atlas_sheets_path, 'hss_path', branch_hss_path, ...
        'hpss_path', branch_hpss_path, 'mpss_path', branch_mpss_path));
else
    fprintf('  No classification file present; nothing to plot (Stage 9 must run first).\n');
end

%% ---- Stage 12: raw PDE trajectory movie -----------------------------------
% Renders the trace's raw spatial state to .mp4 (u,q and v,w,p vs x over time).
% render_movie is passed as an argument, not tested in a top-level if, to avoid
% a checkcode dead-branch warning.
run_stage12_movie(render_movie, trace_path, movie_fps, movie_frame_stride, ...
    movie_resolution_scale, movie_quality, fig_visible, ic_dir);

fprintf('\n============================================================\n');
fprintf(' Done.\n');
fprintf('============================================================\n');

%% ---- local helper --------------------------------------------------------
function name = atlas_run_name(p)
    % Filesystem-safe folder name reflecting the model parameter set, so each
    % parameter set gets its own data/<name>/ results folder. Decimal points
    % become 'p' (0.1 -> 0p1) to avoid dots in directory names.
    f = @(x) strrep(num2str(x, '%g'), '.', 'p');
    name = sprintf('a%s_b%s_c%s_d%s_e%s_f%s_Dm%s_Dc%s_L%s_MBG%s_MP%s_res%d', ...
        f(p.a), f(p.b), f(p.c), f(p.d), f(p.e), f(p.f), ...
        f(p.Dm), f(p.Dc), f(p.L), f(p.M_BG), f(p.M_P), p.res);
end

function print_provenance(label, path)
    if ~exist(path, 'file')
        fprintf('  [%s] MISSING: %s\n', label, path);
        return;
    end
    try
        d = load(path, 'computation_id', 'res', 'M_BG', 'M_P', 'L_total', 'existence_mask');
        fprintf('  [%s] %s\n         id=%s res=%d M_BG=%g M_P=%g L=%g valid=%d/%d\n', ...
            label, path, char(d.computation_id), d.res, d.M_BG, d.M_P, d.L_total, ...
            sum(d.existence_mask(:)), numel(d.existence_mask));
    catch err
        fprintf('  [%s] %s (could not read provenance: %s)\n', label, path, err.message);
    end
end

function print_sheets_provenance(path)
    % atlas_sheets_latest.mat has a different top-level schema (nine named
    % sub-sheets, no flat existence_mask) so this mirrors print_provenance's
    % style rather than reusing it directly.
    if ~exist(path, 'file')
        fprintf('  [SHEETS] MISSING: %s\n', path);
        return;
    end
    try
        d = load(path, 'computation_id', 'res', 'sheet_names', ...
            'source_HSS_id', 'source_HPSS_id', 'source_MPSS_id', 'sheets');
        fprintf('  [SHEETS] %s\n         id=%s res=%d  sources: HSS=%s HPSS=%s MPSS=%s\n', ...
            path, char(d.computation_id), d.res, ...
            char(d.source_HSS_id), char(d.source_HPSS_id), char(d.source_MPSS_id));
        for k = 1:numel(d.sheet_names)
            nm = d.sheet_names{k};
            ex = d.sheets.(nm).sheet_exists;
            fprintf('           %-12s sheet_exists=%d/%d\n', nm, sum(ex(:)), numel(ex));
        end
    catch err
        fprintf('  [SHEETS] %s (could not read provenance: %s)\n', path, err.message);
    end
end

function id = file_id(path)
    % computation_id of a *_latest.mat, or '' if missing/unreadable. Used to
    % fold upstream provenance into the trace staleness key so rebuilding a
    % branch/atlas file invalidates a trace seeded from it.
    id = '';
    if ~exist(path, 'file'), return; end
    try
        d = load(path, 'computation_id');
        if isfield(d, 'computation_id'), id = char(d.computation_id); end
    catch
        id = '';
    end
end

function [need, reason] = decide_trace(trace_path, cfg, force)
    % Rerun the PDE trace iff forced / missing / no stored key / key changed.
    if force
        need = true; reason = 'force_trace=true'; return;
    end
    if ~exist(trace_path, 'file')
        need = true; reason = 'no existing trace file'; return;
    end
    try
        S = load(trace_path, 'run_config');
    catch
        need = true; reason = 'existing trace unreadable'; return;
    end
    if ~isfield(S, 'run_config') || ~isstruct(S.run_config) || isempty(fieldnames(S.run_config))
        need = true; reason = 'existing trace has no stored settings (older file)'; return;
    end
    if ~isequal(S.run_config, cfg)
        need = true; reason = 'trace settings changed'; return;
    end
    need = false; reason = 'settings unchanged';
end

function [need, reason] = decide_class(class_path, trace_path, trace_ran, amb, offtol)
    % Reclassify iff the trace was just regenerated / class file missing /
    % class was built from a different trace / thresholds changed.
    if trace_ran
        need = true; reason = 'trace was regenerated'; return;
    end
    if ~exist(class_path, 'file')
        need = true; reason = 'no existing classification'; return;
    end
    try
        Cp = load(class_path, 'provenance');
        Tr = load(trace_path, 'computation_id');
        pv = Cp.provenance;
    catch
        need = true; reason = 'classification/trace provenance unreadable'; return;
    end
    if ~isfield(pv, 'trace_id') || ~strcmp(char(pv.trace_id), char(Tr.computation_id))
        need = true; reason = 'classification built from a different trace'; return;
    end
    if ~(isfield(pv,'ambiguity_gap') && isequal(pv.ambiguity_gap, amb) && ...
         isfield(pv,'off_atlas_tol') && isequal(pv.off_atlas_tol, offtol))
        need = true; reason = 'classification thresholds changed'; return;
    end
    need = false; reason = 'trace + thresholds unchanged';
end

function run_stage12_movie(do_render, trace_path, fps, frame_stride, resolution_scale, quality, vis, save_dir)
    fprintf('\n=== Stage 12: rendering PDE trajectory movie ===\n');
    if do_render
        if exist(trace_path, 'file')
            plot_atlas_pde_movie(struct('trace_path', trace_path, 'fps', fps, ...
                'frame_stride', frame_stride, 'resolution_scale', resolution_scale, ...
                'quality', quality, 'visible', vis, 'save_dir', save_dir));
        else
            fprintf('  No trace file present; nothing to render (Stage 8 must run first).\n');
        end
    else
        fprintf('  Skipped (render_movie = false).\n');
    end
end
