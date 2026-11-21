function manifest = atlas_build_libraries(model, par, res, out_dir, opts)
%ATLAS_BUILD_LIBRARIES  Build or reuse the half-cell libraries (starting profiles; Stages 1-3 only).
%
%   manifest = ATLAS_BUILD_LIBRARIES(model, par, res)
%   manifest = ATLAS_BUILD_LIBRARIES(model, par, res, out_dir, opts)
%   d = ATLAS_BUILD_LIBRARIES('default_dir', model, par, res)
%   tf = ATLAS_BUILD_LIBRARIES('present', model, out_dir)   all three
%       branch_{HSS,HPSS,MPSS}_latest.mat files exist (no-peak, half-peak, full-peak)
%
%   The library stages of run_atlas_surface_pipeline.m, nothing else (never
%   Stage 2b or Stages 4-12, so the script's force_trace does not apply).
%   The libraries are the isolated half-cell steady states the paper uses as starting profiles.
%   model   'original_polarity' | 'guzzo_noxa'
%   par     the dimensional parameters: original_polarity a b c d e f Dm Dc
%           L M_BG M_P (the script's expressions, e.g. 2.38/4, give the same
%           doubles as the stored libraries); guzzo_noxa the gnx_model
%           'parameters' fields without res (L r_speed dba_factor rho_R rho_B
%           dn_half_range)
%   res     the library grid (the script's p.res; noxa: P.res)
%   out_dir the library folder; '' or omitted = the script's own folder:
%           <repo>/data/<atlas_run_name(p)>/outputs (original_polarity) or
%           <repo>/data/<gnx_run_name(P)>/outputs (guzzo_noxa)
%   opts    force: cellstr of stages to rebuild even when present ('HPSS',
%           'HSS', 'MPSS'; default {}); confirmation: the noxa pre-seeding
%           confirmation file put into the manifest (optional)
%           Higher-mode libraries:
%           max_mode: build the higher-mode libraries n = 3 .. max_mode (n monotone runs; more
%             peaks than the paper's three shapes)
%             (build_mode_branch / gnx_build_mode) after the three, stopping
%             after the first mode whose mask is empty (its file is still
%             written, run_status 'empty'); default 2 = build nothing new.
%             A mode whose branch_MODE<n>_latest.mat exists is reused (opts.force
%             naming 'MODE<n>' rebuilds it).
%           max_reroot_tries: passed to the mode builder (default: the
%             sweeper's res^2, unchanged)
%           pool_workers: CREATE the parallel pool with parpool(pool_workers)
%             before any build; refused when a pool of another size is
%             already open (the size used is in manifest.pool_workers_used)
%           out_dir is required (nonempty) when max_mode > 2
%           mode_opts: options passed to the noxa mode builder (tests)
%
%   original_polarity: a call-for-call transcription of the script's Stages
%   1-3 with the script's p (par, res, its solver-option block in its field
%   order, then out_dir): Stage 1 run_hybrid_sweeper_pdefb(p) +
%   build_branch_library, Stage 2 sweep_hss_multibranch_new(p) +
%   build_branch_library, Stage 3 build_mpss_branch(p, out_dir); each stage
%   runs only when its branch_*_latest.mat is missing or opts.force names it.
%   p.use_parallel follows the script (true): a build opens a parallel pool.
%   guzzo_noxa: gnx_run_libraries(P, repo, struct('out_dir', out_dir,
%   'stages', <missing stages>)) with P = gnx_model('parameters', par + res);
%   stages whose files exist are skipped.
%   Before reuse, the three library files must exist and their stored par
%   must equal the run parameters (the coupled loader checks this again).
%   Returns the coupled pipeline's library manifest struct('directory',
%   out_dir) (+ 'confirmation'); prints what it built or reused.

if ischar(model) && strcmp(model, 'default_dir')
    manifest = default_dir(par, res, out_dir);
    return
end
if ischar(model) && strcmp(model, 'present')
    manifest = all(cellfun(@(t) exist(fullfile(res, sprintf('branch_%s_latest.mat', t)), 'file') == 2, ...
        {'HSS', 'HPSS', 'MPSS'}));
    return
end
if nargin < 4, out_dir = ''; end
if nargin < 5 || isempty(opts), opts = struct(); end
force = {};
if isfield(opts, 'force') && ~isempty(opts.force), force = cellstr(opts.force); end
max_mode = 2;
if isfield(opts, 'max_mode') && ~isempty(opts.max_mode), max_mode = opts.max_mode; end
if ~(isscalar(max_mode) && max_mode == fix(max_mode) && max_mode >= 2)
    error('ATLAS_BUILD_LIBRARIES:MAX_MODE', 'opts.max_mode must be an integer >= 2.');
end
if max_mode > 2 && isempty(out_dir)
    error('ATLAS_BUILD_LIBRARIES:OUT_DIR', ['opts.max_mode > 2 needs an explicit out_dir (the ' ...
        'builder''s default folder is never used for mode builds).']);
end
pool_used = setup_pool(opts);
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
if isempty(out_dir), out_dir = default_dir(model, par, res); end
types = {'HSS', 'HPSS', 'MPSS'};
present = cellfun(@(t) exist(fullfile(out_dir, sprintf('branch_%s_latest.mat', t)), 'file') == 2, types);
switch model
    case 'original_polarity'
        p = script_p(par, res);
        p.out_dir = out_dir;
        hpss_sweep_path  = fullfile(out_dir, 'bvp_sweep_full_results.mat');
        hss_sweep_path   = fullfile(out_dir, 'hss_sweep_multibranch_results.mat');
        branch_hpss_path = fullfile(out_dir, 'branch_HPSS_latest.mat');
        branch_hss_path  = fullfile(out_dir, 'branch_HSS_latest.mat');
        branch_mpss_path = fullfile(out_dir, 'branch_MPSS_latest.mat');
        force_hpss = any(strcmp(force, 'HPSS'));
        force_hss  = any(strcmp(force, 'HSS'));
        force_mpss = any(strcmp(force, 'MPSS'));
        if all(present) && ~force_hpss && ~force_hss && ~force_mpss
            check_reuse(model, out_dir, p);
        end
        if ~exist(out_dir, 'dir') && (force_hpss || force_hss || force_mpss || ~all(present))
            mkdir(out_dir);
        end
        if force_hpss || ~exist(branch_hpss_path, 'file')
            fprintf('\n[HPSS] building (force=%d, exists=%d)...\n', force_hpss, exist(branch_hpss_path,'file')>0);
            run_hybrid_sweeper_pdefb(p);
            build_branch_library(hpss_sweep_path, hss_sweep_path, out_dir);
        else
            fprintf('\n[HPSS] using existing %s\n', branch_hpss_path);
        end
        if force_hss || ~exist(branch_hss_path, 'file')
            fprintf('\n[HSS] building (force=%d, exists=%d)...\n', force_hss, exist(branch_hss_path,'file')>0);
            sweep_hss_multibranch_new(p);
            build_branch_library(hpss_sweep_path, hss_sweep_path, out_dir);
        else
            fprintf('[HSS]  using existing %s\n', branch_hss_path);
        end
        if force_mpss || ~exist(branch_mpss_path, 'file')
            fprintf('\n[MPSS] building via continuation (force=%d, exists=%d)...\n', force_mpss, exist(branch_mpss_path,'file')>0);
            build_mpss_branch(p, out_dir);
        else
            fprintf('[MPSS] using existing %s\n', branch_mpss_path);
        end
    case 'guzzo_noxa'
        P = noxa_P(par, res);
        missing = types(~present | ismember(types, force));
        if isempty(missing)
            check_reuse(model, out_dir, P);
            fprintf('[guzzo_noxa] using existing libraries in %s\n', out_dir);
        else
            fprintf('[guzzo_noxa] building %s in %s\n', strjoin(missing, ', '), out_dir);
            gnx_run_libraries(P, repo, struct('out_dir', out_dir, 'stages', {missing}));
        end
    otherwise
        error('ATLAS_BUILD_LIBRARIES:MODEL', 'Unknown model %s.', model);
end
modes = build_modes(model, par, res, out_dir, opts, max_mode, force);
manifest = struct('directory', out_dir);
if max_mode > 2 || ~isnan(pool_used)
    manifest.modes = modes;
    manifest.pool_workers_used = pool_used;
end
if isfield(opts, 'confirmation') && ~isempty(opts.confirmation)
    manifest.confirmation = opts.confirmation;
end
end

% =========================================================================

function p = script_p(par, res)
% the script's p: the parameter block in its order, res, then its
% solver-option block verbatim (run_atlas_surface_pipeline.m)
p = struct();
for f = {'a', 'b', 'c', 'd', 'e', 'f', 'Dm', 'Dc', 'L', 'M_BG', 'M_P'}
    p.(f{1}) = par.(f{1});
end
p.res = res;
p.seeding_mode   = 'auto';
p.use_parallel   = true;
p.skip_stability = true;
p.mpss_strategy  = 'continuation';
p.phase0_construct = true;
p.osc_fallback   = true;
p.pde_T_chunk    = 2000;
p.pde_max_chunks = 12;
p.pde_T_end      = 1e5;
p.pde_fallback_T = 1e4;
end

function used = setup_pool(opts)
% opts.pool_workers: create the pool (parpool(n)); refuse a pool of another
% size. opts.pool_fns (tests): struct(current, create) replacing
% gcp('nocreate') and parpool.
used = NaN;
if ~isfield(opts, 'pool_workers') || isempty(opts.pool_workers), return; end
n = opts.pool_workers;
if ~(isscalar(n) && n == fix(n) && n >= 1)
    error('ATLAS_BUILD_LIBRARIES:POOL', 'opts.pool_workers must be a positive integer.');
end
cur = @() gcp('nocreate'); make = @parpool;
if isfield(opts, 'pool_fns')
    cur = opts.pool_fns.current; make = opts.pool_fns.create;
end
pool = cur();
if ~isempty(pool)
    if pool.NumWorkers ~= n
        error('ATLAS_BUILD_LIBRARIES:POOL', ['A parallel pool of %d workers is already open; ' ...
            'opts.pool_workers = %d (close it first).'], pool.NumWorkers, n);
    end
else
    pool = make(n);
end
used = pool.NumWorkers;
fprintf('[pool] %d worker(s)\n', used);
end

function modes = build_modes(model, par, res, out_dir, opts, max_mode, force)
% modes 3 .. max_mode, stopping after the first empty one
modes = struct('n', {}, 'path', {}, 'run_status', {}, 'n_valid', {}, 'built', {});
for n = 3:max_mode
    f = fullfile(out_dir, sprintf('branch_MODE%d_latest.mat', n));
    if exist(f, 'file') == 2 && ~any(strcmp(force, sprintf('MODE%d', n)))
        Z = load(f, 'run_status', 'existence_mask');
        fprintf('[MODE%d] using existing %s (run_status %s)\n', n, f, Z.run_status);
        modes(end+1) = struct('n', n, 'path', f, 'run_status', Z.run_status, ...
            'n_valid', nnz(Z.existence_mask), 'built', false); %#ok<AGROW>
    else
        fprintf('\n[MODE%d] building\n', n);
        bo = struct();
        if isfield(opts, 'max_reroot_tries'), bo.max_reroot_tries = opts.max_reroot_tries; end
        switch model
            case 'original_polarity'
                p = script_p(par, res);
                p.out_dir = out_dir;
                R = build_mode_branch(p, out_dir, n, bo);
            case 'guzzo_noxa'
                mo = struct();
                if isfield(opts, 'mode_opts'), mo = opts.mode_opts; end
                R = gnx_build_mode(noxa_P(par, res), n, out_dir, mo);
        end
        modes(end+1) = struct('n', n, 'path', R.path, 'run_status', R.run_status, ...
            'n_valid', R.n_valid, 'built', true); %#ok<AGROW>
    end
    if strcmp(modes(end).run_status, 'empty')
        fprintf('[MODE%d] empty: no higher mode is built (the stop rule is a cost rule).\n', n);
        break
    end
end
end

function P = noxa_P(par, res)
s = par;
s.res = res;
P = gnx_model('parameters', s);
end

function d = default_dir(model, par, res)
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
switch model
    case 'original_polarity'
        d = fullfile(repo, 'data', atlas_run_name(script_p(par, res)), 'outputs');
    case 'guzzo_noxa'
        d = fullfile(repo, 'data', gnx_run_name(noxa_P(par, res)), 'outputs');
    otherwise
        error('ATLAS_BUILD_LIBRARIES:MODEL', 'Unknown model %s.', model);
end
end

function check_reuse(model, out_dir, p)
switch model
    case 'original_polarity'
        names = {'a', 'b', 'c', 'd', 'e', 'f', 'Dm', 'Dc', 'L', 'M_BG', 'M_P'};
    otherwise
        names = {'kR', 'kRB', 'kB', 'kBB', 'dR', 'dB', 'dBA', 'K', 'KA', 'DM', 'DC', 'r_speed', ...
            'dba_factor', 'L', 'rho_R', 'rho_B'};
end
for t = {'HSS', 'HPSS', 'MPSS'}
    f = fullfile(out_dir, sprintf('branch_%s_latest.mat', t{1}));
    Z = load(f, 'par');
    for k = 1:numel(names)
        if ~isfield(Z.par, names{k}) || ~isequal(Z.par.(names{k}), p.(names{k}))
            error('ATLAS_BUILD_LIBRARIES:PAR', ['%s was built at %s = %s, but the run parameters say %s: ' ...
                'another library folder, or a rebuild (opts.force).'], f, names{k}, ...
                num2str(getf(Z.par, names{k}), 17), num2str(p.(names{k}), 17));
        end
    end
end
end

function v = getf(s, f)
if isfield(s, f), v = s.(f); else, v = NaN; end
end

function name = atlas_run_name(p)
    % Filesystem-safe folder name reflecting the model parameter set, so each
    % parameter set gets its own data/<name>/ results folder. Decimal points
    % become 'p' (0.1 -> 0p1) to avoid dots in directory names.
    f = @(x) strrep(num2str(x, '%g'), '.', 'p');
    name = sprintf('a%s_b%s_c%s_d%s_e%s_f%s_Dm%s_Dc%s_L%s_MBG%s_MP%s_res%d', ...
        f(p.a), f(p.b), f(p.c), f(p.d), f(p.e), f(p.f), ...
        f(p.Dm), f(p.Dc), f(p.L), f(p.M_BG), f(p.M_P), p.res);
end
