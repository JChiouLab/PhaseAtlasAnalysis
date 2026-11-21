function report = gnx_run_libraries(P, pipeline_dir, opts)
%GNX_RUN_LIBRARIES  M. xanthus (Guzzo_noXA) stages 1-3: no-peak (HSS), half-peak (HPSS) and full-peak (MPSS)
%   half-cell libraries (starting profiles).
%
%   report = GNX_RUN_LIBRARIES(P, pipeline_dir)
%   report = GNX_RUN_LIBRARIES(P, pipeline_dir, opts)
%   P = gnx_model('parameters', ...). Writes branch_{HSS,HPSS,MPSS}_latest.mat
%   (gnx_write_library) to <pipeline_dir>/data/<gnx_run_name(P)>/outputs, or
%   to opts.out_dir when given. opts (optional): out_dir, stages (cell, any of
%   'HSS','HPSS','MPSS'; default all), sweep (struct passed to gnx_sweep_hpss,
%   e.g. max_calls, solve_opts). Checkpoints go next to the libraries.
%   Calls this run made (gnx_solve counter) are printed as GNX_CALLS=<n>.
%   Called from the guzzo_noxa branch of run_atlas_surface_pipeline.m.

if nargin < 3 || isempty(opts), opts = struct(); end
[name, sha8] = gnx_run_name(P);
if isfield(opts, 'out_dir') && ~isempty(opts.out_dir)
    out_dir = opts.out_dir;
else
    out_dir = fullfile(pipeline_dir, 'data', name, 'outputs');
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
stages = {'HSS', 'HPSS', 'MPSS'};
if isfield(opts, 'stages'), stages = opts.stages; end
sweep = struct();
if isfield(opts, 'sweep'), sweep = opts.sweep; end

here = fileparts(mfilename('fullpath'));
meta = struct();
meta.run_id = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
meta.run_name = name;
meta.sha8 = sha8;
meta.model_version = gnx_model('version', {here});
meta.matlab_version = version;
meta.stages = {stages};
calls0 = gnx_solve('count');
report = struct('out_dir', out_dir, 'run_name', name, 'meta', meta, 'paths', struct());
fprintf('\n[guzzo_noxa] run %s -> %s\n', name, out_dir);

if any(strcmp(stages, 'HSS'))
    H = gnx_sweep_hss(P);
    report.paths.HSS = gnx_write_library(out_dir, 'HSS', H, P, meta);
    fprintf('[guzzo_noxa] HSS: %d/%d cells\n', nnz(H.existence_mask), P.res^2);
end
if any(strcmp(stages, 'HPSS'))
    so = sweep; so.ckpt_path = fullfile(out_dir, sprintf('checkpoint_hpss_%s.mat', meta.run_id));
    S = gnx_sweep_hpss(P, P.L / 2, 'hpss', so);
    S.existence_mask = S.ok_grid;
    S.SOL = S.SOLS;
    report.paths.HPSS = gnx_write_library(out_dir, 'HPSS', S, P, meta);
    report.hpss = struct('n', nnz(S.ok_grid), 'calls', S.n_calls, 'nmax', S.n_nmax, ...
        'cap', S.cap_reached);
    fprintf('[guzzo_noxa] HPSS: %d/%d cells, %d calls, %d at NMax, cap %d\n', ...
        nnz(S.ok_grid), P.res^2, S.n_calls, S.n_nmax, S.cap_reached);
end
if any(strcmp(stages, 'MPSS'))
    so = sweep; so.ckpt_path = fullfile(out_dir, sprintf('checkpoint_quarter_%s.mat', meta.run_id));
    M = gnx_build_mpss(P, so);
    report.paths.MPSS = gnx_write_library(out_dir, 'MPSS', M, P, meta);
    report.mpss = struct('n', nnz(M.existence_mask), 'n_quarter', nnz(M.quarter.ok_grid), ...
        'calls', M.n_calls, 'nmax', M.n_nmax, 'cap', M.quarter.cap_reached);
    fprintf('[guzzo_noxa] MPSS: %d/%d cells (%d quarters), %d calls, %d at NMax, cap %d\n', ...
        nnz(M.existence_mask), P.res^2, nnz(M.quarter.ok_grid), M.n_calls, M.n_nmax, ...
        M.quarter.cap_reached);
end
report.calls = gnx_solve('count') - calls0;
fprintf('GNX_CALLS=%d\n', report.calls);
end
