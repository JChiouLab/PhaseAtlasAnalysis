function rows = run_portable_tests(which_tests)
%RUN_PORTABLE_TESTS  The tests that need only the package code and the stub
%   fixtures: they run from a fresh copy of the package, with none of the
%   authors' development data.
%
%   rows = RUN_PORTABLE_TESTS                   all of them
%   rows = RUN_PORTABLE_TESTS({'sa_t01_settings', ...})   a subset
%
%   First runs setup_atlas('quiet') (the traced bvp4c copy, made when missing
%   or stale). Every test output goes under
%   <package>/data/experiments/portable_tests_<stamp>/ (git-ignored; the
%   wrapper test sa_t86_wrapper_root needs its launch folder inside the package's
%   data/experiments). E.am14_real = false (sa_t52_am14 skips the parts that read
%   the authors' development run folders). sa_t62_deval (the deval replica) runs only when this
%   MATLAB's deval mode is 'replica'; sa_t82_tmp_files skips its chflags case off macOS.
%   No test calls the real corrector on a coupled problem: pipeline runs use
%   the stub model and corrector; sa_t61p_traced_solver and sa_t102_release_selfcheck make tiny
%   bvp4c calls in the self-checks of scratch copies (about 140 in all). Prints PASS, FAIL or
%   SKIP per test and a summary line; the rows also go to a results file in
%   the scratch folder. Do not run it while a run is in flight in the same
%   package folder.

E = sa_test_env();
addpath(E.root); addpath(E.pkg); addpath(E.tests);
setup_atlas('quiet');
addpath(fullfile(E.pkg, 'solver')); addpath(fullfile(E.pkg, 'post'));
E.scratch = fullfile(E.root, 'data', 'experiments', ['portable_tests_' char(datetime('now', ...
    'Format', 'yyyyMMdd_HHmmss'))]);
E.am14_real = false;
if exist(E.scratch, 'dir') ~= 7, mkdir(E.scratch); end
Eg = sag_test_env(E.scratch);
% Tests that need the authors' development data (for example a Guzzo_noXA library under
% data/experiments) are not part of this list.
all_tests = {'sa_t01_settings', 'sa_t04_source_scan', 'sa_t11_stub_model', ...
    'sa_t24_determinism', 'sa_t25_substeps_classes', 'sa_t26_step_rule', 'sa_t27_rescue', 'sa_t28_lineages', ...
    'sa_t29_resume', 'sa_t30_resume_refusals', 'sa_t31_am1', 'sa_t38_39_budget_checkpoint', ...
    'sa_t40_ring', 'sa_t41_am7', 'sa_t42_am8', 'sa_t43_am9', 'sa_t44_settings_modes', 'sa_t45_targets', ...
    'sa_t46_mode_equivalence', 'sa_t47_pilot_report', 'sa_t48_am10', 'sa_t49_am11_rescue_ring', 'sa_t50_am12', ...
    'sa_t51_am13', 'sa_t52_am14', 'sa_t56_rescue_trim', 'sa_t57_cold_cap', 'sa_t60_reports', ...
    'sa_t61p_traced_solver', 'sa_t62_deval', 'sa_t64_cutset', 'sa_t65_replay', 'sa_t67_reports', ...
    'sa_t68_lattice_axes', 'sa_t71_traces', 'sa_t76_screen', 'sa_t77_trace_at', 'sa_t78_grid_version', ...
    'sa_t79_name_roots', 'sa_t80_render_cache', 'sa_t81_shard_repair', 'sa_t82_tmp_files', ...
    'sa_t83_start_lookups', 'sa_t86_wrapper_root', 'sa_t102_release_selfcheck', ...
    'sag_t14_module', 'sag_t20_identities', 'sag_t20b_linear_stub', 'sag_t22_classifiers'};
if nargin < 1 || isempty(which_tests), which_tests = all_tests; end
mode = sa_bvp4c_setup('deval_mode');
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
status = cell(1, numel(which_tests));
for k = 1:numel(which_tests)
    name = which_tests{k};
    if strcmp(name, 'sa_t62_deval') && ~strcmp(mode, 'replica')
        status{k} = 'SKIP';
        fprintf('[SKIP] %s  (deval mode %s: the replica is not used on this MATLAB)\n', name, mode);
        continue
    end
    if startsWith(name, 'sag_') && ~strcmp(name, 'sag_t20b_linear_stub'), Ek = Eg; else, Ek = E; end
    t0 = tic;
    d0 = pwd; p0 = path;
    try
        r = feval(name, Ek);
    catch ME
        r = sa_check(struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {}), name, ...
            'test ran without error', false, sprintf('%s: %s (%s:%d)', ME.identifier, ME.message, ...
            ME.stack(1).name, ME.stack(1).line));
    end
    cd(d0); path(p0);
    rows = [rows, r]; %#ok<AGROW>
    status{k} = tern(all([r.pass]), 'PASS', 'FAIL');
    fprintf('[%s] %s  (%.1f s)\n', status{k}, name, toc(t0));
end
lines = arrayfun(@(r) sprintf('%s | %s | %s | %s', tern(r.pass, 'PASS', 'FAIL'), r.test, r.assertion, ...
    r.detail), rows, 'UniformOutput', false);
txt = strjoin(lines, newline);
summary = sprintf(['portable tests: %d PASS, %d FAIL, %d SKIP of %d tests; %d of %d assertions passed ' ...
    '(deval mode %s; %s)'], sum(strcmp(status, 'PASS')), sum(strcmp(status, 'FAIL')), sum(strcmp(status, 'SKIP')), ...
    numel(which_tests), sum([rows.pass]), numel(rows), mode, version);
fprintf('%s\n%s\n', txt, summary);
fid = fopen(fullfile(E.scratch, 'test_results_portable.txt'), 'w');
fprintf(fid, '%s\n%s\n', txt, summary); fclose(fid);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
