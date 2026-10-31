function rows = sa_t30_resume_refusals(E)
%SA_T30_RESUME_REFUSALS  Milestone-3 T30 (RC-2, RC-5): resume is refused on a
%   changed scientific setting, software or MATLAB version, thread mode,
%   library or code hash; a code-hash change passes only with a journaled
%   override reason (both hash sets recorded); budget leaves may change and
%   are recorded in resume_<n>.mat; a corrupt attempt file is refused; J7b (a
%   stored check relabelled 'record') is an integrity stop on replay.
T = 'T30 resume refusals';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't30');
fail = @(dd, b, ex) sa_m3_wall(dd, ex);
box = [1 1];
B = fullfile(d, 'base');
rB = sa_m3_run(B, box, 'fail', fail);
cases = {
  'changed scientific setting (solver.reltol)', @(K) 0, struct('solver', struct('reltol', 2e-6)), 'SA:RESUME:E_SETTINGS'
  'stored software version differs', @(K) eff(K, 'software_version', 'x/0.1'), struct(), 'SA:RESUME:E_VERSION'
  'stored MATLAB version differs', @(K) eff(K, 'matlab_version', '9.9'), struct(), 'SA:RESUME:E_VERSION'
  'stored thread mode differs', @(K) eff(K, 'thread_mode', 'multi:8'), struct(), 'SA:RESUME:E_VERSION'
  'library SHA-256 differs', @(K) lib(K), struct(), 'SA:RESUME:E_LIBRARY'
  'code hash differs (no override)', @(K) code(K), struct(), 'SA:RESUME:E_CODE'
  'corrupt attempt file', @(K) corrupt(K), struct(), 'SA:RESUME:E_CORRUPT'
  'J7b: stored check relabelled ''record'' with a failing value', @(K) j7b(K), struct(), 'SA:INTEGRITY'
  'tampered stored opts', @(K) topts(K), struct(), 'SA:INTEGRITY'
  };
for c = 1:size(cases, 1)
    K = fullfile(d, sprintf('case%02d', c));
    copyfile(B, K);
    cases{c, 2}(K);
    o = cases{c, 3}; o.resume = true;
    got = '';
    try
        sa_m3_run(K, box, 'fail', fail, 'settings', o);
    catch ME
        got = ME.identifier;
    end
    rows = sa_check(rows, T, sprintf('refused: %s', cases{c, 1}), strcmp(got, cases{c, 4}), got);
end
% code-hash change with a journaled override reason
K = fullfile(d, 'override'); copyfile(B, K); code(K);
r = sa_m3_run(K, box, 'fail', fail, 'settings', struct('resume', true, ...
    'resume_code_override', 'T30 lead-approved test override'));
J = fileread(fullfile(sa_layout('state', K), 'journal.txt'));
Z = load(fullfile(sa_layout('state', K), 'resume_1.mat'));
rows = sa_check(rows, T, 'code-hash change with override: resumes; journal holds the reason and both hash sets', ...
    strcmp(r.stop_reason, 'completed') && contains(J, 'override code-hash change; reason: T30 lead-approved') && ...
    contains(J, 'stored: ') && contains(J, 'current: ') && strcmp(Z.override, 'T30 lead-approved test override'), '');
% budget leaves may change; recorded
K = fullfile(d, 'budget'); copyfile(B, K);
r = sa_m3_run(K, box, 'fail', fail, 'settings', struct('resume', true, 'budget', struct('max_calls', 777), ...
    'checkpoint', struct('every_calls', 3)));
Z = load(fullfile(sa_layout('state', K), 'resume_1.mat'));
rows = sa_check(rows, T, 'budget/checkpoint leaves may change on resume and are recorded in resume_1.mat', ...
    strcmp(r.stop_reason, 'completed') && Z.budget.max_calls == 777 && Z.budget_stored.max_calls == 5000 && ...
    Z.checkpoint.every_calls == 3 && isfield(Z, 'pid') && isfield(Z, 'time') && Z.attempt_count == rB.calls, '');
got = '';
try
    sa_m3_run(fullfile(d, 'nofolder'), box, 'fail', fail, 'settings', struct('resume', true));
catch ME
    got = ME.identifier;
end
rows = sa_check(rows, T, 'resume of a folder that is not a run is refused', strcmp(got, 'SA:RESUME:E_NO_RUN'), got);
rmpath(E.stub);
end

function eff(K, f, v)
Z = load(fullfile(sa_layout('state', K), 'settings_effective.mat'));
Z.provenance.(f) = v;
save(fullfile(sa_layout('state', K), 'settings_effective.mat'), '-struct', 'Z');
end

function lib(K)
Z = load(fullfile(sa_layout('state', K), 'settings_effective.mat'));
Z.library_files(end+1).sha256 = repmat('0', 1, 64);
save(fullfile(sa_layout('state', K), 'settings_effective.mat'), '-struct', 'Z');
end

function code(K)
Z = load(fullfile(sa_layout('state', K), 'settings_effective.mat'));
Z.code_hashes(1).sha256(1) = char('a' + ('0' == Z.code_hashes(1).sha256(1)));
save(fullfile(sa_layout('state', K), 'settings_effective.mat'), '-struct', 'Z');
end

function corrupt(K)
fid = fopen(fullfile(sa_layout('state', K), 'attempts', 'A000003.mat'), 'w'); fwrite(fid, 'not a mat file'); fclose(fid);
end

function j7b(K)
f = fullfile(sa_layout('state', K), 'attempts', 'A000002.mat');
Z = load(f);
k = strcmp({Z.validity.checks.name}, 'stub_noflux');
Z.validity.checks(k).value = 1; Z.validity.checks(k).comparison = 'record';
save(f, '-struct', 'Z');
end

function topts(K)
f = fullfile(sa_layout('state', K), 'attempts', 'A000002.mat');
Z = load(f); Z.opts.nmax = 123;
save(f, '-struct', 'Z');
end
