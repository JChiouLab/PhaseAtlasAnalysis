function rows = sa_t86_wrapper_root(E)
%SA_T86_WRAPPER_ROOT  The wrapper runs seeded_atlas_run with the repository
%   root as the current folder. Zero bvp4c calls (stub model and corrector,
%   the base stub-run configuration). The wrapper is launched from a folder under
%   data/experiments/ (runtime.clean's 'pwd' rule would flag it):
%   (1) an absolute run_dir: the run's recorded runtime.pwd is the
%       repository root, no 'pwd' reason, and the caller's folder is current
%       again afterwards;
%   (2) a relative run_dir (resolved against the caller's folder): the
%       wrapper keeps the caller's folder, says so, and the run records the
%       'pwd' reason as before; apart from pwd, clean and that one reason,
%       its recorded runtime equals run (1)'s, and its records equal run
%       (1)'s (sa_m3_equiv);
%   (3) a resume through the wrapper of a folder stopped at a budget,
%       launched from another folder: resumable, completes, records equal
%       run (1)'s, no integrity note (the launch folder is not part of a
%       run's resume identity);
%   (4) a run whose corrector throws (the stub's kill at call 3): the error reaches the caller and the
%       caller's folder is current again;
%   (5) launched from the package folder with the package not on the path
%       (found only as the current folder): the wrapper keeps the current
%       folder, says so, and the run completes.
T = 'T86 wrapper runs from the repository root';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d0 = sa_m3_testdir(E, 't86');
here0 = fullfile(d0, 'launch'); mkdir(here0);
other = fullfile(d0, 'elsewhere'); mkdir(other);
old = pwd;
back = onCleanup(@() cd(old));
s0 = sa_m3_stub('settings', 'x', [1 1]);
ov = rmfield(s0, {'model', 'run_dir', 'budget'});
base = struct('model', 'stub_two_pool', 'par', sa_m3_stub('params'), 'mode', 'smoke', ...
    'library_manifest', struct(), 'settings_overrides', ov, 'budget', struct('max_calls', 5000), ...
    'make_plot', false);

% (1) absolute run_dir
cd(here0);
c = base; c.run_dir = fullfile(d0, 'abs'); c.corrector = sa_m3_stub('corrector', sa_m3_stub('state'));
out1 = evalc('W1 = seeded_atlas_wrapper(c);');
p1 = pwd;
R1 = runtime_of(c.run_dir);
s1 = first_line(c.run_dir);
ok1 = strcmp(W1.result.stop_reason, 'completed') && strcmp(p1, here0) && strcmp(R1.pwd, E.root) && ...
    ~any(startsWith(R1.clean_reasons, 'pwd ')) && ~contains(s1, 'pwd ') && ~contains(out1, 'relative path');
rows = sa_check(rows, T, ['launched from a folder under data/experiments with an absolute run_dir: the run''s ' ...
    'runtime.pwd is the repository root, no pwd reason (summary line too), the caller''s folder is current again'], ...
    ok1, sprintf('pwd after %s | runtime.pwd %s | %s', p1, R1.pwd, s1));

% (2) relative run_dir: the caller's folder is kept
cd(here0);
c2 = base; c2.run_dir = 'rel'; c2.corrector = sa_m3_stub('corrector', sa_m3_stub('state'));
out2 = evalc('W2 = seeded_atlas_wrapper(c2);');
p2 = pwd;
R2 = runtime_of(fullfile(here0, 'rel'));
pr = R2.clean_reasons(startsWith(R2.clean_reasons, 'pwd '));
f = {'pwd', 'clean', 'clean_reasons'};
same_rest = isequaln(rmfield(R1, f), rmfield(R2, f)) && ...
    isequal(R1.clean_reasons, R2.clean_reasons(~startsWith(R2.clean_reasons, 'pwd ')));
D2 = sa_m3_equiv(c.run_dir, fullfile(here0, 'rel'));
ok2 = strcmp(W2.result.stop_reason, 'completed') && strcmp(p2, here0) && strcmp(R2.pwd, here0) && ...
    isequal(pr, {['pwd ' here0]}) && contains(out2, 'relative path') && same_rest && D2.pass;
rows = sa_check(rows, T, ['a relative run_dir: the caller''s folder is kept (one line says so), the pwd reason is ' ...
    'recorded as before; apart from pwd, clean and that reason the runtime record equals the absolute run''s, ' ...
    'and so do the records (sa_m3_equiv)'], ok2, sprintf('pwd after %s | reasons %s | same rest %d | %s', p2, ...
    strjoin(R2.clean_reasons, '; '), same_rest, strjoin(D2.diffs, ' | ')));

% (3) a resume through the wrapper, launched from another folder
cd(here0);
c3 = base; c3.run_dir = fullfile(d0, 'resumed'); c3.corrector = sa_m3_stub('corrector', sa_m3_stub('state'));
c3.budget = struct('max_calls', 5);
evalc('W3a = seeded_atlas_wrapper(c3);');
cd(other);
c3.budget = struct('max_calls', 5000);
c3.corrector = sa_m3_stub('corrector', sa_m3_stub('state'));
evalc('W3 = seeded_atlas_wrapper(c3);');
p3 = pwd;
J = strsplit(strtrim(fileread(fullfile(sa_layout('state', c3.run_dir), 'journal.txt'))), newline);
D3 = sa_m3_equiv(c.run_dir, c3.run_dir);
ok3 = strcmp(W3a.result.stop_reason, 'budget_calls') && strcmp(W3.status.state, 'resumable') && ...
    strcmp(W3.result.stop_reason, 'completed') && strcmp(p3, other) && D3.pass && ...
    ~any(startsWith(J, 'note integrity')) && any(startsWith(J, 'note resume 1'));
rows = sa_check(rows, T, ['a folder stopped at its budget, resumed through the wrapper from another folder: ' ...
    'resumable, completes, records equal the uninterrupted run''s, no integrity note, the caller''s folder is ' ...
    'current again'], ok3, sprintf('%s -> %s %s | pwd after %s | %s', W3a.result.stop_reason, W3.status.state, ...
    W3.result.stop_reason, p3, strjoin(D3.diffs, ' | ')));

% (4) an error inside the run
cd(here0);
c4 = base; c4.run_dir = fullfile(d0, 'error');
c4.corrector = sa_m3_stub('corrector', sa_m3_stub('state', 'kill_at', 3));   % throws SA:TEST:KILL at call 3
id = '';
try
    evalc('seeded_atlas_wrapper(c4);');
catch ME
    id = ME.identifier;
end
p4 = pwd;
rows = sa_check(rows, T, 'a corrector error inside the run reaches the caller and the caller''s folder is current again', ...
    strcmp(id, 'SA:TEST:KILL') && strcmp(p4, here0), sprintf('%s | pwd after %s', id, p4));

% (5) the package found only as the current folder
pkg = fileparts(which('seeded_atlas_wrapper'));
p0 = path;
back_path = onCleanup(@() path(p0));
rmpath(pkg);
cd(pkg);
c5 = base; c5.run_dir = fullfile(d0, 'pkgcwd'); c5.corrector = sa_m3_stub('corrector', sa_m3_stub('state'));
id5 = '';
out5 = '';
try
    out5 = evalc('W5 = seeded_atlas_wrapper(c5);');
catch ME
    id5 = ME.identifier;
end
p5 = pwd;
path(p0);
clear back_path
ok5 = isempty(id5) && strcmp(W5.result.stop_reason, 'completed') && strcmp(p5, pkg) && ...
    contains(out5, 'not on the path');
if ok5
    D5 = sa_m3_equiv(c.run_dir, c5.run_dir);
    ok5 = D5.pass;
end
rows = sa_check(rows, T, ['launched from the package folder with the package not on the path: the current folder ' ...
    'is kept (one line says so), the run completes and its records equal the absolute run''s'], ok5, ...
    sprintf('%s | pwd after %s', id5, p5));
cd(old);
rmpath(E.stub);
end

function R = runtime_of(run_dir)
Z = load(fullfile(sa_layout('state', run_dir), 'settings_effective.mat'), 'runtime');
R = Z.runtime;
end

function s = first_line(run_dir)
t = strsplit(fileread(fullfile(run_dir, 'summary.txt')), newline);
s = t{1};
end
