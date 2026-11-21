function rows = sa_t76_screen(E)
%SA_T76_SCREEN  bvp4c's mesh-limit warning stays
%   enabled and recorded but is not printed (report.solver_warnings brings it
%   back); other solve output is printed unchanged; progress stage and
%   status lines; neither changes what is solved or the order. Zero bvp4c
%   calls: stub corrector on the stub model; the warning is MATLAB's own
%   RelTolNotMet message, raised by warning(message(...)).
T = 'T76 screen output';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d0 = sa_m3_testdir(E, 't76');
wid = 'MATLAB:bvp4c:RelTolNotMet';

% ---- the strip rule on the two printed forms -------------------------------------------------------
msg = getString(message(wid, 10000, 10001, '0.0123', '1e-06'));
desk = sprintf('Warning: %s\n', msg);                                   % desktop form
batch = ['[Warning: Unable to meet the tolerance without using more than 10000 mesh points.' newline ...
    ' The last mesh of 9986 points and the solution are available in the output' newline ' argument.' newline ...
    ' The maximum residual is 0.0123, while requested accuracy is 1e-06.] ' newline];   % -batch form (wrapped)
ill = sprintf('Warning: %s\n', getString(message('MATLAB:bvp4c:IllCondJac', '1e-17')));
w0 = warning; warning('off', 'backtrace');
real = evalc('warning(message(wid, 10000, 10001, ''0.0123'', ''1e-06''));');   % this MATLAB's own printed form
warning(w0);
[r4, n4] = sa_progress('strip', ['x' newline real 'y' newline]);
[r1, n1] = sa_progress('strip', desk);
[r2, n2] = sa_progress('strip', batch);
[r3, n3] = sa_progress('strip', ['before' newline batch 'between' newline desk ill 'after' newline]);
rows = sa_check(rows, T, ['the strip rule removes the mesh-limit warning in the printed forms (plain, wrapped, ' ...
    'and the form this MATLAB prints, with its marker characters) and keeps every other line, including the ill-conditioned-Jacobian warning'], ...
    n1 == 1 && isempty(strtrim(r1)) && n2 == 1 && isempty(strtrim(r2)) && n3 == 2 && ...
    strcmp(r3, ['before' newline 'between' newline ill 'after' newline]) && n4 == 1 && ...
    strcmp(r4, ['x' newline 'y' newline]), sprintf('%d %d %d %d; this MATLAB prints %s', n1, n2, n3, n4, ...
    mat2str(double(real(1:3)))));

% ---- stub runs: hidden / shown / recorded -------------------------------------------------------------------
cnt = containers.Map({'n', 'seen', 'msg'}, {0, {}, {}});
mk = @(st) warn_corrector(sa_m3_stub('corrector', st), cnt);
dA = fullfile(d0, 'A'); dB = fullfile(d0, 'B'); dC = fullfile(d0, 'C');
stA = sa_m3_stub('state');
oA = evalc('rA = seeded_atlas_run(sa_m3_stub(''params''), struct(), sa_m3_stub(''settings'', dA, [1 1]), mk(stA));');
cnt('n') = 0; cnt('seen') = {}; cnt('msg') = {};
sB = sa_m3_stub('settings', dB, [1 1]); sB.report = struct('solver_warnings', true);
stB = sa_m3_stub('state');
oB = evalc('rB = seeded_atlas_run(sa_m3_stub(''params''), struct(), sB, mk(stB));');
seenB = cnt('seen'); msgB = cnt('msg');
nw = floor((rA.calls + 1) / 3);                           % the corrector warns on calls 2, 5, 8, ...
rows = sa_check(rows, T, ['default: no mesh-limit warning text on the screen, the one-line note once, the other ' ...
    'solve output on every call'], ~contains(oA, 'Unable to meet the tolerance') && ...
    numel(strfind(oA, 'mesh-limit warnings are not printed')) == 1 && ...
    numel(strfind(oA, 'stub solve output')) == rA.calls, ...
    sprintf('%d calls, %d warned; on screen: %d warning texts, %d notes, %d other lines', rA.calls, nw, ...
    numel(strfind(oA, 'Unable to meet the tolerance')), numel(strfind(oA, 'mesh-limit warnings are not printed')), ...
    numel(strfind(oA, 'stub solve output'))));
rows = sa_check(rows, T, 'report.solver_warnings = true prints every mesh-limit warning again (and no note)', ...
    numel(strfind(oB, 'Unable to meet the tolerance')) == nw && ~contains(oB, 'mesh-limit warnings are not printed') ...
    && numel(strfind(oB, 'stub solve output')) == rB.calls, sprintf('%d shown', numel(strfind(oB, 'Unable to meet'))));
aA = rows_wall0(rA.attempts); aB = rows_wall0(rB.attempts);
fa = arrayfun(@(a) load(fullfile(sa_layout('state', dA), a.file), 'warning_id', 'warning_msg'), rA.attempts);
fb = arrayfun(@(a) load(fullfile(sa_layout('state', dB), a.file), 'warning_id', 'warning_msg'), rB.attempts);
warned = strcmp({fa.warning_id}, wid);
seen_ok = all(cellfun(@(c) strcmp(c, wid), seenB(2:3:end))) && all(cellfun(@isempty, seenB(1:3:end))) && ...
    all(contains(msgB(2:3:end), 'Unable to meet the tolerance'));
rows = sa_check(rows, T, ['hidden or shown, the records are identical (attempt rows with wall_s cleared, every ' ...
    'attempt file''s warning_id / warning_msg); the warned calls record RelTolNotMet; lastwarn inside the ' ...
    'corrector holds the warning and its message, captured or not'], isequaln(aA, aB) && isequal(fa, fb) && ...
    sum(warned) == nw && isequal(find(warned), 2:3:rA.calls) && seen_ok, sprintf('%d warned records', sum(warned)));

% ---- progress lines -----------------------------------------------------------------------------------------
cnt('n') = 0; cnt('seen') = {}; cnt('msg') = {};
sC = sa_m3_stub('settings', dC, [1 1]); sC.report = struct('progress_interval_s', 1e-3);
stC = sa_m3_stub('state');
oC = evalc('rC = seeded_atlas_run(sa_m3_stub(''params''), struct(), sC, mk(stC));');
need = {'stage: inputs', 'stage: cold search', 'stage: continuation (interior phase, search pass 1)', ...
    'stage: rescue (interior phase, search pass 1)', 'stage: final reports (stop: completed)'};
st_lines = regexp(oC, '\[progress\] status:[^\n]*', 'match');
fields_ok = ~isempty(st_lines) && all(cellfun(@(l) contains(l, 'calls ') && contains(l, ' of 5000') && ...
    contains(l, 'roots ') && contains(l, 'elapsed '), st_lines));
rows = sa_check(rows, T, ['progress: a stage line at each stage start; with a 1 ms interval a status line after ' ...
    'calls, each with the stage, calls of the budget, roots and the elapsed time'], ...
    all(cellfun(@(s) contains(oC, s), need)) && numel(st_lines) >= 5 && fields_ok, ...
    sprintf('%d status lines; first: %s', numel(st_lines), first(st_lines)));
rows = sa_check(rows, T, 'progress and hidden warnings change nothing solved: rows of the 1 ms-interval run equal the default run''s', ...
    isequaln(rows_wall0(rC.attempts), aA) && isequal(rC.stop_reason, rA.stop_reason));
cold_line = regexp(oC, '\[progress\] stage: cold search[^\n]*', 'match', 'once');
rows = sa_check(rows, T, 'the cold stage line gives the progress inside the stage (targets called of the total)', ...
    contains(cold_line, '0 of 1 cold targets called'), cold_line);
ic = strfind(oC, 'runtime.clean:'); ip = strfind(oC, '[progress] stage: cold search');
rows = sa_check(rows, T, 'the runtime.clean line is printed before the search starts (after the inputs stage line)', ...
    ~isempty(ic) && ~isempty(ip) && ic(1) < ip(1));

% ---- resume-free: both leaves may change on resume -------------------------------------------------------------------
dD = fullfile(d0, 'D');
cnt('n') = 0; cnt('seen') = {}; cnt('msg') = {};
sD = sa_m3_stub('settings', dD, [1 1]); sD.budget.max_calls = 6;
sD.report = struct('solver_warnings', true, 'progress_interval_s', 1e-3);
stD = sa_m3_stub('state');
evalc('rD1 = seeded_atlas_run(sa_m3_stub(''params''), struct(), sD, mk(stD));');
sD.budget.max_calls = 5000; sD.resume = true; sD.report = struct('solver_warnings', false, 'progress_interval_s', Inf);
oD = evalc('rD2 = seeded_atlas_run(sa_m3_stub(''params''), struct(), sD, mk(stD));');
rows = sa_check(rows, T, ['report.solver_warnings and report.progress_interval_s are resume-free: a resume with both ' ...
    'changed runs to the same records; Inf prints stage lines only'], strcmp(rD1.stop_reason, 'budget_calls') && ...
    strcmp(rD2.stop_reason, 'completed') && isequaln(rows_wall0(rD2.attempts), aA) && ~contains(oD, 'status:') && ...
    contains(oD, 'replayed and checked first') && contains(oD, 'replaying stored calls'), rD2.stop_reason);
bad = '';
try
    sa_settings('resolve', setf(sa_m3_stub('settings', 'x', [1 1]), 'report', struct('progress_interval_s', 0)), ...
        sa_model('stub_two_pool'));
catch ME
    bad = ME.identifier;
end
rows = sa_check(rows, T, 'report.progress_interval_s must be positive (0 refused)', strcmp(bad, 'SA:SETTINGS:E_FIELD'), bad);
end

% =========================================================================

function c = warn_corrector(c0, cnt)
c = @(th, g, t, o) warn_call(c0, cnt, th, g, t, o);
end

function [sol, info] = warn_call(c0, cnt, th, g, t, o)
% Every third call raises MATLAB's own RelTolNotMet warning (enabled,
% backtrace off, as the real correctors set it) and records lastwarn, like
% the real adapters; every call also prints one other line.
n = cnt('n') + 1; cnt('n') = n;
prev = warning; cl = onCleanup(@() warning(prev));
warning('on', 'MATLAB:bvp4c:RelTolNotMet'); warning('off', 'backtrace');
lastwarn('', '');
if mod(n, 3) == 2
    warning(message('MATLAB:bvp4c:RelTolNotMet', 10000, 10001, '0.0123', '1e-06'));
end
fprintf('stub solve output %d\n', n);
[sol, info] = c0(th, g, t, o);
[m, i] = lastwarn;
s = cnt('seen'); s{end+1} = i; cnt('seen') = s;
q = cnt('msg'); q{end+1} = m; cnt('msg') = q;
info.warning_id = i;
if ~isempty(sol) && ~isempty(i), sol.stats.maxres = 1; end   % a warned call is refused, as bvp4c's would be
end

function a = rows_wall0(a)
for k = 1:numel(a), a(k).wall_s = 0; end
end

function s = first(c)
if isempty(c), s = ''; else, s = c{1}; end
end

function s = setf(s, f, v)
s.(f) = v;
end
