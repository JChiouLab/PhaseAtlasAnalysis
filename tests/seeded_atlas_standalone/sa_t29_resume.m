function rows = sa_t29_resume(E)
%SA_T29_RESUME  Kill, then resume == uninterrupted (strict
%   resume equivalence) for kills mid-call (journal start, no
%   file), during an attempt save (planted .tmp_ file, no attempt file),
%   between the attempt file and its 'end' line, and at EVERY call boundary;
%   'interrupted' records count toward budget.max_calls only, never toward
%   edge, lineage, node or rescue caps.
T = 'T29 kill and resume';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't29');
fail = @(dd, b, ex) sa_m3_wall(dd, ex);
box = [2 1];
U = fullfile(d, 'U');
rU = sa_m3_run(U, box, 'fail', fail);
N = rU.calls;
% ---- (a) mid-call kills ---------------------------------------------------
ok = true; det = {};
for k = unique([2, 7, round(N / 2), N - 1])
    K = fullfile(d, sprintf('mid_%03d', k));
    kill_run(K, box, fail, k);
    r = resume(K, box, fail);
    D = sa_m3_equiv(U, K);
    ok = ok && D.pass && D.n_interrupted == 1 && isequal(r.replay_interrupted, k);
    det{end+1} = sprintf('k=%d %d', k, D.pass); %#ok<AGROW>
end
rows = sa_check(rows, T, sprintf('(a) kill mid-call at 4 points of %d calls: resumed == uninterrupted', N), ok, strjoin(det, ' '));
% ---- (b) kill during an attempt save --------------------------------------
k = 9; K = fullfile(d, 'save');
kill_run(K, box, fail, k);
fid = fopen(fullfile(sa_layout('state', K), 'attempts', sprintf('A%06d.mat.tmp_99999.mat', k)), 'w'); fwrite(fid, 'partial'); fclose(fid);
r = resume(K, box, fail);
D = sa_m3_equiv(U, K);
Z = load(fullfile(sa_layout('state', K), 'resume_1.mat'));
rows = sa_check(rows, T, '(b) kill during a save (planted .tmp_, no file): listed, never read; resumed == uninterrupted', ...
    D.pass && numel(Z.stray_tmp) == 1 && isequal(r.replay_interrupted, k), ...
    sprintf('stray %s; %s', strjoin(Z.stray_tmp, ','), strjoin(D.diffs, ' | ')));
% ---- (c) kill between the attempt file and its end line --------------------
k = 12; K = fullfile(d, 'fileend');
kill_run(K, box, fail, k + 1);
J = strsplit(strtrim(fileread(fullfile(sa_layout('state', K), 'journal.txt'))), newline);
J = J(~startsWith(J, sprintf('start A%06d', k + 1)) & ~startsWith(J, sprintf('end A%06d', k)));
write_lines(fullfile(sa_layout('state', K), 'journal.txt'), J);
r = resume(K, box, fail);
D = sa_m3_equiv(U, K);
J2 = fileread(fullfile(sa_layout('state', K), 'journal.txt'));
rows = sa_check(rows, T, '(c) attempt file present without its end line: complete, end appended, no interrupted record', ...
    D.pass && D.n_interrupted == 0 && isempty(r.replay_interrupted) && ...
    contains(J2, sprintf('end A%06d', k)), strjoin(D.diffs, ' | '));
% ---- (d) every call boundary -----------------------------------------------
bad = [];
for k = 1:N - 1
    K = fullfile(d, sprintf('bnd_%03d', k));
    kill_run(K, box, fail, k + 1);
    J = strsplit(strtrim(fileread(fullfile(sa_layout('state', K), 'journal.txt'))), newline);
    J = J(~startsWith(J, sprintf('start A%06d', k + 1)));
    write_lines(fullfile(sa_layout('state', K), 'journal.txt'), J);
    resume(K, box, fail);
    D = sa_m3_equiv(U, K);
    if ~(D.pass && D.n_interrupted == 0), bad(end+1) = k; end %#ok<AGROW>
    rmdir(K, 's');
end
rows = sa_check(rows, T, sprintf('(d) stop at every one of the %d call boundaries, resume: == uninterrupted', N - 1), ...
    isempty(bad), sprintf('failing boundaries %s', mat2str(bad)));
% ---- interrupted records and caps ---------------------------------------------
% a lineage cap binds; the kill hits a call of that lineage; the resumed
% lineage makes exactly as many calls as uninterrupted (the interrupted record
% is not counted toward the lineage cap)
o = struct('continuation', struct('max_calls_per_lineage', 6));
U2 = fullfile(d, 'U2');
r2 = sa_m3_mc_run(U2, o);
l2 = find(arrayfun(@(a) a.lineage == 2, r2.attempts));
kk = l2(3);
K2 = fullfile(d, 'K2');
try, sa_m3_mc_run(K2, o, 'kill_at', kk); catch, end
s = struct('resume', true); oo = o; oo.resume = true; %#ok<NASGU>
rK = sa_m3_mc_run(K2, oo);
D = sa_m3_equiv(U2, K2);
rows = sa_check(rows, T, ['interrupted record excluded from the lineage cap: lineage 2 makes the same ' ...
    '6 calls (cap) after a kill inside it'], D.pass && rK.Q.lin(2).calls == 6 && r2.Q.lin(2).calls == 6 && ...
    strcmp(rK.Q.lin(2).status, 'lineage_cap'), sprintf('calls %d/%d; %s', r2.Q.lin(2).calls, ...
    rK.Q.lin(2).calls, strjoin(D.diffs, ' | ')));
% ... but counted toward budget.max_calls: resume with max_calls = N
K3 = fullfile(d, 'K3'); kill_run(K3, box, fail, 5);
r3 = sa_m3_run(K3, box, 'fail', fail, 'settings', struct('resume', true, 'budget', struct('max_calls', N)));
rows = sa_check(rows, T, 'interrupted record counted toward budget.max_calls: max_calls = N stops one unit short', ...
    strcmp(r3.stop_reason, 'budget_calls') && r3.calls == N && sum(~strcmp({r3.attempts.outcome}, 'interrupted')) == N - 1, ...
    sprintf('%s after %d records', r3.stop_reason, r3.calls));
rmpath(E.stub);
end

function kill_run(K, box, fail, k)
try
    sa_m3_run(K, box, 'fail', fail, 'kill_at', k);
catch ME
    assert(strcmp(ME.identifier, 'SA:TEST:KILL'), ME.message);
end
end

function r = resume(K, box, fail)
r = sa_m3_run(K, box, 'fail', fail, 'settings', struct('resume', true));
end

function write_lines(f, L)
fid = fopen(f, 'w'); fprintf(fid, '%s\n', L{:}); fclose(fid);
end
