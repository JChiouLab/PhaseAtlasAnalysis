function rows = sa_t81_shard_repair(E)
%SA_T81_SHARD_REPAIR  A stop between an accepted call's attempt file and its
%   node's shard write, at a node that gets no later discovery: the replay
%   of the next resume rewrites that shard (sa_search repair_shards).
%   Zero bvp4c calls (stub model, the multi-seed run of sa_t28_lineages).
%   The stop is built exactly: the run stops (budget) after call k - 1; the
%   folder is copied; a resume makes call k live; then the folder is put
%   back to the copy except for call k's attempt file and the journal lines
%   up to its 'start' (what a kill right after the attempt write leaves:
%   no shard flush, no 'end' line, checkpoints as before). Two cases: the
%   node had a shard file before call k (it stays stamped below k; the multi-seed
%   run, two roots per node), and call k was the node's only discovery (no
%   shard file; the one-family run of sa_t29_resume, box [2 1] with the sa_m3_wall
%   refusals).
%   Checks: after the resume the shard is stamped k and byte-equal (after
%   the 128-byte MAT header, which holds the creation time) to the
%   uninterrupted run's; sa_m3_equiv passes; exactly one 'shard ...
%   rewritten' journal note; seeded_atlas_roots lists the roots of
%   roots_index.csv; a second resume and the resume of the clean
%   uninterrupted folder write no such note and change no shard byte.
%   Left as integrity stops (never rewritten): a shard stamped above its
%   node's last discovery, and a shard file of a node without a discovery.
T = 'T81 stale-shard repair at the end of the replay';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't81');
wall = @(dd, b, ex) sa_m3_wall(dd, ex);
one = @(dd, x) sa_m3_run(dd, [2 1], 'fail', wall, 'settings', x);
cases = struct('name', {'stamped_below', 'missing'}, 'run', {@(dd, x) sa_m3_mc_run(dd, x), one}, ...
    'earlier', {true, false});
for c = 1:numel(cases)
    run = cases(c).run;
    U = fullfile(d, ['U_' cases(c).name]);
    rU = run(U, struct());
    A = rU.attempts;
    acc = find(strcmp({A.outcome}, 'accepted'));
    later = arrayfun(@(k) any(strcmp({A(acc(acc > k)).target_key}, A(k).target_key)), acc);
    earlier = arrayfun(@(k) any(strcmp({A(acc(acc < k)).target_key}, A(k).target_key)), acc);
    k = acc(find(~later & earlier == cases(c).earlier & acc > 5, 1, 'first'));
    cases(c).k = k; cases(c).U = U; cases(c).rU = rU;
    K = fullfile(d, cases(c).name);
    if isempty(k)
        rows = sa_check(rows, T, sprintf('%s: a suitable call exists in the stub run', cases(c).name), false, '');
        continue
    end
    key = A(k).target_key;
    sf = fullfile(sa_layout('state', K), 'roots', sa_store_io('shard_file', key));
    run(K, struct('budget', struct('max_calls', k - 1)));
    had = exist(sf, 'file') == 2;
    stamp0 = -1;
    if had, Z0 = load(sf, 'attempt_count'); stamp0 = Z0.attempt_count; end
    snap = [K '_snap'];
    copyfile(K, snap);
    run(K, struct('resume', true, 'budget', struct('max_calls', k)));
    af = fullfile(sa_layout('state', K), 'attempts', sprintf('A%06d.mat', k));
    J = strsplit(strtrim(fileread(fullfile(sa_layout('state', K), 'journal.txt'))), newline);
    js = find(startsWith(J, sprintf('start A%06d', k)), 1);
    tmp = [K '_A.mat'];
    copyfile(af, tmp);
    rmdir(K, 's');
    movefile(snap, K);
    movefile(tmp, af);
    write_lines(fullfile(sa_layout('state', K), 'journal.txt'), J(1:js));
    stale = (~had && exist(sf, 'file') ~= 2) || (had && stamp_of(sf) == stamp0);
    r = run(K, struct('resume', true));
    D = sa_m3_equiv(U, K);
    uf = fullfile(sa_layout('state', U), 'roots', sa_store_io('shard_file', key));
    beq = bytes_equal(uf, sf);
    notes = repair_notes(K);
    Rk = seeded_atlas_roots(K);
    [ri_ids, ri_n] = index_ids(K);
    roots_ok = numel(Rk) == ri_n && isequal(sort({Rk.root_id}), sort(ri_ids)) && numel(Rk) == rU.n_roots;
    rows = sa_check(rows, T, sprintf(['%s: stop after A%06d''s attempt file (node %s, shard before: %s) -> after ' ...
        'the resume the shard is stamped %d, byte-equal to U''s after the header, sa_m3_equiv passes, one ' ...
        'journal note'], cases(c).name, k, key, tern(had, sprintf('stamped %d', stamp0), 'none'), k), ...
        stale && strcmp(r.stop_reason, 'completed') && stamp_of(sf) == k && stamp_of(uf) == k && beq && D.pass && ...
        numel(notes) == 1 && contains(notes{1}, key) && contains(notes{1}, sprintf('stamped %d', k)), ...
        sprintf('%s | notes: %s', strjoin(D.diffs, ' | '), strjoin(notes, ' || ')));
    rows = sa_check(rows, T, sprintf(['%s: seeded_atlas_roots lists the %d roots of roots_index.csv (as the ' ...
        'uninterrupted run)'], cases(c).name, ri_n), roots_ok, sprintf('%d post-tool roots, %d index rows', ...
        numel(Rk), ri_n));
    % a second resume: nothing left to repair
    b0 = shard_bytes(K);
    run(K, struct('resume', true));
    rows = sa_check(rows, T, sprintf('%s: a second resume writes no repair note and leaves every shard byte', ...
        cases(c).name), numel(repair_notes(K)) == 1 && isequal(b0, shard_bytes(K)), '');
end
% ---- the clean folder (the multi-seed run) --------------------------------------------------------
U = cases(1).U; rU = cases(1).rU; A = rU.attempts;
C = fullfile(d, 'clean');
copyfile(U, C);
b0 = shard_bytes(C);
rC = sa_m3_mc_run(C, struct('resume', true));
DC = sa_m3_equiv(U, C);
rows = sa_check(rows, T, 'resume of the clean uninterrupted folder: zero repair notes, every shard byte unchanged', ...
    isempty(repair_notes(C)) && isequal(b0, shard_bytes(C)) && DC.pass && strcmp(rC.stop_reason, 'completed') && ...
    contains(fileread(fullfile(sa_layout('state', C), 'journal.txt')), 'replay complete'), strjoin(DC.diffs, ' | '));
% ---- stays an integrity stop: a shard stamped above its node's last discovery -----------
k = cases(1).k; key = A(k).target_key;
H = fullfile(d, 'above');
copyfile(U, H);
hf = fullfile(sa_layout('state', H), 'roots', sa_store_io('shard_file', key));
Z = sa_store_io('read_shard', hf);
sa_store_io('save_shard', Z.store, hf, key, k + 1);
[eid, msg] = resume_error(H);
rows = sa_check(rows, T, sprintf(['a shard stamped above its node''s last discovery (%s stamped %d, last A%06d): ' ...
    'integrity stop "never met", not rewritten'], key, k + 1, k), strcmp(eid, 'SA:INTEGRITY') && ...
    contains(msg, 'never met') && stamp_of(hf) == k + 1 && isempty(repair_notes(H)), msg);
% ---- stays an integrity stop: a shard file of a node without a discovery ------------------
N = fullfile(d, 'nodisc');
copyfile(U, N);
nk = '';
for q = 1:numel(rU.L.key)
    if ~isKey(rU.shards, rU.L.key{q}), nk = rU.L.key{q}; break; end
end
nf = fullfile(sa_layout('state', N), 'roots', sa_store_io('shard_file', nk));
sa_store_io('save_shard', sa_root_store('empty', Z.store.settings), nf, nk, k);
[eid, msg] = resume_error(N);
rows = sa_check(rows, T, sprintf(['a shard file of a node without a discovery (%s): integrity stop "never met", ' ...
    'not rewritten'], nk), ~isempty(nk) && strcmp(eid, 'SA:INTEGRITY') && contains(msg, 'never met') && ...
    isempty(repair_notes(N)), msg);
rmpath(E.stub);
end

function [eid, msg] = resume_error(K)
eid = ''; msg = '';
try
    sa_m3_mc_run(K, struct('resume', true));
catch ME
    eid = ME.identifier; msg = ME.message;
end
end

function s = stamp_of(f)
s = -1;
if exist(f, 'file') == 2, Z = load(f, 'attempt_count'); s = Z.attempt_count; end
end

function L = repair_notes(K)
J = strsplit(strtrim(fileread(fullfile(sa_layout('state', K), 'journal.txt'))), newline);
L = J(startsWith(J, 'note shard ') & contains(J, ' rewritten stamped '));
end

function [ids, n] = index_ids(K)
L = strsplit(strtrim(fileread(fullfile(K, 'roots_index.csv'))), newline);
h = strsplit(L{1}, ',');
c = find(strcmp(h, 'root_id'), 1);
ids = cell(1, numel(L) - 1);
for q = 2:numel(L)
    p = regexp(L{q}, ',(?=(?:[^"]*"[^"]*")*[^"]*$)', 'split');
    ids{q - 1} = strrep(p{c}, '"', '');
end
n = numel(ids);
end

function B = shard_bytes(K)
f = dir(fullfile(sa_layout('state', K), 'roots', '*.mat'));
B = struct();
for q = 1:numel(f)
    b = read_bytes(fullfile(f(q).folder, f(q).name));
    B.(matlab.lang.makeValidName(f(q).name)) = b;
end
end

function tf = bytes_equal(a, b)
x = read_bytes(a); y = read_bytes(b);
tf = numel(x) == numel(y) && numel(x) > 128 && isequal(x(129:end), y(129:end));
end

function b = read_bytes(f)
b = [];
fid = fopen(f, 'r');
if fid < 0, return; end
b = fread(fid, Inf, '*uint8');
fclose(fid);
end

function write_lines(f, L)
fid = fopen(f, 'w'); fprintf(fid, '%s\n', L{:}); fclose(fid);
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
