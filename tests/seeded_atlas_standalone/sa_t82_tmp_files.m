function rows = sa_t82_tmp_files(E)
%SA_T82_TMP_FILES  sa_write's temporary files in the shard listings, and
%   the repair's intent note. Zero bvp4c calls (stub model, the
%   multi-seed run of sa_t28_lineages).
%   A kill inside sa_write leaves a temporary sibling beside the file it was
%   replacing. On MATLAB R2025b save keeps the name as given, so the file is
%   <name>.mat.tmp_<pid>, which *.mat does not match (case sa_write below).
%   The hand-planted cases use <name>.mat.tmp_<pid>.mat, the name a save
%   that appends .mat would give (sa_write's fallback branch), which *.mat
%   does match: they are a GUARD for the listings' '.tmp_' filter, not a
%   name this MATLAB produces. The stale folder is built as in sa_t81_shard_repair (a stop
%   between an accepted call k's attempt file and its shard write, at a node
%   with an earlier discovery and no later one), then:
%     complete   (guard) a complete planted temporary sibling of the stale
%                shard (the uninterrupted run's shard bytes: a kill between
%                the temporary save and the rename)
%     truncated  (guard) half of those bytes (a kill during the temporary save)
%     missing    (guard) the shard file removed, a complete planted sibling left
%     sa_write   the temporary file made by sa_write itself: the stale shard
%                is locked (chflags uchg) so that sa_write's save succeeds and
%                its rename fails, leaving exactly the file a kill between the
%                save and the rename leaves; its name must be
%                <name>.mat.tmp_<pid> with nothing appended (for the resume
%                it is given another pid, as a killed process's file has:
%                this process's own next save would reuse the name)
%   Each resume completes; the shard is stamped k and byte-equal (after the
%   128-byte MAT header) to the uninterrupted run's; no integrity note; one
%   'repair starts' and one 'rewritten' note; the resume note counts one
%   stray temporary file, which is left unread and byte-unchanged;
%   sa_m3_equiv passes once it is set aside.
%   A stop between the intent note and the save (roots/ made read-only, so
%   the save fails right after the note): the folder is as before plus the
%   note; the next resume repairs, with a new 'repair starts' note followed
%   by its 'rewritten' note. A stop after the save, before its 'rewritten'
%   note (the saved file put in place by hand on the same folder): the next
%   resume finds nothing to repair and writes no shard note.
%   A completed folder with planted temporary files (an older version of one
%   shard, a truncated copy of another, an attempt file's copy):
%   seeded_atlas_roots (verified and not), seeded_atlas_plot and
%   seeded_atlas_compact give what they give on the clean folder and leave
%   the planted files in place.
T = 'T82 temporary files in the shard listings; repair intent note';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't82');
run = @(dd, x) sa_m3_mc_run(dd, x);
U = fullfile(d, 'U');
rU = run(U, struct());
A = rU.attempts;
acc = find(strcmp({A.outcome}, 'accepted'));
later = arrayfun(@(k) any(strcmp({A(acc(acc > k)).target_key}, A(k).target_key)), acc);
earlier = arrayfun(@(k) any(strcmp({A(acc(acc < k)).target_key}, A(k).target_key)), acc);
k = acc(find(~later & earlier & acc > 5, 1, 'first'));
if isempty(k)
    rows = sa_check(rows, T, 'a suitable call exists in the stub run', false, '');
    rmpath(E.stub);
    return
end
key = A(k).target_key;
sfn = @(D) fullfile(sa_layout('state', D), 'roots', sa_store_io('shard_file', key));
tmpn = @(D) [sfn(D) '.tmp_99999.mat'];
uf = sfn(U);
ub = read_bytes(uf);
K = fullfile(d, 'K');
make_stale(run, K, k);
s0 = stamp_of(sfn(K));
k0 = read_bytes(sfn(K));

% ---- resume with a temporary sibling in roots/ -------------------------------------------------
cases = {'complete', 'truncated', 'missing'};
for c = 1:numel(cases)
    Kc = fullfile(d, ['P_' cases{c}]);
    copyfile(K, Kc);
    switch cases{c}
        case 'complete',  write_bytes(tmpn(Kc), ub); was = sprintf('was stamped %d', s0);
        case 'truncated', write_bytes(tmpn(Kc), ub(1:floor(numel(ub) / 2))); was = sprintf('was stamped %d', s0);
        case 'missing',   delete(sfn(Kc)); write_bytes(tmpn(Kc), ub); was = 'was missing';
    end
    tb = read_bytes(tmpn(Kc));
    [eid, msg, r] = resume(Kc);
    J = journal(Kc);
    st = shard_notes(J, 'repair starts'); rw = shard_notes(J, ' rewritten stamped ');
    rn = J(startsWith(J, 'note resume '));
    ok = isempty(eid) && strcmp(r.stop_reason, 'completed') && stamp_of(sfn(Kc)) == k && ...
        bytes_equal(uf, sfn(Kc)) && ~any(startsWith(J, 'note integrity')) && numel(st) == 1 && ...
        numel(rw) == 1 && contains(st{1}, key) && contains(rw{1}, key) && contains(rw{1}, was) && ...
        ~isempty(rn) && contains(rn{end}, ', 1 stray temporary file(s)') && isequal(read_bytes(tmpn(Kc)), tb);
    [Rk, e1] = attempt(@() seeded_atlas_roots(Kc));
    if exist(tmpn(Kc), 'file') == 2, delete(tmpn(Kc)); end   % set aside: sa_m3_equiv compares the shard file sets
    [D, e2] = attempt(@() sa_m3_equiv(U, Kc));
    if isempty(D), D = struct('pass', false, 'diffs', {{e2}}); end
    rows = sa_check(rows, T, sprintf(['guard: %s planted temporary sibling of the stale shard %s (stamped %d, last discovery ' ...
        'A%06d): the resume completes, the shard is stamped %d and byte-equal to the uninterrupted run''s, no ' ...
        'integrity note, one repair-starts and one rewritten note, one stray temporary file counted and left ' ...
        'unread; sa_m3_equiv passes; seeded_atlas_roots lists %d roots'], cases{c}, key, s0, k, k, rU.n_roots), ...
        ok && D.pass && isempty(e1) && numel(Rk) == rU.n_roots, sprintf('%s %s | %s | %s | %s', eid, msg, e1, ...
        strjoin(D.diffs, ' | '), strjoin([st, rw], ' || ')));
end

% ---- the temporary file made by sa_write itself ----------------------------------------------------
% chflags is macOS/BSD only; elsewhere this case is skipped with a printed reason
if ~ismac
    fprintf('T82: the sa_write case is skipped on this platform (it locks a file with chflags, macOS/BSD only)\n');
else
Kw = fullfile(d, 'P_sa_write');
copyfile(K, Kw);
fw = sfn(Kw);
[s1, o1] = system(sprintf('chflags uchg "%s"', fw));          % the rename onto the shard fails
ew = '';
try
    sa_write('mat', fw, load(uf));                              % the uninterrupted run's shard variables
catch ME
    ew = ME.identifier;
end
[s2, o2] = system(sprintf('chflags nouchg "%s"', fw));
tw = dir([fw '.tmp_*']);
name_ok = isscalar(tw) && ~isempty(regexp(tw(1).name, ['^' regexptranslate('escape', ...
    sa_store_io('shard_file', key)) '\.tmp_\d+$'], 'once'));
% the killed process is gone and the resume runs under another pid: the file keeps sa_write's form with
% another process's id (this process's own next save of the shard would reuse its name)
twf = ''; tb = [];
if isscalar(tw)
    twf = [fw '.tmp_' sprintf('%d', feature('getpid') + 1)];
    movefile(fullfile(tw(1).folder, tw(1).name), twf);
    tb = read_bytes(twf);
end
[eid, msg, r] = resume(Kw);
J = journal(Kw);
st = shard_notes(J, 'repair starts'); rw = shard_notes(J, ' rewritten stamped ');
rn = J(startsWith(J, 'note resume '));
part = [s1 == 0, s2 == 0, strcmp(ew, 'SA:WRITE'), name_ok, isempty(eid), strcmp(r.stop_reason, 'completed'), ...
    stamp_of(sfn(Kw)) == k, bytes_equal(uf, sfn(Kw)), ~any(startsWith(J, 'note integrity')), numel(st) == 1, ...
    numel(rw) == 1, ~isempty(rn) && contains(rn{end}, ', 1 stray temporary file(s)'), ~isempty(twf) && isequal(read_bytes(twf), tb)];
ok = all(part);
[Rk, e1] = attempt(@() seeded_atlas_roots(Kw));
[D, e2] = attempt(@() sa_m3_equiv(U, Kw));                     % *.mat does not match the file: no set-aside
if isempty(D), D = struct('pass', false, 'diffs', {{e2}}); end
rows = sa_check(rows, T, sprintf(['the temporary file sa_write itself leaves (save done, rename refused) is ' ...
    'named <shard>.mat.tmp_<pid> with nothing appended; the resume completes, repairs the shard (stamped %d, ' ...
    'byte-equal to the uninterrupted run''s), counts one stray temporary file and leaves it unread; ' ...
    'sa_m3_equiv passes with the file in place; seeded_atlas_roots lists %d roots'], k, rU.n_roots), ...
    ok && D.pass && isempty(e1) && numel(Rk) == rU.n_roots, sprintf('%s | %s | %s %s | %s %s | %s | %s | %s', ...
    mat2str(part), strjoin({tw.name}, ','), ew, strtrim([o1 o2]), eid, msg, e1, strjoin(D.diffs, ' | '), strjoin([st, rw], ' || ')));
end

% ---- a stop between the intent note and the save -----------------------------------------------
Ki = fullfile(d, 'intent');
copyfile(K, Ki);
rd = fullfile(sa_layout('state', Ki), 'roots');
fileattrib(rd, '-w', 'a');
[eid1, msg1] = resume(Ki);
fileattrib(rd, '+w', 'u');
J1 = journal(Ki);
left = dir(fullfile(rd, '*.tmp_*'));
stop_ok = ~isempty(eid1) && startsWith(J1{end}, sprintf('note shard %s repair starts', key)) && ...
    isempty(shard_notes(J1, ' rewritten stamped ')) && isequal(read_bytes(sfn(Ki)), k0) && isempty(left);
rows = sa_check(rows, T, sprintf(['a save that fails right after the intent note (roots/ read-only) leaves the ' ...
    'folder as before plus the note "shard %s repair starts", no rewritten note, no temporary file'], key), ...
    stop_ok, sprintf('%s %s | last line: %s', eid1, msg1, J1{end}));
% the save had landed before the stop: the same folder with the file in place
Ks = [Ki '_saved'];
copyfile(Ki, Ks);
write_bytes(sfn(Ks), ub);
[eid3, msg3, r3] = resume(Ks);
J3 = journal(Ks);
D3 = sa_m3_equiv(U, Ks);
rows = sa_check(rows, T, ['a stop after the save, before its rewritten note: the next resume finds nothing to ' ...
    'repair (no new shard note) and completes; sa_m3_equiv passes'], isempty(eid3) && ...
    strcmp(r3.stop_reason, 'completed') && numel(shard_notes(J3, 'repair starts')) == 1 && ...
    isempty(shard_notes(J3, ' rewritten stamped ')) && D3.pass && ~any(startsWith(J3, 'note integrity')), ...
    sprintf('%s %s | %s', eid3, msg3, strjoin(D3.diffs, ' | ')));
% the save had not landed: the next resume repairs
[eid2, msg2, r2] = resume(Ki);
J2 = journal(Ki);
st = find(startsWith(J2, sprintf('note shard %s repair starts', key)));
rw = find(startsWith(J2, sprintf('note shard %s rewritten stamped %d', key, k)));
D2 = sa_m3_equiv(U, Ki);
rows = sa_check(rows, T, ['the next resume repairs: two repair-starts notes, one rewritten note right after the ' ...
    'second; the shard is stamped k and byte-equal to the uninterrupted run''s; sa_m3_equiv passes'], ...
    isempty(eid2) && strcmp(r2.stop_reason, 'completed') && numel(st) == 2 && isscalar(rw) && rw == st(2) + 1 && ...
    stamp_of(sfn(Ki)) == k && bytes_equal(uf, sfn(Ki)) && D2.pass && ~any(startsWith(J2, 'note integrity')), ...
    sprintf('%s %s | starts %s rewritten %s | %s', eid2, msg2, mat2str(st), mat2str(rw), strjoin(D2.diffs, ' | ')));

% ---- the post tools on a completed folder with planted temporary files --------------------------
ks = rU.shards.keys;
k2 = ks{find(~strcmp(ks, key), 1)};
sf2 = @(D) fullfile(sa_layout('state', D), 'roots', sa_store_io('shard_file', k2));
C = fullfile(d, 'clean_tmp');
copyfile(U, C);
Uc = fullfile(d, 'clean');
copyfile(U, Uc);
b2 = read_bytes(sf2(U));
planted = {tmpn(C), [sf2(C) '.tmp_99999.mat'], ...
    fullfile(sa_layout('state', C), 'attempts', 'A000001.mat.tmp_99999.mat')};
write_bytes(planted{1}, k0);                                   % an older version of the node's shard
write_bytes(planted{2}, b2(1:floor(numel(b2) / 2)));           % a truncated copy
copyfile(fullfile(sa_layout('state', C), 'attempts', sprintf('A%06d.mat', k)), planted{3});
pb = cellfun(@read_bytes, planted, 'UniformOutput', false);
Rv = seeded_atlas_roots(Uc); [Rc, e1] = attempt(@() seeded_atlas_roots(C));
Rv0 = seeded_atlas_roots(Uc, struct('verify', false)); [Rc0, e2] = attempt(@() seeded_atlas_roots(C, struct('verify', false)));
rows = sa_check(rows, T, ['guard: seeded_atlas_roots (verified and not) on a folder with planted temporary files ' ...
    '(an older version of one shard, a truncated copy of another, a copy of an attempt file) equals the clean ' ...
    'folder''s'], isempty(e1) && isempty(e2) && isequaln(Rv, Rc) && isequaln(Rv0, Rc0) && numel(Rc) == rU.n_roots, ...
    sprintf('%d / %d roots | %s | %s', numel(Rc), numel(Rv), e1, e2));
Pu = seeded_atlas_plot(Uc); [Pc, e1] = attempt(@() seeded_atlas_plot(C));
ok = false;
if isempty(e1)
    Ku = load(fullfile(Uc, 'plots', 'roots_3d_key.mat')); Kc = load(fullfile(C, 'plots', 'roots_3d_key.mat'));
    ok = Pc.n_roots == Pu.n_roots && isequal(Ku.key.roots, Kc.key.roots) && ...
        strcmp(fileread(fullfile(Uc, 'plots', 'roots_3d.csv')), fileread(fullfile(C, 'plots', 'roots_3d.csv')));
end
rows = sa_check(rows, T, ['guard: seeded_atlas_plot on that folder: the same roots, plot key (shard names and SHA-256) ' ...
    'and roots_3d.csv as the clean folder''s'], ok, sprintf('%d roots | %s', Pu.n_roots, e1));
Ou = seeded_atlas_compact(Uc); [Oc, e1] = attempt(@() seeded_atlas_compact(C));
ok = false; det = e1;
if isempty(e1)
    Mu = load(fullfile(sa_layout('state', Uc), 'compaction.mat'));
    Mc = load(fullfile(sa_layout('state', C), 'compaction.mat'));
    same_rm = isequal({Mu.removed.file}, {Mc.removed.file}) && isequal({Mu.removed.sha256}, {Mc.removed.sha256});
    kept = all(cellfun(@(f, b) isequal(read_bytes(f), b), planted, pb));
    ok = Oc.n_roots == Ou.n_roots && Oc.n_removed == Ou.n_removed && same_rm && kept;
    det = sprintf('removed %d / %d; roots %d / %d; planted kept %d', Oc.n_removed, Ou.n_removed, Oc.n_roots, ...
        Ou.n_roots, kept);
end
rows = sa_check(rows, T, ['guard: seeded_atlas_compact on that folder removes the same files as on the clean folder, ' ...
    'verifies the same roots and leaves the planted files byte-unchanged'], ok, det);
rmpath(E.stub);
end

% =========================================================================

function make_stale(run, K, k)
% as in sa_t81_shard_repair: stop after call k - 1, copy, make call k live, then put the copy
% back except call k's attempt file and the journal up to its 'start'
run(K, struct('budget', struct('max_calls', k - 1)));
snap = [K '_snap'];
copyfile(K, snap);
run(K, struct('resume', true, 'budget', struct('max_calls', k)));
af = fullfile(sa_layout('state', K), 'attempts', sprintf('A%06d.mat', k));
J = journal(K);
js = find(startsWith(J, sprintf('start A%06d', k)), 1);
tmp = [K '_A.mat'];
copyfile(af, tmp);
rmdir(K, 's');
movefile(snap, K);
movefile(tmp, af);
fid = fopen(fullfile(sa_layout('state', K), 'journal.txt'), 'w'); fprintf(fid, '%s\n', J{1:js}); fclose(fid);
end

function [v, err] = attempt(f)
v = []; err = '';
try
    v = f();
catch ME
    err = sprintf('%s: %s', ME.identifier, ME.message);
end
end

function [eid, msg, r] = resume(K)
eid = ''; msg = ''; r = struct('stop_reason', '');
try
    r = sa_m3_mc_run(K, struct('resume', true));
catch ME
    eid = ME.identifier; msg = ME.message;
end
end

function J = journal(K)
J = strsplit(strtrim(fileread(fullfile(sa_layout('state', K), 'journal.txt'))), newline);
end

function L = shard_notes(J, what)
L = J(startsWith(J, 'note shard ') & contains(J, what));
end

function s = stamp_of(f)
s = -1;
if exist(f, 'file') == 2, Z = load(f, 'attempt_count'); s = Z.attempt_count; end
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

function write_bytes(f, b)
fid = fopen(f, 'w'); fwrite(fid, b, 'uint8'); fclose(fid);
end
