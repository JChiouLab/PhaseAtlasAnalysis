function rows = sa_t43_am9(E)
%SA_T43_AM9  Compact shard storage: shard
%   serialization omits a discovery's returned_profile / comparison_profile
%   that is bit-identical to its root's representative_profile, with a flag
%   (format 'sa_store_io/shard-2'); every reader restores it.
%   (a) a stub run ('comparison' mesh) writes shard-2 files; the raw form
%       holds [] exactly where flagged; read_shard gives the final in-memory
%       store (isequaln); load_run_shards (with its re-derivation) equals it;
%   (b) the same shards rewritten in the older shard-1 form load
%       unchanged (read_shard and load_run_shards);
%   (c) a flag that disagrees with the stored field is refused (E_DEDUP);
%   (d) bit-strict: a profile isequal to the representative but differing
%       in the sign of a zero is NOT omitted, and the round trip is exact
%       in every IEEE bit.
%   The checkpoint comparison on the decoded form is that of sa_t29_resume,
%   sa_t38_39_budget_checkpoint and sa_t42_am8 (sa_m3_equiv and the resume compare read through read_shard).
T = 'T43 AM-9 shard dedup';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
M = sa_model('stub_two_pool');
d = sa_m3_testdir(E, 't43');
% ---- (a) --------------------------------------------------------------------------
run_dir = fullfile(d, 'run');
res = sa_m3_mc_run(run_dir, struct('identity', struct('candidate_mesh', 'comparison')));
Fm = sa_store_io('full_map', res);   % result.shards are light (memory); their full form
f = dir(fullfile(sa_layout('state', run_dir), 'roots', '*.mat'));
fmt = true; raw = true; same = true; nflag = 0; nfield = 0; keys = {};
for k = 1:numel(f)
    fk = fullfile(f(k).folder, f(k).name);
    Z = load(fk);
    fmt = fmt && strcmp(Z.format, 'sa_store_io/shard-2') && isfield(Z, 'profile_dedup');
    [ok, nf, nt] = raw_ok(Z);
    raw = raw && ok; nflag = nflag + nf; nfield = nfield + nt;
    S = sa_store_io('read_shard', fk);
    keys{end+1} = S.key; %#ok<AGROW>
    same = same && isKey(Fm, S.key) && isequaln(S.store, Fm(S.key));
end
rows = sa_check(rows, T, '(a) every shard is shard-2; raw [] exactly where flagged; read_shard == final in-memory store', ...
    numel(f) == res.shards.Count && fmt && raw && same && nflag > 0, ...
    sprintf('%d shards, %d of %d profile fields omitted', numel(f), nflag, nfield));
L = sa_store_io('load_run_shards', sa_layout('state', run_dir), M);
ok = L.Count == res.shards.Count;
for k = keys, ok = ok && isequaln(L(k{1}), Fm(k{1})); end
rows = sa_check(rows, T, '(a) load_run_shards (re-derivation + kernel replay) of shard-2 files == in-memory stores', ok);
% ---- (b) --------------------------------------------------------------------------
old_dir = fullfile(d, 'old');
copyfile(run_dir, old_dir);
for k = 1:numel(f)
    fk = fullfile(sa_layout('state', old_dir), 'roots', f(k).name);
    S = sa_store_io('read_shard', fk);
    sa_write('mat', fk, struct('store', S.store, 'format', 'sa_store_io/shard-1', ...
        'key', S.key, 'attempt_count', S.attempt_count));
end
L1 = sa_store_io('load_run_shards', sa_layout('state', old_dir), M);
ok = L1.Count == res.shards.Count;
for k = keys, ok = ok && isequaln(L1(k{1}), Fm(k{1})); end
S1 = sa_store_io('read_shard', fullfile(sa_layout('state', old_dir), 'roots', f(1).name));
Z1 = load(fullfile(sa_layout('state', old_dir), 'roots', f(1).name));
rows = sa_check(rows, T, '(b) older shard-1 files load unchanged (load_run_shards; read_shard returns the file as stored)', ...
    ok && isequaln(S1, Z1) && ~isfield(Z1, 'profile_dedup'));
% ---- (c) --------------------------------------------------------------------------
fk = fullfile(sa_layout('state', run_dir), 'roots', f(1).name);
fc = fullfile(d, 'tamper.mat');
Z = load(fk);
Zt = Z; Zt.profile_dedup{1}(2, 1) = false;           % flag cleared, field still []
sa_write('mat', fc, Zt);
g1 = errid(@() sa_store_io('read_shard', fc));
Zt = Z; q = find(~Zt.profile_dedup{1}(1, :), 1);
note = 'no unflagged returned profile in root 1';
if isempty(q)                                        % flag a stored field instead: refuse too
    Zt.store.roots(1).discoveries(1).comparison_profile = Z.store.roots(1).representative_profile;
    note = 'flagged comparison profile stored';
else
    Zt.profile_dedup{1}(1, q) = true;
end
sa_write('mat', fc, Zt);
g2 = errid(@() sa_store_io('read_shard', fc));
Zt = Z; Zt = rmfield(Zt, 'profile_dedup');
sa_write('mat', fc, Zt);
g3 = errid(@() sa_store_io('read_shard', fc));
rows = sa_check(rows, T, '(c) flag/field disagreement and a shard-2 file without flags are refused', ...
    strcmp(g1, 'SA:STORE_IO:E_DEDUP') && strcmp(g2, 'SA:STORE_IO:E_DEDUP') && strcmp(g3, 'SA:STORE_IO:E_FORMAT'), ...
    sprintf('%s / %s (%s) / %s', g1, g2, note, g3));
% ---- (d) --------------------------------------------------------------------------
p = sa_m3_stub('params'); th = M.normalize(p);
S = sa_settings('resolve', struct('model', 'stub_two_pool', 'run_dir', 'unused'), M);
K = sa_identity('store_settings', M, th, S, M.field_scale(th, [], S));
cm = K.comparison_mesh;
p0 = repmat(2 + 0.1 * cos(pi / th.Lh * cm), 4, 1);
p0(1, 5) = 0; pm = p0; pm(1, 5) = -0;
mk = @(prof, id) struct('target_key', 'k:0,0', 'coordinate', [0 0], 'mesh', cm, 'profile', prof, ...
    'seed_family', 'X_X', 'observed_family', 'X_X', 'route', 'cold', 'attempt_id', id, 'valid', true, ...
    'validity', struct('valid', true));
st = sa_root_store('empty', K);
st = sa_root_store('add', st, mk(p0, 'A000001'));
[st, o2] = sa_root_store('add', st, mk(pm, 'A000002'));
fz = fullfile(d, 'signed_zero.mat');
sa_store_io('save_shard', st, fz, 'k:0,0', 2);
Z = load(fz);
R = sa_store_io('read_shard', fz);
L = sa_store_io('load_shard', fz, M, @(dd) deal(true, {}));
d2 = st.roots(1).discoveries(2);
rows = sa_check(rows, T, ['(d) bit-strict: a returned profile isequal to the representative but with -0 for +0 is ' ...
    'stored, not omitted; read and load are exact in every bit'], ...
    strcmp(o2.status, 'matched_existing') && isequal(d2.returned_profile, st.roots(1).representative_profile) && ...
    isequal(Z.profile_dedup{1}(:, 1), [true; true]) && ~Z.profile_dedup{1}(1, 2) && ...
    Z.profile_dedup{1}(2, 2) == bits(d2.comparison_profile, st.roots(1).representative_profile) && ...
    isequaln(R.store, st) && isequaln(L, st) && all_bits(R.store, st) && all_bits(L, st), ...
    sprintf('flags %s; returned(1,5) sign bit after read %d', mat2str(Z.profile_dedup{1}), ...
    signbit_of(R.store.roots(1).discoveries(2).returned_profile(1, 5))));
rmpath(E.stub);
end

function [ok, nf, nt] = raw_ok(Z)
names = {'returned_profile', 'comparison_profile'};
ok = iscell(Z.profile_dedup) && numel(Z.profile_dedup) == numel(Z.store.roots);
nf = 0; nt = 0;
for r = 1:numel(Z.store.roots)
    F = Z.profile_dedup{r};
    ok = ok && F(2, 1);                              % the first discovery's comparison IS the representative
    for q = 1:numel(Z.store.roots(r).discoveries)
        for k = 1:2
            ok = ok && (F(k, q) == isempty(Z.store.roots(r).discoveries(q).(names{k})));
            nf = nf + F(k, q); nt = nt + 1;
        end
    end
end
end

function id = errid(fn)
id = '';
try, fn(); catch ME, id = ME.identifier; end
end

function tf = bits(a, b)
tf = isequal(size(a), size(b)) && isequal(typecast(a(:), 'uint64'), typecast(b(:), 'uint64'));
end

function tf = all_bits(A, B)
tf = true;
for r = 1:numel(A.roots)
    tf = tf && bits(A.roots(r).representative_profile, B.roots(r).representative_profile);
    for q = 1:numel(A.roots(r).discoveries)
        tf = tf && bits(A.roots(r).discoveries(q).returned_profile, B.roots(r).discoveries(q).returned_profile) && ...
            bits(A.roots(r).discoveries(q).comparison_profile, B.roots(r).discoveries(q).comparison_profile);
    end
end
end

function s = signbit_of(x)
s = double(bitshift(typecast(x, 'uint64'), -63));
end
