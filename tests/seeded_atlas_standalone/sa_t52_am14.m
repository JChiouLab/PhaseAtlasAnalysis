function rows = sa_t52_am14(E)
%SA_T52_AM14  Milestone-4 T52 (packet v0.3 AM-14): load_shard recomputes every
%   discovery's returned mesh, returned profile and comparison profile from
%   its attempt's sol and compares bits.
%   (a) a CONSISTENTLY tampered shard (re-built through the kernel, so the
%       AM-1 replay passes) is refused by load_shard alone (SA:STORE_IO:
%       E_PROFILE), with the profile check built from the run folder;
%   (b) the untampered shards of the same run pass, every row bit-equal;
%   (c) real shards (copies) of K2' (G3), K6 and a deterministic sample of K3
%       pass with the run's re-derivation and profile functions; they were
%       written under software 0.2, so the check is version-scoped: reported
%       (largest difference), never failed;
%   (d) pre-AM-9 shard-1 files pass;
%   (e) an older-version record whose profile differs is reported, not failed.
%   E.am14_real = false skips (c)-(e) (they read milestone-3 folders).
T = 'T52 AM-14 shard profile check';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't52');
U = fullfile(d, 'U');
sa_m3_run(U, [2 1]);
M = sa_model('stub_two_pool');
fn = sa_store_io('run_rederive', sa_layout('state', U), M);
f = dir(fullfile(sa_layout('state', U), 'roots', '*.mat'));
ok = true; nrow = 0; alleq = true;
for k = 1:numel(f)
    [~, ~, ~, pr] = sa_store_io('load_shard', fullfile(f(k).folder, f(k).name), M, fn);
    ok = ok && pr.checked; nrow = nrow + numel(pr.rows);
    alleq = alleq && all([pr.rows.equal]) && all([pr.rows.same_version]);
end
rows = sa_check(rows, T, '(b) untampered shards: checked from the run folder, every discovery bit-equal to its recomputation', ...
    ok && alleq && nrow > 0, sprintf('%d shards, %d discoveries', numel(f), nrow));
% (a) consistent tamper
file = fullfile(f(1).folder, f(1).name);
Z = sa_store_io('read_shard', file);
st = tamper(Z.store, 1.001);
sa_store_io('save_shard', st, file, Z.key, Z.attempt_count);
got = '';
try
    sa_store_io('load_shard', file, M, fn);
catch ME
    got = ME.identifier;
end
got_noprof = '';
try
    sa_store_io('load_shard', file, M, fn, []);
catch ME
    got_noprof = ME.identifier;
end
rows = sa_check(rows, T, ['(a) a consistently tampered shard (kernel replay passes) is refused by load_shard alone ' ...
    '(SA:STORE_IO:E_PROFILE); with the check explicitly disabled it would load'], ...
    strcmp(got, 'SA:STORE_IO:E_PROFILE') && isempty(got_noprof), [got ' / ' got_noprof]);
rmpath(E.stub);
if isfield(E, 'am14_real') && ~E.am14_real, return; end
% (c) real shards (copies)
m3 = fullfile(E.root, 'data', 'experiments', 'standalone_seeded_atlas_20261012', 'm3');
sets = {'K2p', fullfile(m3, 'partG3_20261013_125459', 'K2p_A'), 'guzzo_noxa', 1; ...
    'K6', fullfile(m3, 'partK6_20261013_134209', 'K6'), 'original_polarity', 1; ...
    'K3', fullfile(m3, 'partK3_20261013_151020', 'K3'), 'guzzo_noxa', 12};
for s = 1:size(sets, 1)
    [ok, det] = real_set(d, sets{s, 1}, sets{s, 2}, sets{s, 3}, sets{s, 4});
    rows = sa_check(rows, T, sprintf('(c) real %s shards (copies; every %d-th shard by name): load with the profile check; version-scoped report', ...
        sets{s, 1}, sets{s, 4}), ok, det);
end
% (d) shard-1 files
sd = fullfile(m3, 'partK2_20261013_113221', 'A');
[ok, det] = real_set(d, 'K2_shard1', sd, 'guzzo_noxa', 1);
f1 = dir(fullfile(sa_layout('state', sd), 'roots', '*.mat'));
Zf = load(fullfile(f1(1).folder, f1(1).name), 'format');
rows = sa_check(rows, T, '(d) pre-AM-9 shard-1 files (K2 run A) load with the profile check', ...
    ok && strcmp(Zf.format, 'sa_store_io/shard-1'), [Zf.format '; ' det]);
% (e) an older-version record with a differing profile: reported, not failed
src = sets{1, 2};
Mn = sa_model('guzzo_noxa');
g = dir(fullfile(sa_layout('state', src), 'roots', '*.mat'));
cp = fullfile(d, 'old_tamper.mat');
Z = sa_store_io('read_shard', fullfile(g(1).folder, g(1).name));
sa_store_io('save_shard', tamper(Z.store, 1.001), cp, Z.key, Z.attempt_count);
got = '';
try
    [~, ~, ~, pr] = sa_store_io('load_shard', cp, Mn, @(dd) deal(true, {}), sa_store_io('run_profiles', sa_layout('state', src), Mn));
catch ME
    got = ME.identifier;
end
rows = sa_check(rows, T, '(e) older-version record (software 0.2) with a differing profile: reported with its largest difference, not failed', ...
    isempty(got) && ~any([pr.rows.same_version]) && ~all([pr.rows.equal]) && max([pr.rows.max_abs_diff]) > 0, ...
    sprintf('%s max diff %.3g', got, max([pr.rows.max_abs_diff])));
end

function st2 = tamper(st, factor)
% rebuild the store through the kernel with every returned profile scaled:
% the AM-1 replay (rebuild) then matches; only AM-14 sees the change
disc = [];
for r = 1:numel(st.roots)
    for q = 1:numel(st.roots(r).discoveries)
        dd = st.roots(r).discoveries(q);
        disc = [disc; struct('n', dd.discovery_number, 'key', st.roots(r).target_key, 'd', dd)]; %#ok<AGROW>
    end
end
[~, o] = sort([disc.n]); disc = disc(o);
st2 = sa_root_store('empty', st.settings);
for k = 1:numel(disc)
    dd = disc(k).d;
    c = struct('target_key', disc(k).key, 'coordinate', dd.coordinate, 'mesh', dd.returned_mesh, ...
        'profile', dd.returned_profile * factor, 'seed_family', dd.seed_family, ...
        'observed_family', dd.observed_family, 'route', dd.route, 'attempt_id', dd.attempt_id, ...
        'valid', true, 'validity', dd.validity);
    st2 = sa_root_store('add', st2, c);
end
end

function [ok, det] = real_set(d, name, src, model, every)
M = sa_model(model);
fn = sa_store_io('run_rederive', sa_layout('state', src), M);
pf = sa_store_io('run_profiles', sa_layout('state', src), M);
f = dir(fullfile(sa_layout('state', src), 'roots', '*.mat'));
names = sort({f.name});
names = names(1:every:end);
cd0 = fullfile(d, ['copies_' name]); mkdir(cd0);
ok = true; n = 0; neq = 0; same = 0; mx = 0; bad = '';
for k = 1:numel(names)
    c = fullfile(cd0, names{k});
    copyfile(fullfile(sa_layout('state', src), 'roots', names{k}), c);
    try
        [~, ~, ~, pr] = sa_store_io('load_shard', c, M, fn, pf);
        n = n + numel(pr.rows); neq = neq + sum([pr.rows.equal]); same = same + sum([pr.rows.same_version]);
        mx = max([mx, pr.rows.max_abs_diff]);
        ok = ok && pr.checked;
    catch ME
        ok = false; bad = [names{k} ': ' ME.identifier ' ' ME.message];
        break
    end
end
rmdir(cd0, 's');
det = sprintf('%d shards, %d discoveries, %d bit-equal, %d same-version, largest difference %.3g %s', ...
    numel(names), n, neq, same, mx, bad);
end
