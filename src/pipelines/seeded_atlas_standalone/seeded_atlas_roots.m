function R = seeded_atlas_roots(run_dir, opts)
%SEEDED_ATLAS_ROOTS  The roots (quasi-steady states) of a run folder, loaded the
%   verified way (see docs/USER_GUIDE.md).
%
%   R = SEEDED_ATLAS_ROOTS(run_dir)
%   R = SEEDED_ATLAS_ROOTS(run_dir, struct('verify', false, 'sol', true))
%
%   verify (default true): every shard through sa_store_io('load_shard')
%   with the run's re-derivation and profile functions (verdict
%   re-derivation, kernel replay, profile check), as 'load_run_shards'
%   does; false reads the shards decoded but unverified (read_shard;
%   analysis only). Never raw-load a shard: deduplicated profiles are
%   stored empty.
%   One shard at a time: each shard is loaded (and verified), reduced at
%   once to the fields listed under R below, and released before the next
%   one is read, so memory holds one full shard plus the reduced list, not
%   every full store of the run (a production run's shards do not fit in
%   memory together). The shards are read in the same order and through
%   the same calls as a load of every shard into one map, so the result
%   and any error are the same; a store that cannot be reduced (a
%   malformed root) is kept whole and meets its error where it always did.
%   sol (default false): also load each representative attempt's iterate
%   (after every shard is read; one file per root, all kept in R).
%   run_dir is a run folder of either layout (sa_layout): its state folder
%   is read; nothing is written. A compacted folder loads unchanged (every
%   discovery's attempt file is kept).
%
%   R: struct array, lattice order, one element per root: root_id, key,
%   index, coordinate, routes, observed_family, n_discoveries,
%   representative_attempt_id, comparison_mesh, representative_profile
%   (identity rows on the comparison mesh), sol (if requested).

if nargin < 2, opts = struct(); end
run_dir = sa_layout('state', run_dir);           % internal/ or a flat folder
verify = ~isfield(opts, 'verify') || opts.verify;
want_sol = isfield(opts, 'sol') && opts.sol;
E = load(fullfile(run_dir, 'settings_effective.mat'), 'settings', 'lattice');
E.settings = sa_settings('legacy_fill', E.settings);   % older folders load with neutral leaves
M = sa_model(E.settings.model);
A = load(fullfile(run_dir, 'attempts.mat'), 'attempts');
ids = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:numel(A.attempts), ids(A.attempts(k).attempt_id) = k; end
if verify                                         % load_run_shards' order: re-derivation, profiles, listing
    fn = sa_store_io('run_rederive', run_dir, M);
    pf = sa_store_io('run_profiles', run_dir, M);
end
Sh = containers.Map('KeyType', 'char', 'ValueType', 'any');   % node key -> reduced roots (or the whole store)
Rd = containers.Map('KeyType', 'char', 'ValueType', 'logical');  % node key -> reduced
f = sa_store_io('shard_files', run_dir);           % skips sa_write temporary files
for k = 1:numel(f)
    if verify
        [st, key] = sa_store_io('load_shard', fullfile(f(k).folder, f(k).name), M, fn, pf);
    else
        Z = sa_store_io('read_shard', fullfile(f(k).folder, f(k).name));
        key = Z.key; st = Z.store;
        Z = []; %#ok<NASGU>
    end
    [v, ok] = reduce(st);
    Sh(key) = v;                                   % a later file of the same key replaces it, as in a map
    Rd(key) = ok;
    v = []; %#ok<NASGU>
    st = []; %#ok<NASGU>                           % the full store is released before the next read
end
keys = Sh.keys;
pos = cellfun(@(k) find(strcmp(E.lattice.key, k), 1), keys);
[~, o] = sort(pos);
R = struct('root_id', {}, 'key', {}, 'index', {}, 'coordinate', {}, 'routes', {}, ...
    'observed_family', {}, 'n_discoveries', {}, 'representative_attempt_id', {}, ...
    'comparison_mesh', {}, 'representative_profile', {}, 'sol', {});
for k = o
    st = Sh(keys{k});
    q = pos(k);
    if Rd(keys{k})                                 % reduced: the fields below, read from the full store
        for r = 1:numel(st)
            rr = st(r);
            ob = '';
            if isKey(ids, rr.representative_attempt_id)
                ob = A.attempts(ids(rr.representative_attempt_id)).observed_family;
            end
            sol = [];
            if want_sol
                Z = load(fullfile(run_dir, 'attempts', [rr.representative_attempt_id '.mat']), 'sol');
                sol = Z.sol;
            end
            R(end+1) = struct('root_id', sprintf('%s/%s', keys{k}, rr.root_id), 'key', keys{k}, ...
                'index', E.lattice.index(q, :), 'coordinate', E.lattice.coordinate(q, :), ...
                'routes', {rr.routes}, 'observed_family', ob, ...
                'n_discoveries', rr.n_discoveries, ...
                'representative_attempt_id', rr.representative_attempt_id, ...
                'comparison_mesh', rr.comparison_mesh, 'representative_profile', rr.representative_profile, ...
                'sol', sol); %#ok<AGROW>
        end
        continue
    end
    for r = 1:numel(st.roots)                      % a store kept whole (it did not reduce)
        rr = st.roots(r);
        ob = '';
        if isKey(ids, rr.representative_attempt_id)
            ob = A.attempts(ids(rr.representative_attempt_id)).observed_family;
        end
        sol = [];
        if want_sol
            Z = load(fullfile(run_dir, 'attempts', [rr.representative_attempt_id '.mat']), 'sol');
            sol = Z.sol;
        end
        R(end+1) = struct('root_id', sprintf('%s/%s', keys{k}, rr.root_id), 'key', keys{k}, ...
            'index', E.lattice.index(q, :), 'coordinate', E.lattice.coordinate(q, :), ...
            'routes', {unique({rr.discoveries.route})}, 'observed_family', ob, ...
            'n_discoveries', numel(rr.discoveries), ...
            'representative_attempt_id', rr.representative_attempt_id, ...
            'comparison_mesh', rr.comparison_mesh, 'representative_profile', rr.representative_profile, ...
            'sol', sol); %#ok<AGROW>
    end
end
end

function [v, ok] = reduce(st)
% The fields the root list reads, root by root (ok = true): the
% values are the store's own (routes = unique of the discovery routes, as
% below), so the list built from them is the one built from the full store.
% A store that does not reduce (a missing field, a route that unique
% refuses) is returned whole (ok = false): the loop above then meets the same error at
% the same root as it always did.
try
    red = struct('root_id', {}, 'representative_attempt_id', {}, 'routes', {}, 'n_discoveries', {}, ...
        'comparison_mesh', {}, 'representative_profile', {});
    for r = 1:numel(st.roots)
        rr = st.roots(r);
        red(r).root_id = rr.root_id;
        red(r).representative_attempt_id = rr.representative_attempt_id;
        red(r).routes = unique({rr.discoveries.route});
        red(r).n_discoveries = numel(rr.discoveries);
        red(r).comparison_mesh = rr.comparison_mesh;
        red(r).representative_profile = rr.representative_profile;
    end
    v = red; ok = true;
catch
    v = st; ok = false;
end
end
