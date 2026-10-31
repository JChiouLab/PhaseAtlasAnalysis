function rows = sa_t80_render_cache(E)
%SA_T80_RENDER_CACHE  The cached rendering of the human-facing files under
%   report.files_interval_s = 1e-6 (attempts.csv, summary.txt and
%   roots_index.csv rewritten after every live call from the render cache),
%   on stub runs that exercise multi-root nodes (two families cold at four
%   targets, box [3 2], the per-root pass), the rescue ladder (the sa_m3_wall
%   refusals at k1 = 2, box [3 1]) and the per-root pass alone (branch B
%   spreads only through it). Each configuration runs uninterrupted (U), and
%   killed twice inside a call (SA:TEST:KILL: journal start, no attempt
%   file -> an interrupted record) and resumed each time (K). Zero bvp4c calls.
%   Per configuration: no report self-check note in either journal (the
%   cached render equalled the full render at every stop); the files were
%   written after the last live call before each kill; the light result
%   stores made full (sa_store_io 'full_map') equal the shard files;
%   result.shards is light; K's interrupted rows in attempts.csv name the
%   call that re-ran them; K == U (sa_m3_equiv, interrupted records
%   permitted); U's roots_index.csv equals the one a run with the default
%   interval writes at its stop (render only at the stop).
T = 'T80 render cache (files after every call, two kills)';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't80');
rep = struct('report', struct('files_interval_s', 1e-6));
wall = @(dd, b, ex) sa_m3_wall(dd, ex);
noforce = @(dd, b) '';
o_rp = struct('families', {{'X_X', 'Y_X'}}, 'targets', struct('list', [0 0]), ...
    'continuation', struct('starts', 'anchor'), 'routes', struct('rescue', false), ...
    'root_pass', struct('enabled', true));
o_mrn = struct('families', {{'X_X', 'Y_X'}}, 'targets', struct('list', [0 0; 2 0; -2 1; 1 1]), ...
    'root_pass', struct('enabled', true));
cfg = {
    'multiroot', @(dd, x, varargin) sa_m3_run(dd, [3 2], 'settings', mrg(o_mrn, x), varargin{:})
    'rescue',    @(dd, x, varargin) sa_m3_run(dd, [3 1], 'settings', x, 'fail', wall, varargin{:})
    'rootpass',  @(dd, x, varargin) sa_m3_run(dd, [2 2], 'settings', mrg(o_rp, x), 'force', noforce, varargin{:})
    };
for c = 1:size(cfg, 1)
    nm = cfg{c, 1}; run = cfg{c, 2};
    U = fullfile(d, ['U_' nm]); U0 = fullfile(d, ['U0_' nm]); K = fullfile(d, ['K_' nm]);
    rU = run(U, rep);
    rU0 = run(U0, struct());
    A = rU.attempts;
    nmr = 0;
    for k = rU.shards.keys, nmr = nmr + (numel(rU.shards(k{1}).roots) >= 2); end
    ex = sprintf(['%d calls, %d roots, %d multi-root nodes; accepted: continuation %d, rescue %d, per-root ' ...
        'pass %d, matched %d'], rU.calls, rU.n_roots, nmr, cnt(A, 'continuation'), ...
        sum(strcmp({A.outcome}, 'accepted') & startsWith({A.route}, 'rescue')), cnt(A, 'root_pass'), ...
        sum(strcmp({A.store_status}, 'matched_existing')));
    % ---- two kills inside a call, a resume after each ------------------------
    k1 = round(rU.calls * 0.4); k2 = round(rU.calls * 0.75);
    e1 = kill(run, K, rep, k1);
    n1 = csv_rows(K);                       % written after call k1 - 1
    x = rep; x.resume = true;
    e2 = kill(run, K, x, k2 - k1 + 1);      % live calls restart at 1 = attempt k1 + 1 (the retry of k1):
    n2 = csv_rows(K);                       % killed inside attempt k2 + 1, written after attempt k2
    rK = run(K, x);
    D = sa_m3_equiv(U, K);
    % ---- the checks ------------------------------------------------------------
    sc = self_check_notes(U) + self_check_notes(K) + self_check_notes(U0);
    [ok_fU, nfU] = refill_equals_shards(rU);
    [ok_fK, nfK] = refill_equals_shards(rK);
    light = all_light(rU) && all_light(rK);
    [ok_int, det_int] = interrupted_text(K, rK.attempts);
    ri_eq = strcmp(fileread(fullfile(U, 'roots_index.csv')), fileread(fullfile(U0, 'roots_index.csv')));
    has_variety = (c ~= 1 || nmr >= 2) && (c ~= 2 || any(startsWith({A.route}, 'rescue') & ...
        strcmp({A.outcome}, 'accepted'))) && (c ~= 3 || cnt(A, 'root_pass') > 0);
    rows = sa_check(rows, T, sprintf(['%s: no report self-check note (U, K, default-interval run); files written ' ...
        'after the last live call before each kill (attempts.csv rows %d = k1 - 1, and %d = k2 with the first ' ...
        'interrupted row); the run exercises its feature'], nm, n1, n2), ...
        sc == 0 && strcmp(e1, 'SA:TEST:KILL') && strcmp(e2, 'SA:TEST:KILL') && n1 == k1 - 1 && n2 == k2 && ...
        has_variety, sprintf('%s; self-check notes %d; kills %s, %s at %d and %d', ex, sc, e1, e2, k1, k2));
    rows = sa_check(rows, T, sprintf(['%s: the light result stores made full equal the shard files (U %d nodes, ' ...
        'K %d nodes); result.shards light'], nm, nfU, nfK), ok_fU && ok_fK && light && nfU > 0, '');
    rows = sa_check(rows, T, sprintf(['%s: K''s interrupted rows in attempts.csv name the call that re-ran ' ...
        'them; K == U (sa_m3_equiv: %d units, %d interrupted)'], nm, D.n_units, D.n_interrupted), ...
        ok_int && D.pass && D.n_interrupted == 2 && strcmp(rK.stop_reason, rU.stop_reason), ...
        sprintf('%s; %s', det_int, strjoin(D.diffs, ' | ')));
    rows = sa_check(rows, T, sprintf(['%s: roots_index.csv rendered after every call equals the one rendered ' ...
        'only at the stop (default interval)'], nm), ri_eq && rU0.calls == rU.calls, '');
end
rmpath(E.stub);
end

function id = kill(run, K, x, k)
id = '';
try
    run(K, x, 'kill_at', k);
catch ME
    id = ME.identifier;
end
end

function n = csv_rows(K)
L = strsplit(strtrim(fileread(fullfile(K, 'attempts.csv'))), newline);
n = numel(L) - 1;
end

function n = self_check_notes(d)
n = numel(regexp(fileread(fullfile(sa_layout('state', d), 'journal.txt')), 'note report self-check', 'match'));
end

function [ok, n] = refill_equals_shards(r)
F = sa_store_io('full_map', r);
ok = true; n = 0;
for k = F.keys
    Z = sa_store_io('read_shard', fullfile(sa_layout('state', r.run_dir), 'roots', sa_store_io('shard_file', k{1})));
    ok = ok && isequaln(Z.store, F(k{1}));
    n = n + 1;
end
f = dir(fullfile(sa_layout('state', r.run_dir), 'roots', '*.mat'));
ok = ok && numel(f) == F.Count;
end

function tf = all_light(r)
tf = true;
for k = r.shards.keys
    s = r.shards(k{1});
    if ~isempty(s.roots), tf = tf && sa_store_io('is_light', s); end
    for q = 1:numel(s.roots)
        tf = tf && all(arrayfun(@(dd) isempty(dd.returned_mesh) && isempty(dd.returned_profile), ...
            s.roots(q).discoveries));
    end
end
end

function [ok, det] = interrupted_text(K, A)
L = strsplit(strtrim(fileread(fullfile(K, 'attempts.csv'))), newline);
ints = find(strcmp({A.outcome}, 'interrupted'));
ok = numel(ints) == 2; det = {};
for q = ints
    id = A(q).attempt_id;
    by = find(strcmp({A.retry_of}, id));
    line = L{1 + q};
    good = numel(by) == 1 && startsWith(line, ['"' id '"']) && contains(line, sprintf('it was run again as %s', A(by).attempt_id)) && ...
        ~contains(line, 'the next call of the same work unit');
    ok = ok && good;
    if numel(by) == 1, det{end+1} = sprintf('%s run again as %s: %d', id, A(by).attempt_id, good); end %#ok<AGROW>
end
det = strjoin(det, ', ');
end

function n = cnt(A, route)
n = sum(strcmp({A.outcome}, 'accepted') & strcmp({A.route}, route));
end

function a = mrg(a, b)
f = fieldnames(b);
for k = 1:numel(f)
    if isstruct(b.(f{k})) && isfield(a, f{k}) && isstruct(a.(f{k}))
        a.(f{k}) = mrg(a.(f{k}), b.(f{k}));
    else
        a.(f{k}) = b.(f{k});
    end
end
end
