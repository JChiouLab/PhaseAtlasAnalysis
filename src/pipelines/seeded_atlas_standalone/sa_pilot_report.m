function varargout = sa_pilot_report(run_dir, varargin)
%SA_PILOT_REPORT  The pilot report: a pure
%   function of a run folder (a cost and diagnostic summary of a small run,
%   used to plan a larger one).
%   Terms: NMax-bound = the call reached the bvp4c mesh-point limit NMax;
%   reach map = the lattice drawn by node class; frontier table = where each
%   refused continuation step failed; MATLAB-h = MATLAB hours; edge = one
%   continuation step between neighbouring grid points; ring = grid points
%   outside the library's native range; cut_by_cap = refused on a mesh bound
%   while the cold mesh cap was tighter than the run's solver. Zero solves; rebuildable; deterministic; runs
%   on any stop reason (budget_calls, budget_wall, disk_floor, completed).
%
%   SA_PILOT_REPORT(run_dir)   writes <run_dir>/pilot_report.txt and the
%                              section tables pilot/*.csv in the state folder
%                              (sa_layout: <run_dir>/internal/pilot, or
%                              <run_dir>/pilot for a flat folder); refused on
%                              a compacted folder (it counts attempt files)
%   R = SA_PILOT_REPORT(run_dir)   also returns the computed numbers
%   E = SA_PILOT_REPORT(run_dir, 'large_change')   only the edge_large_change section,
%       returned, NOTHING written (works on older folders read only:
%       step distances are recomputed from the attempt files when an edge
%       carries no d_sub). SA_PILOT_REPORT(run_dir, 'large_change', Q) uses
%       the scheduler state Q instead of the folder's queue_final/queue.mat.
%
%   Sections (each also a CSV under pilot/): 1 settings (non-default leaves
%   and source, stop reason, calls, process time); 2 validity and refusals
%   by check (first-failing and all-failing per route class; ring split;
%   no_return by guard; nodal-only positivity); 3 cost (solver s per call for valid,
%   invalid and NMax-bound calls, median and p90; process s per live call;
%   NMax share by route; MB per attempt file by class and per root shard;
%   replay s per record); 4 returned mesh by call-chain depth (valid and
%   invalid) and guess sizes; 5 route yield (nodes reached per route and
%   rescue rung, calls per reached node, calls per refused edge, rescue
%   calls and yield by region and rung, repeated calls = attempts whose
%   guess fingerprint equals an earlier one's); 6 reach map (ASCII lattice
%   of node classes, ring in lower case; reach_map.csv) and the frontier
%   table (per refused edge: coordinate of the last valid substep and of the
%   failing call, the failing check); 7 identity (near-threshold band,
%   same-target distinct pairs, nodes with >= 2 roots;
%   edge_large_change with the attempt-level shape switch and the
%   second-difference diagnostic); 8 production projection (calls,
%   MATLAB-h, disk: central / high / cap, formulas printed).
%
%   Report-only rules (no gate, stored value or verdict): the
%   reach fraction f is taken over the whole planned region (reached / region
%   nodes; refused and never-attempted nodes are both unreached); a report of
%   a run that stopped before completing is marked PARTIAL, counts its
%   pending region nodes as unreached and, while no refused node has been
%   measured, charges each unreached node the per-node bound (2h + 1) +
%   rescue.max_calls_per_node; after a resume the process/solver factor
%   is measured on the live calls only, with the replayed records' time
%   (records x report.replay_s_per_record, a model setting) removed;
%   section 1 shows the current segment's budget and checkpoint leaves
%   (resume_<n>.mat); no_return is tabulated by guard and error identifier
%   (hit counts summed, not part of the key); all 11 tables are written at
%   every stop, including 0 calls (values not measured print NaN); the
%   projected roots use the same f.
%
%   Report-only: when the cold targets are
%   not all_seedable (rule: settings targets.policy is anything but
%   'all_seedable'), section 7 prints MULTIPLICITY NOT TESTED (roots are
%   those connected to the cold-seeded nodes only; coexisting roots
%   elsewhere are not looked for; the multi-root count is a lower bound;
%   per_seed: NOT TESTED BY COLD; only continuation and the pass look),
%   and section 8's central_nine alternative carries the same caveat.
%   Section 8's cold count is the per_seed seed count (the production
%   default; census.seed_counts) unless the folder's policy is all_seedable
%   or the folder was written without the seed counts (then all_seedable,
%   thinned).
%   Repeated calls exclude the retry of an interrupted record: a call with
%   retry_of set (or a killed retry, itself an interrupted record) whose
%   earlier equal-fingerprint calls are all interrupted records is counted
%   separately, never as a repeat.
%
%   Report-only, no formula involved: section 2
%   adds the reason class cut_by_cap per route (cold calls refused on a mesh
%   bound under cold.mesh_cap); section 3 splits NMax-bound calls into "at
%   the cold cap" and "at solver.nmax"; section 5 adds the per-root pass
%   (route root_pass: calls, new roots, matched, refused, calls per new
%   root; pilot/root_pass.csv = the rounds table, Q.root_pass), the rescue
%   no_basis count, and NOTE lines (the pass posed n targets in the
%   pilot region: a cross pilot underprices it; the cap column omits it);
%   section 8 uses the thinned cold-target count when targets.thin is set
%   and says that b_bound is conservative under secant-only rescue (the
%   secant-only per-node bound is 2 x n_pools calls).
%   The section-3 split uses the recorded
%   largest proposed mesh P (P > the call's nmax) where present, and says
%   for how many calls; rows from older versions keep the convergence class.
%   Older folders load through sa_settings('legacy_fill') (neutral values).
%
%   edge_large_change: d = d(parent, child) for a 1-call edge and
%   max(d(parent, child), max over substeps of d_sub / s_sub) for a
%   multi-call edge; edge_large_change := d > report.large_change_factor x
%   the median d of the run's completed 1-call interior continuation edges
%   (at least 30 of them, else not_computed). Report-only: the stored
%   large_change stays as recorded. Second difference on collinear lineage
%   chains: d(child, 2 parent - grandparent) / d(parent, child).

mode = 'all';
if ~isempty(varargin), mode = varargin{1}; end
sa_report('refuse_if_compacted', run_dir, 'sa_report(''pilot'')');
X = load_run(sa_layout('state', run_dir), ~strcmp(mode, 'large_change'));
if strcmp(mode, 'large_change')
    if numel(varargin) >= 2 && ~isempty(varargin{2}), X.Q = varargin{2}; end   % e.g. a stored final queue
    varargout{1} = large_change(X);
    return
end
R = struct();
out = {};
out{end+1} = sprintf('PILOT REPORT  %s', run_dir);
out{end+1} = sprintf('(sa_pilot_report, %s; zero solves; a pure function of the run folder)', ...
    X.S.software_version);
[out, R.settings] = sec_settings(X, out);
[out, R.validity] = sec_validity(X, out);
[out, R.cost] = sec_cost(X, out);
[out, R.mesh] = sec_mesh(X, out);
[out, R.yield] = sec_yield(X, out);
[out, R.reach] = sec_reach(X, out);
[out, R.identity, R.large_change] = sec_identity(X, out);
[out, R.projection] = sec_projection(X, R, out);
sa_write('text', fullfile(run_dir, 'pilot_report.txt'), [strjoin(out, newline) newline]);
if nargout > 0, varargout{1} = R; end
end

% =========================================================================
% inputs

function X = load_run(d, write)
X = struct('dir', d, 'write', write);
E = load(fullfile(d, 'settings_effective.mat'));
E.settings = sa_settings('legacy_fill', E.settings);   % neutral values for older folders
X.E = E; X.S = E.settings;
if ~isfield(X.S, 'report'), X.S.report = struct('large_change_factor', 3); end   % older folders
if isfield(E, 'settings_sources'), X.SRC = E.settings_sources; else, X.SRC = struct(); end
A = load(fullfile(d, 'attempts.mat'));
X.A = A.attempts;
n = numel(X.A);
X.gn = nan(1, n); X.fp = repmat({''}, 1, n);
if isfield(A, 'guess_nodes'), m = min(n, numel(A.guess_nodes)); X.gn(1:m) = A.guess_nodes(1:m); end
if isfield(A, 'guess_fingerprint'), m = min(n, numel(A.guess_fingerprint)); X.fp(1:m) = A.guess_fingerprint(1:m); end
X.stop_reason = A.stop_reason; X.calls = A.calls;
X.process_s = NaN;
if exist(fullfile(d, 'process_time.mat'), 'file') == 2
    Zp = load(fullfile(d, 'process_time.mat')); X.process_s = Zp.process_s;
end
if exist(fullfile(d, 'queue_final.mat'), 'file') == 2
    Z = load(fullfile(d, 'queue_final.mat')); X.Q = Z.Q; X.queue_file = 'queue_final.mat';
elseif exist(fullfile(d, 'queue.mat'), 'file') == 2
    Z = load(fullfile(d, 'queue.mat')); X.Q = Z.Q; X.queue_file = 'queue.mat';
else
    X.Q = []; X.queue_file = '';
end
X.ids = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:n, X.ids(X.A(k).attempt_id) = k; end
X.L = E.lattice;
% roots (decoded shards; analysis read, never re-validated here)
X.roots = struct('key', {}, 'root_id', {}, 'routes', {}, 'n_disc', {}, 'rep', {});
X.shard_bytes = [];
f = sa_store_io('shard_files', d);            % roots/*.mat without sa_write's temporary files
for k = 1:numel(f)
    Z = sa_store_io('read_shard', fullfile(f(k).folder, f(k).name));
    X.shard_bytes(end+1) = f(k).bytes;
    for r = 1:numel(Z.store.roots)
        rr = Z.store.roots(r);
        X.roots(end+1) = struct('key', Z.key, 'root_id', sprintf('%s/%s', Z.key, rr.root_id), ...
            'routes', {unique({rr.discoveries.route})}, 'n_disc', numel(rr.discoveries), ...
            'rep', rr.representative_attempt_id);
    end
end
X.nroots = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:numel(X.roots)
    if isKey(X.nroots, X.roots(k).key), X.nroots(X.roots(k).key) = X.nroots(X.roots(k).key) + 1;
    else, X.nroots(X.roots(k).key) = 1; end
end
fa = dir(fullfile(d, 'attempts', 'A*.mat'));
X.file_bytes = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:numel(fa), X.file_bytes(fa(k).name(1:end-4)) = fa(k).bytes; end
% resumes: the records replayed by the last resume (resume_<n>.mat)
X.resumes = numel(dir(fullfile(d, 'resume_*.mat')));
X.replay = [];
if X.resumes > 0
    Zr = load(fullfile(d, sprintf('resume_%d.mat', X.resumes)), 'attempt_count');
    X.replay = Zr.attempt_count;
end
% the current segment's budget and checkpoint leaves: a
% resume may change them; settings_effective.mat holds the first segment's
X.seg = struct('budget', X.S.budget, 'checkpoint', X.S.checkpoint, 'n', 0);
if X.resumes > 0
    Zs = load(fullfile(d, sprintf('resume_%d.mat', X.resumes)), 'budget', 'checkpoint');
    if isfield(Zs, 'budget'), X.seg.budget = Zs.budget; end
    if isfield(Zs, 'checkpoint'), X.seg.checkpoint = Zs.checkpoint; end
    X.seg.n = X.resumes;
end
% replay seconds per stored record: the report leaf
% report.replay_s_per_record (model layer);
% a folder written without the leaf takes the current model layer's
% value, else the package default
[X.rate, X.rate_src] = replay_rate(X);
X.pairs = {};
if exist(fullfile(d, 'identity_pairs.csv'), 'file') == 2
    X.pairs = csv_rows(fullfile(d, 'identity_pairs.csv'));
end
if write && exist(fullfile(d, 'pilot'), 'dir') ~= 7, mkdir(fullfile(d, 'pilot')); end
end

function [r, src] = replay_rate(X)
if isfield(X.S.report, 'replay_s_per_record')
    r = X.S.report.replay_s_per_record; src = 'settings';
    if isfield(X.SRC, 'report') && isfield(X.SRC.report, 'replay_s_per_record')
        src = X.SRC.report.replay_s_per_record;
    end
    return
end
D = sa_settings('defaults');
r = D.report.replay_s_per_record; src = 'package default (leaf absent in this folder)';
M = sa_model(X.S.model);                % as large_change does
if isfield(M, 'setting_defaults') && isfield(M.setting_defaults, 'report') && ...
        isfield(M.setting_defaults.report, 'replay_s_per_record')
    r = M.setting_defaults.report.replay_s_per_record;
    src = 'model layer of the current code (leaf absent in this folder)';
end
end

function C = csv_rows(f)
t = strsplit(strtrim(fileread(f)), newline);
C = cellfun(@(s) strrep(regexp(s, ',(?=(?:[^"]*"[^"]*")*[^"]*$)', 'split'), '"', ''), t, ...
    'UniformOutput', false);
end

function write_csv(X, name, header, rows)
if ~X.write, return; end
lines = [{strjoin(header, ',')}, cellfun(@(r) strjoin(cellfun(@cell_text, r, 'UniformOutput', false), ','), ...
    rows, 'UniformOutput', false)];
sa_write('text', fullfile(X.dir, 'pilot', name), [strjoin(lines, newline) newline]);
end

function t = cell_text(v)
if ischar(v)
    if any(v == ',') || any(v == '"'), t = ['"' strrep(v, '"', '""') '"']; else, t = v; end
elseif islogical(v), t = sprintf('%d', v);
elseif isempty(v), t = '';
else, t = sprintf('%.6g', v);
end
end

function c = route_class(a)
% cold | refinement | import | continuation | substep | rescue_<rung>
% (substeps of a rescue edge are rescue_<rung>_substep; import: the import route)
if strcmp(a.kind, 'cold') || strcmp(a.kind, 'refinement') || strcmp(a.kind, 'import'), c = a.kind; return; end
r = a.route; if isempty(r), r = a.kind; end
if isfinite(a.fraction) && a.fraction < 1
    if strcmp(r, 'continuation'), c = 'substep'; else, c = [r '_substep']; end
else
    c = r;
end
end

function tf = is_nmax(a)
tf = strcmp(a.convergence_class, 'UNRESOLVED_RELTOL_NMAX');
end

function tf = returned(a)
tf = ~any(strcmp(a.outcome, {'no_return', 'interrupted'}));
end

function v = q90(x)
x = sort(x(isfinite(x)));
if isempty(x), v = NaN; return; end
v = x(max(1, ceil(0.9 * numel(x))));
end

function v = med(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = median(x); end
end

% =========================================================================
% 1 settings

function [out, R] = sec_settings(X, out)
out{end+1} = ''; out{end+1} = '1. SETTINGS';
S = X.S;
seg = '';
if X.seg.n > 0
    seg = sprintf('; segment %d after resume %d; first segment''s budget %d', X.seg.n + 1, X.seg.n, ...
        S.budget.max_calls);
end
out{end+1} = sprintf('model %s  mode %s  stop reason %s  calls %d (budget %d this segment%s)  process %.1f s (this segment)', ...
    S.model, S.mode, X.stop_reason, X.calls, X.seg.budget.max_calls, seg, X.process_s);
if isfield(S.lattice, 'points_per_axis') && ~isempty(S.lattice.points_per_axis)
    out{end+1} = sprintf('lattice step %s (derived from points_per_axis %s)', mat2str(S.lattice.step, 6), ...
        mat2str(S.lattice.points_per_axis));
end
rows = {};
if ~isempty(fieldnames(X.SRC))
    % budget.* and checkpoint.* as in force in this segment (resume_<n>.mat);
    % a leaf a resume changed is sourced 'resume_<n>'
    SRC = X.SRC;
    for g = {'budget', 'checkpoint'}
        for f = fieldnames(X.seg.(g{1})).'
            if isfield(S.(g{1}), f{1}) && ~isequaln(S.(g{1}).(f{1}), X.seg.(g{1}).(f{1}))
                S.(g{1}).(f{1}) = X.seg.(g{1}).(f{1});
                SRC.(g{1}).(f{1}) = sprintf('resume_%d', X.seg.n);
            end
        end
    end
    L = sa_settings('effective_lines', S, SRC);
    for k = 1:numel(L)
        t = regexp(L{k}, '^(\S+)\s+(.*\S)\s+(\S+)$', 'tokens', 'once');
        if isempty(t) || any(strcmp(t{3}, {'package', 'derived'})), continue; end
        rows{end+1} = t; %#ok<AGROW>
        out{end+1} = sprintf('  %-40s %-30s %s', t{1}, t{2}, t{3}); %#ok<AGROW>
    end
else
    out{end+1} = '  (no settings_sources: the folder was written before setting sources were recorded)';
end
write_csv(X, 'settings.csv', {'leaf', 'value', 'source'}, rows);
R = struct('n_nondefault', numel(rows));
end

% =========================================================================
% 2 validity and refusals

function [out, R] = sec_validity(X, out)
out{end+1} = ''; out{end+1} = '2. VALIDITY AND REFUSALS BY CHECK (returned calls; first-failing / all-failing)';
A = X.A;
keys = {}; first = []; allc = [];
for k = 1:numel(A)
    a = A(k);
    if ~returned(a) || isempty(a.failing), continue; end
    c = strsplit(a.failing, ';');
    rc = route_class(a); ring = strcmp(a.phase, 'ring');
    for j = 1:numel(c)
        key = sprintf('%s|%s|%d', c{j}, rc, ring);
        i = find(strcmp(keys, key), 1);
        if isempty(i), keys{end+1} = key; first(end+1) = 0; allc(end+1) = 0; i = numel(keys); end %#ok<AGROW>
        allc(i) = allc(i) + 1;
        if j == 1, first(i) = first(i) + 1; end
    end
end
[keys, o] = sort(keys); first = first(o); allc = allc(o);
rows = cell(1, numel(keys));
for i = 1:numel(keys)
    p = strsplit(keys{i}, '|');
    rows{i} = {p{1}, p{2}, str2double(p{3}), first(i), allc(i)};
    out{end+1} = sprintf('  %-28s %-24s ring %s  first %5d  all %5d', p{1}, p{2}, p{3}, first(i), allc(i)); %#ok<AGROW>
end
if isempty(keys), out{end+1} = '  none'; end
% Reason class cut_by_cap per route (cold calls refused on a mesh bound
% while cold.mesh_cap was tighter; a superset of the calls the cap changed)
if isfield(A, 'cut_by_cap')
    cut = A(arrayfun(@(a) isequal(a.cut_by_cap, true), A));
    rcu = arrayfun(@(a) route_class(a), cut, 'UniformOutput', false);
    [u, ~, j] = unique(rcu);
    out{end+1} = sprintf(['  reason class cut_by_cap (refused on a mesh bound under the cold cap, cold.mesh_cap = %g; ' ...
        'not evidence that no root exists): %d call(s)'], X.S.cold.mesh_cap, numel(cut));
    for i = 1:numel(u)
        out{end+1} = sprintf('    %-24s %5d', u{i}, sum(j == i)); %#ok<AGROW>
        rows{end+1} = {'cut_by_cap', u{i}, 0, sum(j == i), sum(j == i)}; %#ok<AGROW>
    end
end
write_csv(X, 'validity_by_check.csv', {'check', 'route_class', 'ring', 'first_failing', 'all_failing'}, rows);
% no_return by guard
nr = A(strcmp({A.outcome}, 'no_return') | (startsWith({A.outcome}, 'substep') & ...
    strcmp({A.convergence_class}, 'NO_RETURN')));
% one row per guard and error identifier: the per-call hit
% count is summed, never part of the key
[g, hits] = cellfun(@(r) guard_key(r), {nr.corrector_reason}, 'UniformOutput', false);
hits = [hits{:}];
[u, ~, j] = unique(g);
out{end+1} = sprintf('  no_return (incl. substeps): %d calls; by guard and error identifier (calls, guard hits summed):', ...
    numel(nr));
rows = {};
for i = 1:numel(u)
    out{end+1} = sprintf('    %-64s calls %5d  hits %d', u{i}, sum(j == i), sum(hits(j == i))); %#ok<AGROW>
    rows{end+1} = {u{i}, sum(j == i), sum(hits(j == i))}; %#ok<AGROW>
end
write_csv(X, 'no_return_by_guard.csv', {'guard', 'calls', 'hits'}, rows);
rc5 = sum(arrayfun(@(a) a.nodal_min > 0 && ~(a.positivity_min > 0), A));
out{end+1} = sprintf('  positive at the mesh nodes but not between them (between-node positivity check): %d', rc5);
R = struct('n_keys', numel(keys), 'no_return', numel(nr), 'rc5', rc5);
end

function [t, hits] = guard_key(r)
% '<guard>: <n> hit(s); <id> <message>' -> '<guard> <id>' and n;
% 'bvp4c_error: <id> <message>' -> 'bvp4c_error <id>'; a reason without a
% guard prefix -> its text with numbers replaced by '#', at most 64 chars
hits = 0;
if isempty(r), t = '(no reason recorded)'; return; end
r = strtrim(regexprep(r, '\s+', ' '));
h = regexp(r, '(\d+) hit\(s\)', 'tokens', 'once');
if ~isempty(h), hits = str2double(h{1}); end
gname = regexp(r, '^([A-Za-z]\w*):\s', 'tokens', 'once');
eid = regexp(r, '[A-Za-z]\w*(?::\w+)+', 'match');
if ~isempty(gname), eid = eid(~strcmp(eid, gname{1})); end
if ~isempty(gname)
    t = gname{1};
    if ~isempty(eid), t = [t ' ' eid{1}]; end
else
    t = regexprep(r, '\d+(\.\d+)?', '#');
    if numel(t) > 64, t = t(1:64); end
end
end

% =========================================================================
% 3 cost

function [out, R] = sec_cost(X, out)
out{end+1} = ''; out{end+1} = '3. COST';
A = X.A;
ret = arrayfun(@returned, A);
valid = [A.valid] & ret;
nmax = arrayfun(@is_nmax, A);
inval = ret & ~valid & ~nmax;
w = [A.wall_s];
cls = {'valid', valid; 'invalid (not NMax-bound)', inval; 'NMax-bound', nmax; ...
    'no_return', ~ret & ~strcmp({A.outcome}, 'interrupted')};
rows = {};
R = struct();
for c = 1:size(cls, 1)
    m = cls{c, 2};
    mb = arrayfun(@(a) file_mb(X, a.attempt_id), A(m));
    out{end+1} = sprintf('  %-26s calls %5d  corrector s median %.3g  p90 %.3g  sum %.1f;  MB per attempt file mean %.3g', ...
        cls{c, 1}, sum(m), med(w(m)), q90(w(m)), sum(w(m), 'omitnan'), mean_or_nan(mb)); %#ok<AGROW>
    rows{end+1} = {cls{c, 1}, sum(m), med(w(m)), q90(w(m)), sum(w(m), 'omitnan'), mean_or_nan(mb)}; %#ok<AGROW>
end
R.s_valid = mean_or_nan(w(valid | inval | ~ret)); R.s_nmax = mean_or_nan(w(nmax));
R.mb_valid = mean_or_nan(arrayfun(@(a) file_mb(X, a.attempt_id), A(~nmax)));
R.mb_nmax = mean_or_nan(arrayfun(@(a) file_mb(X, a.attempt_id), A(nmax)));
R.nmax_share = sum(nmax) / max(1, numel(A));
% the process/solver factor on this segment's LIVE calls:
% after a resume the segment's process time also holds the replay of the
% stored records, removed as records x report.replay_s_per_record (the
% scheduler reads no clock during replay); never below the live solver s
nrep = 0; if ~isempty(X.replay), nrep = X.replay; end
live = numel(A) - nrep;
wl = sum(w(nrep + 1:end), 'omitnan');
R.replay_removed_s = 0; R.replay_capped = false;
proc_live = X.process_s;
if nrep > 0
    R.replay_capped = nrep * X.rate > X.process_s - wl;
    R.replay_removed_s = min(nrep * X.rate, max(0, X.process_s - wl));
    proc_live = X.process_s - R.replay_removed_s;
end
R.process_per_live = NaN; R.overhead = NaN;
if live > 0
    R.process_per_live = proc_live / live;
    if wl > 0, R.overhead = proc_live / wl; end
end
R.corrector_per_call = NaN;
if numel(A) > 0, R.corrector_per_call = sum(w, 'omitnan') / numel(A); end
if nrep > 0
    out{end+1} = sprintf(['  process s per live call %.3g (this segment: %.1f s, less %.1f s replay = %d records x %.3g s ' ...
        '[report.replay_s_per_record, %s], over %d live calls); process/corrector on live calls %.3g; corrector s per call %.3g'], ...
        R.process_per_live, X.process_s, R.replay_removed_s, nrep, X.rate, X.rate_src, live, R.overhead, ...
        R.corrector_per_call);
    if R.replay_capped
        out{end+1} = sprintf(['  NOTE: %d records x %.3g s exceeds this segment''s process time less its live corrector ' ...
            'time (%.1f s): the replay removed is capped there, the factor is 1 and the live overhead of this ' ...
            'segment is not identifiable (the setting is above this segment''s replay rate)'], nrep, X.rate, ...
            X.process_s - wl);
    end
else
    out{end+1} = sprintf('  process s per live call %.3g (this segment: %.1f s over %d live calls); process/corrector %.3g; corrector s per call %.3g', ...
        R.process_per_live, X.process_s, live, R.overhead, R.corrector_per_call);
end
% NMax-bound calls at the cold cap (the call's own nmax below
% solver.nmax) and at solver.nmax
if isfield(A, 'nmax')
    % where the call's largest proposed mesh P
    % is recorded, "NMax-bound" in this split is P > the call's own nmax (the
    % exact form); elsewhere the convergence class
    nb = nmax; ptxt = '';
    if isfield(A, 'max_proposed_mesh')
        hp = arrayfun(@(a) isfinite(a.max_proposed_mesh), A);
        if any(hp)
            nb(hp) = arrayfun(@(a) a.max_proposed_mesh > a.nmax, A(hp));
            ptxt = sprintf('; NMax-bound from the recorded proposals (P > the call''s nmax) for %d of %d calls', ...
                sum(hp), numel(A));
        end
    end
    atcap = nb & arrayfun(@(a) isfinite(a.nmax) && a.nmax < X.S.solver.nmax, A);
    out{end+1} = sprintf(['  NMax-bound calls split: at the cold cap (nmax %g) %d, corrector s sum %.1f; at solver.nmax ' ...
        '(%d) %d, corrector s sum %.1f%s'], X.S.cold.mesh_cap, sum(atcap), sum(w(atcap), 'omitnan'), X.S.solver.nmax, ...
        sum(nb & ~atcap), sum(w(nb & ~atcap), 'omitnan'), ptxt);
    rows{end+1} = {'nmax_bound at cold cap', sum(atcap), med(w(atcap)), q90(w(atcap)), sum(w(atcap), 'omitnan'), NaN};
    rows{end+1} = {'nmax_bound at solver.nmax', sum(nb & ~atcap), med(w(nb & ~atcap)), q90(w(nb & ~atcap)), ...
        sum(w(nb & ~atcap), 'omitnan'), NaN};
    R.nmax_at_cap = sum(atcap); R.nmax_at_solver = sum(nb & ~atcap);
end
rc = cellfun(@(a) route_class(a), num2cell(A), 'UniformOutput', false);
[u, ~, j] = unique(rc);
out{end+1} = '  NMax share by route class:';
for i = 1:numel(u)
    out{end+1} = sprintf('    %-26s %5d calls  NMax %4d (%.1f %%)', u{i}, sum(j == i), sum(nmax(j == i)), ...
        100 * sum(nmax(j == i)) / sum(j == i)); %#ok<AGROW>
    rows{end+1} = {['nmax_share ' u{i}], sum(j == i), sum(nmax(j == i)), NaN, NaN, NaN}; %#ok<AGROW>
end
nr = max(1, numel(X.roots));
R.mb_per_root = sum(X.shard_bytes) / 2^20 / nr;
out{end+1} = sprintf('  MB per root shard: %.3g MB per root (%d shards, %d roots); NMax share %.1f %%', ...
    R.mb_per_root, numel(X.shard_bytes), numel(X.roots), 100 * R.nmax_share);
if isempty(X.replay)
    out{end+1} = '  replay s per record: not resumed';
    R.replay_s = NaN;
else
    % the scheduler reads no clock during replay, so the replay cost is
    % bounded from the segment: (process s - live solver s) / replayed records
    R.replay_s = (X.process_s - wl) / max(1, X.replay);
    out{end+1} = sprintf(['  replay s per record <= %.3g (upper bound: last segment %.1f s minus %.1f s live corrector ' ...
        'time, over %d replayed records; includes set-up and live overhead); resumes %d; the setting ' ...
        'report.replay_s_per_record is %.3g s'], R.replay_s, X.process_s, wl, X.replay, X.resumes, X.rate);
end
write_csv(X, 'cost.csv', {'class', 'calls', 'median_s', 'p90_s', 'sum_s', 'mb_per_file'}, rows);
end

function v = file_mb(X, id)
if isKey(X.file_bytes, id), v = X.file_bytes(id) / 2^20; else, v = NaN; end
end

function v = mean_or_nan(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = mean(x); end
end

% =========================================================================
% 4 mesh by depth

function [out, R] = sec_mesh(X, out)
out{end+1} = ''; out{end+1} = '4. RETURNED MESH BY CALL-CHAIN DEPTH (depth 0 = cold / refinement)';
A = X.A; n = numel(A);
depth = nan(1, n);
for k = 1:n
    depth(k) = depth_of(X, k, depth);
end
rows = {};
for dd = unique(depth(isfinite(depth)))
    m = depth == dd;
    v = m & [A.valid]; iv = m & ~[A.valid] & arrayfun(@returned, A);
    mv = [A(v).returned_mesh]; mi = [A(iv).returned_mesh];
    out{end+1} = sprintf(['  depth %2d: valid %4d mesh median %6.0f max %6.0f | invalid %4d mesh median %6.0f ' ...
        'max %6.0f | guess nodes median %6.0f max %6.0f'], dd, sum(v), med(mv), max_or_nan(mv), sum(iv), ...
        med(mi), max_or_nan(mi), med(X.gn(m)), max_or_nan(X.gn(m))); %#ok<AGROW>
    rows{end+1} = {dd, sum(v), med(mv), max_or_nan(mv), sum(iv), med(mi), max_or_nan(mi), med(X.gn(m)), ...
        max_or_nan(X.gn(m))}; %#ok<AGROW>
end
write_csv(X, 'mesh_by_depth.csv', {'depth', 'n_valid', 'mesh_median_valid', 'mesh_max_valid', ...
    'n_invalid', 'mesh_median_invalid', 'mesh_max_invalid', 'guess_median', 'guess_max'}, rows);
R = struct('max_depth', max_or_nan(depth));
end

function dd = depth_of(X, k, depth)
dd = 0; a = X.A(k); seen = 0;
while ~isempty(a.parent_attempt) && isKey(X.ids, a.parent_attempt) && ~strcmp(a.kind, 'cold') && ...
        ~strcmp(a.kind, 'refinement')
    j = X.ids(a.parent_attempt);
    if isfinite(depth(j)), dd = dd + 1 + depth(j); return; end
    dd = dd + 1; a = X.A(j); seen = seen + 1;
    if seen > numel(X.A), break; end
end
end

function v = max_or_nan(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = max(x); end
end

% =========================================================================
% 5 route yield

function [out, R] = sec_yield(X, out)
out{end+1} = ''; out{end+1} = '5. ROUTE YIELD';
A = X.A;
rows = {};
nodes = unique({X.roots.key});
rts = unique([X.roots.routes]);
R = struct();
for i = 1:numel(rts)
    nr = numel(unique({X.roots(cellfun(@(r) any(strcmp(r, rts{i})), {X.roots.routes})).key}));
    calls = sum(strcmp({A.route}, rts{i}) | (strcmp(rts{i}, 'cold') & strcmp({A.kind}, 'cold')));
    out{end+1} = sprintf('  route %-16s nodes reached %4d  calls %5d  calls per reached node %.3g', ...
        rts{i}, nr, calls, calls / max(nr, 1)); %#ok<AGROW>
    rows{end+1} = {['route ' rts{i}], nr, calls, calls / max(nr, 1)}; %#ok<AGROW>
end
out{end+1} = sprintf('  nodes holding >= 1 root: %d', numel(nodes));
R.nodes_reached = numel(nodes);
if ~isempty(X.Q) && isfield(X.Q, 'edges') && ~isempty(X.Q.edges)
    E = X.Q.edges;
    for kd = unique({E.kind})
        e = E(strcmp({E.kind}, kd{1}));
        rf = e(strcmp({e.status}, 'refused'));
        cp = e(strcmp({e.status}, 'completed'));
        out{end+1} = sprintf('  edges %-16s %4d: completed %4d (calls %5d), refused %4d (calls %5d, %.3g per refused edge)', ...
            kd{1}, numel(e), numel(cp), sum([cp.calls]), numel(rf), sum([rf.calls]), ...
            sum([rf.calls]) / max(1, numel(rf))); %#ok<AGROW>
        rows{end+1} = {['edges ' kd{1}], numel(cp), sum([cp.calls]), sum([rf.calls]) / max(1, numel(rf))}; %#ok<AGROW>
    end
    % rescue by region and rung
    ring = X.Q.ring;
    for reg = {'interior', 'ring'}
        for kd = {'rescue_diagonal', 'rescue_secant', 'rescue_deep'}
            e = E(strcmp({E.kind}, kd{1}) & (ring([E.q]).' == strcmp(reg{1}, 'ring')));
            if isempty(e), continue; end
            nodes_y = numel(unique([e(strcmp({e.status}, 'completed')).q]));
            out{end+1} = sprintf('  rescue %-8s %-16s edges %4d calls %5d nodes rescued %3d', reg{1}, kd{1}, ...
                numel(e), sum([e.calls]), nodes_y); %#ok<AGROW>
            rows{end+1} = {sprintf('rescue %s %s', reg{1}, kd{1}), nodes_y, sum([e.calls]), numel(e)}; %#ok<AGROW>
        end
    end
end
% The per-root pass (route root_pass) and the rescue no_basis count
R.root_pass = struct('enabled', X.S.root_pass.enabled, 'targets', 0, 'calls', 0, 'new_roots', 0, ...
    'matched', 0, 'refused', 0, 'calls_per_new_root', NaN);
rp_rows = {};
if ~isempty(X.Q) && isfield(X.Q, 'root_pass')
    RP = X.Q.root_pass;
    E = X.Q.edges;
    if isempty(E), pe = E; else, pe = E(strcmp({E.kind}, 'root_pass')); end
    R.root_pass.targets = sum([RP.targets]);
    R.root_pass.calls = sum([pe.calls]);
    R.root_pass.new_roots = sum([RP.new_roots]);
    R.root_pass.matched = sum([RP.matched]);
    R.root_pass.refused = sum(strcmp({pe.status}, 'refused'));
    R.root_pass.calls_per_new_root = R.root_pass.calls / max(1, R.root_pass.new_roots);
    out{end+1} = sprintf(['  per-root pass (root_pass.enabled %d): targets posed %d, calls %d, new roots %d, matched %d, ' ...
        'refused edges %d, calls per new root %.3g'], R.root_pass.enabled, R.root_pass.targets, R.root_pass.calls, ...
        R.root_pass.new_roots, R.root_pass.matched, R.root_pass.refused, R.root_pass.calls_per_new_root);
    rows{end+1} = {'root_pass', R.root_pass.new_roots, R.root_pass.calls, R.root_pass.calls_per_new_root};
    if R.root_pass.enabled
        out{end+1} = sprintf(['  NOTE: the per-root pass is on; it posed %d target(s) in this run''s region. A cross ' ...
            'pilot poses few or none (pass targets need a root next to a node of the region), so the pilot ' ...
            'underprices the pass in production; pass edges into already-reached nodes enter a; the cap ' ...
            'column omits the pass.'], R.root_pass.targets);
    end
    for k = 1:numel(RP)
        r = RP(k);
        rp_rows{end+1} = {r.search_pass, r.phase, r.round, r.start_roots, r.candidates, r.known_a, ...
            r.known_b_only, r.targets, r.edges, r.calls, r.completed, r.new_roots, ...
            r.matched, r.refused, r.unattempted, r.max_calls_edge, r.status}; %#ok<AGROW>
    end
end
write_csv(X, 'root_pass.csv', {'search_pass', 'phase', 'round', 'start_roots', 'candidates', 'known_a', ...
    'known_b_only', 'targets', 'edges', 'calls', 'completed', 'new_roots', 'matched', ...
    'refused', 'unattempted', 'max_calls_edge', 'status'}, rp_rows);
R.no_basis = 0;
if ~isempty(X.Q) && isfield(X.Q, 'rescue_result')
    R.no_basis = sum(strcmp(X.Q.rescue_result, 'no_basis') & ~X.Q.rescue_done);
    out{end+1} = sprintf(['  rescue rungs {%s}, ring %d; nodes with a zero-call ladder (no_basis: no lineage ' ...
        'pair at q-e, q-2e): %d'], strjoin(X.S.rescue.rungs, ','), X.S.rescue.ring, R.no_basis);
    rows{end+1} = {'rescue_no_basis', R.no_basis, 0, NaN};
end
% repeated calls (guess fingerprint). The retry of an
% interrupted record (retry_of set; every earlier call with its fingerprint
% is an interrupted record) is the same unit re-run after a kill, not a
% repeat; it is counted separately.
[R.repeated, R.retries, R.fingerprinted] = repeats(X);
out{end+1} = sprintf(['  repeated calls (equal guess fingerprint to an earlier call; the retry of an interrupted ' ...
    'record excluded): %d of %d fingerprinted calls; retries of interrupted records (retry_of set): %d'], ...
    R.repeated, R.fingerprinted, R.retries);
rows{end+1} = {'repeated_calls', R.repeated, R.fingerprinted, NaN};
rows{end+1} = {'retries_of_interrupted', R.retries, R.fingerprinted, NaN};
write_csv(X, 'route_yield.csv', {'item', 'nodes_or_completed', 'calls', 'per_unit'}, rows);
end

function [nrep, nretry, nfp] = repeats(X)
% A fingerprinted call whose fingerprint equals an earlier call's is a
% repeat, unless it is the retry of an interrupted record: retry_of set and
% every earlier call with that fingerprint is an interrupted record.
% A retry that was itself killed is an interrupted record without retry_of
% (sa_resume writes it blank); it is counted as a retry on the same rule.
nrep = 0; nretry = 0; nfp = 0;
seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');   % fp -> any non-interrupted earlier
hasr = isfield(X.A, 'retry_of');
for k = 1:numel(X.A)
    f = X.fp{k};
    if isempty(f), continue; end
    nfp = nfp + 1;
    intr = strcmp(X.A(k).outcome, 'interrupted');
    if isKey(seen, f)
        retry = ((hasr && ~isempty(X.A(k).retry_of)) || intr) && ~seen(f);
        if retry, nretry = nretry + 1; else, nrep = nrep + 1; end
        seen(f) = seen(f) || ~intr;
    else
        seen(f) = ~intr;
    end
end
end

% =========================================================================
% 6 reach map and frontier

function [out, R] = sec_reach(X, out)
out{end+1} = ''; out{end+1} = '6. REACH MAP (R reached, x refused, . unattempted, s policy-skipped, u unfinished; ring in lower case; blank = outside the region)';
R = struct('classes', struct(), 'classes_ring', struct(), 'n_region', 0, 'n_ring_region', 0);
for c = {'reached', 'refused', 'unattempted', 'policy_skipped', 'unfinished'}
    R.classes.(c{1}) = 0; R.classes_ring.(c{1}) = 0;
end
hr = {'key', 'index', 'coordinate', 'ring', 'in_region', 'class', 'n_roots'};
hf = {'edge', 'kind', 'p', 'q', 'ring', 'fraction_valid', 'last_valid_coordinate', 'failing_coordinate', ...
    'failing_check'};
if isempty(X.Q)
    % 0 calls (e.g. disk_floor at the start): the tables are written with
    % their headers
    out{end+1} = '  (no scheduler state)';
    write_csv(X, 'reach_map.csv', hr, {});
    write_csv(X, 'frontier.csv', hf, {});
    return
end
Q = X.Q; L = X.L; E = Q.edges;
n = Q.n;
cls = repmat({''}, n, 1);
for q = 1:n
    if ~Q.box(q), continue; end
    key = L.key{q};
    if isKey(X.nroots, key), cls{q} = 'reached'; continue; end
    if isempty(E), into = []; else, into = E([E.q] == q & [E.calls] > 0); end
    if ~isempty(into) && any(strcmp({into.status}, 'unfinished')), cls{q} = 'unfinished';
    elseif ~isempty(into), cls{q} = 'refused';
    elseif X.S.search.ring_last && Q.ring(q) && ~any(strcmp(Q.phases_run, 'ring')), cls{q} = 'policy_skipped';
    else, cls{q} = 'unattempted';
    end
end
idx = L.index;
if size(idx, 2) == 2
    k1 = unique(idx(:, 1)); k2 = unique(idx(:, 2));
    out{end+1} = sprintf('  axis 1 (rows) %d..%d, axis 2 (columns) %d..%d', k1(1), k1(end), k2(1), k2(end));
    for a = k1(end:-1:1).'
        line = blanks(numel(k2));
        for b = 1:numel(k2)
            q = find(idx(:, 1) == a & idx(:, 2) == k2(b), 1);
            c = sym_of(cls{q});
            if Q.ring(q), c = lower(c); end
            line(b) = c;
        end
        out{end+1} = sprintf('  %4d %s', a, line); %#ok<AGROW>
    end
end
rows = cell(1, n);
for q = 1:n
    rows{q} = {L.key{q}, mat2str(idx(q, :)), mat2str(L.coordinate(q, :), 10), Q.ring(q), Q.box(q), cls{q}, ...
        nroots_of(X, L.key{q})};
end
write_csv(X, 'reach_map.csv', hr, rows);
for c = {'reached', 'refused', 'unattempted', 'policy_skipped', 'unfinished'}
    R.classes.(c{1}) = sum(strcmp(cls, c{1}));
    R.classes_ring.(c{1}) = sum(strcmp(cls, c{1}) & Q.ring(:));
end
R.n_region = sum(Q.box); R.n_ring_region = sum(Q.box & Q.ring);
out{end+1} = sprintf('  region nodes %d (ring %d): reached %d, refused %d, unattempted %d, policy-skipped %d, unfinished %d', ...
    R.n_region, R.n_ring_region, R.classes.reached, R.classes.refused, R.classes.unattempted, ...
    R.classes.policy_skipped, R.classes.unfinished);
% frontier table
out{end+1} = '  Frontier (refused edges): last valid coordinate -> failing coordinate, failing check';
rows = {};
if ~isempty(E)
    for e = E(strcmp({E.status}, 'refused'))
        kp = Q.index(e.p, :); kq = Q.index(e.q, :);
        lv = L.anchor + sa_lattice('steps', L) .* (kp + e.fraction_valid * (kq - kp));
        ff = NaN(size(lv));
        if ~isempty(e.attempts) && isKey(X.ids, e.attempts{end})
            a = X.A(X.ids(e.attempts{end}));
            ff = a.coordinate;
            chk = strsplit(a.failing, ';'); chk = chk{1};
            if isempty(chk), chk = a.outcome; end
        else
            chk = e.failing;
        end
        out{end+1} = sprintf('    edge %4d %-10s %s -> %s: valid to %s, failed at %s on %s', e.id, e.kind, ...
            L.key{e.p}, L.key{e.q}, mat2str(lv, 4), mat2str(ff, 4), chk); %#ok<AGROW>
        rows{end+1} = {e.id, e.kind, L.key{e.p}, L.key{e.q}, Q.ring(e.q), e.fraction_valid, mat2str(lv, 10), ...
            mat2str(ff, 10), chk}; %#ok<AGROW>
    end
end
write_csv(X, 'frontier.csv', hf, rows);
end

function n = nroots_of(X, key)
if isKey(X.nroots, key), n = X.nroots(key); else, n = 0; end
end

% =========================================================================
% 7 identity and edge_large_change

function [out, R, LC] = sec_identity(X, out)
out{end+1} = ''; out{end+1} = '7. IDENTITY';
R = struct();
P = X.pairs;
band = X.S.identity.near_threshold;
nb = 0; nd = 0; np = 0;
if numel(P) > 1
    h = P{1};
    ck = find(strcmp(h, 'kind')); cr = find(strcmp(h, 'scaled_rms_v1'));
    ca = find(strcmp(h, 'root_a')); cb = find(strcmp(h, 'root_b'));
    for k = 2:numel(P)
        r = P{k};
        v = str2double(r{cr});
        np = np + 1;
        nb = nb + (v >= band(1) && v < band(2));
        nd = nd + (strcmp(r{ck}, 'same_target') && ~strcmp(r{ca}, r{cb}));
    end
end
multi = sum(cell2mat(values(X.nroots)) >= 2);
out{end+1} = sprintf(['  identity pairs %d; near-threshold band [%g, %g): %d; same-target distinct-root pairs %d; ' ...
    'nodes with >= 2 roots %d'], np, band(1), band(2), nb, nd, multi);
% Multiplicity is tested only with all_seedable cold targets (per_seed:
% not by cold; only continuation and the pass look; rescue fills empty nodes)
note = sa_report('multiplicity_note', X.S);
if ~isempty(note), out{end+1} = ['  ' note]; end
R = struct('pairs', np, 'band', nb, 'distinct_pairs', nd, 'multi_root_nodes', multi, ...
    'multiplicity_tested', isempty(note));
LC = large_change(X);
out = [out, LC.lines];
write_csv(X, 'identity.csv', {'item', 'value'}, {{'pairs', np}, {'near_threshold_band', nb}, ...
    {'same_target_distinct_pairs', nd}, {'nodes_with_2plus_roots', multi}, ...
    {'multiplicity_tested', double(isempty(note))}});
end

function LC = large_change(X)
% edge_large_change, zero solve.
LC = struct('lines', {{}}, 'median', NaN, 'n_basis', 0, 'threshold', NaN, 'status', 'not_computed', ...
    'edges', []);
hl = {'edge', 'kind', 'phase', 'calls', 'd_parent_child', 'max_dsub_over_ssub', 'd', 'edge_large_change', ...
    'attempt_shape_switch', 'second_difference'};
if isempty(X.Q) || ~isfield(X.Q, 'edges') || isempty(X.Q.edges)
    LC.lines{end+1} = '  edge_large_change: not_computed (no edges)';
    write_csv(X, 'large_change.csv', hl, {});   % every table at every stop
    return
end
E = X.Q.edges;
M = sa_model(X.S.model);
K = X.E.store_settings;
fs = K.field_scale(:);
rows_id = [M.identity_fields.row];
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
prof = @(id) profile_of(X, id, K.comparison_mesh, rows_id, cache);
comp = find(strcmp({E.status}, 'completed'));
T = struct('id', {}, 'kind', {}, 'phase', {}, 'calls', {}, 'd', {}, 'd_pc', {}, 'd_sub_norm', {}, ...
    'shape_switch', {}, 'second_diff', {}, 'dmax_check', {});
for k = comp
    e = E(k);
    if isfield(e, 'd_sub') && ~isempty(e.d_sub)
        dn = max(e.d_sub ./ e.s_sub); dchk = max(e.d_sub);
    else
        [dsub, ssub] = substeps_of(X, e, prof, fs);
        if isempty(dsub), dn = NaN; dchk = NaN; else, dn = max(dsub ./ ssub); dchk = max(dsub); end
    end
    d = e.d_parent_child;
    if e.calls > 1 && isfinite(dn), d = max(d, dn); end
    sw = NaN;
    if isKey(X.ids, e.parent_attempt) && ~isempty(e.attempts) && isKey(X.ids, e.attempts{end})
        sw = ~strcmp(X.A(X.ids(e.parent_attempt)).observed_family, X.A(X.ids(e.attempts{end})).observed_family);
    end
    T(end+1) = struct('id', e.id, 'kind', e.kind, 'phase', e.phase, 'calls', e.calls, 'd', d, ...
        'd_pc', e.d_parent_child, 'd_sub_norm', dn, 'shape_switch', sw, ...
        'second_diff', second_diff(X, E, e, prof, fs), 'dmax_check', abs(dchk - e.d_max_substep)); %#ok<AGROW>
end
basis = T(strcmp({T.kind}, 'continuation') & strcmp({T.phase}, 'interior') & [T.calls] == 1);
LC.n_basis = numel(basis);
cint = T(strcmp({T.kind}, 'continuation') & strcmp({T.phase}, 'interior'));
if numel(basis) >= 30
    LC.median = median([basis.d_pc]);
    LC.threshold = X.S.report.large_change_factor * LC.median;
    LC.status = 'computed';
    fl = [T.d] > LC.threshold;
    for k = 1:numel(T), T(k).flag = fl(k); end
    fli = [cint.d] > LC.threshold;
    LC.lines{end+1} = sprintf(['  edge_large_change: d > %g x median %.4g = %.4g (basis: %d completed 1-call ' ...
        'interior continuation edges): flagged %d / %d completed interior continuation edges (%d of them ' ...
        'with an attempt-level shape switch; switches in total %d); all completed edges flagged %d / %d'], ...
        X.S.report.large_change_factor, LC.median, LC.threshold, numel(basis), sum(fli), numel(cint), ...
        sum(fli & [cint.shape_switch] == 1), sum([cint.shape_switch] == 1), sum(fl), numel(T));
    LC.flagged_ids = [T(fl).id];
    LC.flagged_interior = sum(fli); LC.n_interior = numel(cint);
    LC.flag_switch = sum(fli & [cint.shape_switch] == 1);
else
    for k = 1:numel(T), T(k).flag = NaN; end
    LC.lines{end+1} = sprintf('  edge_large_change: not_computed (%d completed 1-call interior continuation edges < 30)', ...
        numel(basis));
    LC.flagged_ids = [];
end
sd = [T.second_diff]; sd = sd(isfinite(sd));
if ~isempty(sd)
    qs = quantile_set(sd);
    LC.lines{end+1} = sprintf(['  second difference d(child, 2 parent - grandparent)/d(parent, child) on %d collinear ' ...
        'chains: min %.3g q25 %.3g median %.3g q75 %.3g max %.3g'], numel(sd), qs);
end
LC.dmax_mismatch = sum([T.dmax_check] > 0);
LC.edges = T;
rows = arrayfun(@(t) {t.id, t.kind, t.phase, t.calls, t.d_pc, t.d_sub_norm, t.d, t.flag, t.shape_switch, ...
    t.second_diff}, T, 'UniformOutput', false);
write_csv(X, 'large_change.csv', hl, rows);
end

function q = quantile_set(x)
x = sort(x); n = numel(x);
pick = @(p) x(max(1, min(n, round(p * (n - 1)) + 1)));
q = [x(1) pick(0.25) median(x) pick(0.75) x(end)];
end

function P = profile_of(X, id, cm, rows, cache)
if isKey(cache, id), P = cache(id); return; end
Z = load(fullfile(X.dir, 'attempts', [id '.mat']), 'sol');
P = deval(Z.sol, cm, rows);
cache(id) = P; %#ok<NASGU>
end

function d = dist(a, b, fs)
dd = (a - b) ./ fs;
d = sqrt(sum(dd(:).^2) / numel(dd));
end

function [dsub, ssub] = substeps_of(X, e, prof, fs)
% the valid chain parent -> valid substeps -> child (recomputed, older folders)
dsub = []; ssub = [];
if ~isKey(X.ids, e.parent_attempt), return; end
ids = e.attempts;
if strcmp(e.rule, 'myxo_restart')
    k0 = find(~arrayfun(@(i) X.A(X.ids(ids{i})).valid, 1:numel(ids)), 1, 'last');
    if isempty(k0), k0 = 0; end
    ids = ids(k0+1:end);
end
prev = prof(e.parent_attempt); fprev = 0;
if isfield(e, 'start_fraction') && e.start_fraction > 0, fprev = e.start_fraction; end
for i = 1:numel(ids)
    a = X.A(X.ids(ids{i}));
    if ~a.valid, continue; end
    cur = prof(ids{i});
    dsub(end+1) = dist(prev, cur, fs); %#ok<AGROW>
    ssub(end+1) = a.fraction - fprev; %#ok<AGROW>
    prev = cur; fprev = a.fraction;
end
end

function r = second_diff(X, E, e, prof, fs)
r = NaN;
if e.lineage <= 0 || isempty(e.attempts), return; end
kp = X.Q.index(e.p, :); kq = X.Q.index(e.q, :);
g = row_of(X.Q, 2 * kp - kq);
if g == 0, return; end
m = find([E.lineage] == e.lineage & [E.q] == e.p & [E.p] == g & strcmp({E.status}, 'completed'), 1);
if isempty(m), return; end
eg = E(m);
if ~strcmp(eg.attempts{end}, e.parent_attempt), return; end
Pg = prof(eg.parent_attempt); Pp = prof(e.parent_attempt); Pq = prof(e.attempts{end});
den = dist(Pp, Pq, fs);
if den > 0, r = dist(Pq, 2 * Pp - Pg, fs) / den; end
end

function r = row_of(Q, idx)
if any(idx < Q.kmin) || any(idx > Q.kmax), r = 0; return; end
r = 1 + sum((idx - Q.kmin) .* Q.strides);
end

% =========================================================================
% 8 production projection

function [out, P] = sec_projection(X, R, out)
out{end+1} = ''; out{end+1} = '8. PRODUCTION PROJECTION (central / high / cap)';
S = X.S; P = struct();
cen = X.E.census;
if isfield(cen, 'computed') && ~cen.computed
    % focus: a focus run computes no census (no library)
    out{end+1} = '  not computed: a focus run has no library census and no cold targets';
    return
end
fam = X.E.targets.families;
fi = find(ismember(cen.names, fam));
if isfield(cen, 'seedable')
    cold_all = 0;
    for f = fi(:).'
        keep = sa_targets('thin_keep', S, cen.names{f}, X.L.index);   % thinned count (an input)
        cold_all = cold_all + sum(cen.seedable(:, f) & keep(:));
    end
    c9 = sa_modes('central_nine', X.L);
    cold_c9 = sum(sum(cen.seedable(c9, fi)));
else
    cold_all = sum(cen.counts(fi)); cold_c9 = NaN;
end
% The cold count of the projected production run: per_seed (the production
% default) when the census holds the seed counts and the
% folder's own policy is not all_seedable; else all_seedable (thinned)
cold = cold_all; cold_pol = 'all_seedable';
if isfield(cen, 'seed_counts') && ~strcmp(S.targets.policy, 'all_seedable')
    cold = sum(cen.seed_counts(fi)); cold_pol = 'per_seed';
end
N = numel(X.L.key);
rc = R.reach.classes;
nreg = getf(R.reach, 'n_region', 0);
% The reach fraction over the WHOLE planned region; refused
% and never-attempted region nodes are both unreached and both charged b
f = NaN; if nreg > 0, f = rc.reached / nreg; end
unreached = nreg - rc.reached;
completed = strcmp(X.stop_reason, 'completed');
E = []; if ~isempty(X.Q) && isfield(X.Q, 'edges'), E = X.Q.edges; end
a = NaN; b = NaN; nref = 0;
if ~isempty(E)
    reached_q = unique([E(strcmp({E.status}, 'completed')).q]);
    if ~isempty(reached_q)
        a = sum([E(ismember([E.q], reached_q)).calls]) / numel(reached_q);
    end
    refused_q = setdiff(unique([E([E.calls] > 0).q]), reached_q);
    nref = numel(refused_q);
    if nref > 0, b = sum([E(ismember([E.q], refused_q)).calls]) / nref; end
end
h = S.continuation.max_halvings;
% No refused node measured yet (typically a partial pilot):
% b is the per-node bound of one refused edge, (2h + 1) + rescue calls
b_bound = (2 * h + 1) + S.rescue.max_calls_per_node * S.routes.rescue;
b_src = 'measured';
if nref == 0, b = b_bound; b_src = 'bound: no refused node measured, (2h + 1) + rescue.max_calls_per_node'; end
if S.routes.rescue && isequal(S.rescue.rungs, {'secant'}) && nref == 0
    b_src = [b_src sprintf(['; conservative under secant-only rescue (its per-node bound is 2 x n_pools = %d ' ...
        'calls, not max_calls_per_node = %d)'], 2 * size(X.L.index, 2), S.rescue.max_calls_per_node)];
end
per_node_cap = (S.continuation.max_lineages + 1) * S.continuation.max_edges_per_node * (2 * h + 1) + ...
    S.rescue.max_calls_per_node;
% with the lineage slots at or above the number of box nodes
% the per-node cap bounds nothing: the cap column is unbounded (Inf)
nbox = nreg; if nbox == 0, nbox = N; end
unbounded = S.continuation.max_lineages >= nbox;
C = [cold + N * f * a + N * (1 - f) * b, cold + N * a + N * (1 - f) * b, cold + N * per_node_cap];
if unbounded, C(3) = Inf; end
if ~completed && nreg == 0
    out{end+1} = sprintf(['  PARTIAL: stop reason %s before the scheduler held a region (cold pass unfinished or 0 ' ...
        'calls): f, a and the central and high calls are not computed; resume the pilot to completion'], X.stop_reason);
elseif ~completed
    out{end+1} = sprintf(['  PARTIAL: stop reason %s before the search finished. The %d region nodes not reached so far ' ...
        '(incl. nodes not yet attempted) count as unreached and are charged b; this is not the pilot''s central ' ...
        'estimate: resume the pilot to completion for that'], X.stop_reason, unreached);
end
if strcmp(cold_pol, 'per_seed')
    cold_txt = sprintf('cold targets (per_seed, one per library seed) %d [all_seedable %d; central_nine %d]', ...
        cold, cold_all, cold_c9);
else
    cold_txt = sprintf('cold targets (all_seedable) %d [central_nine %d]', cold_all, cold_c9);
end
out{end+1} = sprintf(['  calls = cold + N f a + N (1 - f) b: %s; N = %d ' ...
    'admissible nodes; f = %.3g = reached %d / %d region nodes (refused %d, not attempted %d, policy-skipped %d, ' ...
    'unfinished %d: every unreached node is charged b); a = %.3g calls per reached node; b = %.3g calls per ' ...
    'refused node (%s)'], cold_txt, N, f, rc.reached, nreg, rc.refused, rc.unattempted, ...
    rc.policy_skipped, rc.unfinished, a, b, b_src);
if unbounded
    out{end+1} = sprintf(['  central %.0f; high (every node reached: cold + N a + N (1 - f) b) %.0f; cap unbounded ' ...
        '(lineage caps above the box size: max_lineages %g >= %d box nodes)'], C(1), C(2), ...
        S.continuation.max_lineages, nbox);
else
    out{end+1} = sprintf(['  central %.0f; high (every node reached: cold + N a + N (1 - f) b) %.0f; cap (cold + N x ' ...
        '[(max_lineages + 1) x max_edges_per_node x (2h + 1) + rescue.max_calls_per_node = %d]) %.0f'], C(1), C(2), ...
        per_node_cap, C(3));
end
out{end+1} = sprintf(['  with central_nine cold instead: central %.0f, high %.0f (branch connected to the central nine ' ...
    'only; multiplicity not tested: coexisting roots elsewhere are not looked for)'], C(1) - cold + cold_c9, ...
    C(2) - cold + cold_c9);
% hours: solver s per call by the NMax share x the live process/solver
% factor + replay at report.replay_s_per_record
sh = R.cost.nmax_share;
if sh == 0, s_call = R.cost.s_valid; else, s_call = sh * R.cost.s_nmax + (1 - sh) * R.cost.s_valid; end
if ~isfinite(s_call), s_call = R.cost.corrector_per_call; end
ovh = R.cost.overhead; ovh_note = '';
if ~isfinite(ovh), ovh = 1; ovh_note = ' (not measured: 1 assumed)';
elseif ovh < 1, ovh = 1; ovh_note = ' (below 1: 1 used)';
end
wall = X.seg.budget.wall_s;
H_live = C * s_call * ovh / 3600;
nres = 0; if isfinite(wall) && isfinite(H_live(1)), nres = floor(H_live(1) * 3600 / wall); end
Hh = H_live + nres * X.rate * C / 2 / 3600;
out{end+1} = sprintf(['  MATLAB-h = calls x (NMax share %.3g x %.3g s + (1 - share) x %.3g s) x process/corrector %.3g%s ' ...
    '+ replay %.3g s x calls/2 x resumes (%d at budget.wall_s %g): central %.2f h, high %.2f h, cap %.2f h'], ...
    sh, R.cost.s_nmax, R.cost.s_valid, ovh, ovh_note, X.rate, nres, wall, Hh);
seg_T = [8 16]; seg_n = floor(H_live(1) ./ seg_T); seg_add = seg_n * X.rate * C(1) / 2 / 3600;
out{end+1} = sprintf(['  replay at %.3g s per stored record (report.replay_s_per_record, %s): in segments of 8 h the ' ...
    'central needs %d resume(s), +%.2f h; of 16 h, %d, +%.2f h'], X.rate, X.rate_src, seg_n(1), seg_add(1), ...
    seg_n(2), seg_add(2));
mbc = R.cost.nmax_share * R.cost.mb_nmax + (1 - R.cost.nmax_share) * R.cost.mb_valid;
if ~isfinite(mbc), mbc = R.cost.mb_valid; end
rpn = numel(X.roots) / max(1, R.yield.nodes_reached);
roots_c = N * rpn * [f 1 1];            % the same f as the calls
Dg = (C * mbc + roots_c * R.cost.mb_per_root) / 1024;
free = sa_search('free_gib', X.dir) - X.seg.checkpoint.min_free_gib;
out{end+1} = sprintf(['  disk GB = calls x MB per call (%.3g) + roots (N f x %.3g roots per reached node = %.0f; high and cap %.0f) ' ...
    'x MB per root (%.3g): central %.2f, high %.2f, cap %.2f; free above the %g GiB floor now %.2f GiB%s'], mbc, rpn, ...
    roots_c(1), roots_c(2), R.cost.mb_per_root, Dg, X.seg.checkpoint.min_free_gib, free, warn(Dg(1), free));
P = struct('cold_all', cold_all, 'cold', cold, 'cold_policy', cold_pol, 'cold_c9', cold_c9, 'N', N, 'f', f, 'a', a, 'b', b, 'calls', C, ...
    'hours', Hh, 'disk_gb', Dg, 'free_gib', free, 'n_region', nreg, 'reached', rc.reached, ...
    'unreached', unreached, 'n_refused_measured', nref, 'b_source', b_src, 'b_bound', b_bound, ...
    'partial', ~completed, 's_call', s_call, 'overhead', ovh, 'resumes', nres, 'replay_rate', X.rate, ...
    'replay_rate_source', X.rate_src, 'hours_live', H_live, 'segment_h', seg_T, 'segment_resumes', seg_n, ...
    'segment_replay_h', seg_add, 'roots', roots_c);
write_csv(X, 'projection.csv', {'quantity', 'central', 'high', 'cap'}, {{'calls', C(1), C(2), C(3)}, ...
    {'matlab_h', Hh(1), Hh(2), Hh(3)}, {'disk_gb', Dg(1), Dg(2), Dg(3)}, ...
    {'reach_fraction_f', f, 1, NaN}, {'calls_per_reached_node_a', a, a, NaN}, ...
    {'calls_per_refused_node_b', b, b, NaN}, {'partial', double(~completed), NaN, NaN}});
end

function t = warn(dg, free)
t = '';
if dg > 0.8 * free, t = '  WARNING: central disk above 80 % of the free space above the floor'; end
end

function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function c = sym_of(cl)
switch cl
    case 'reached', c = 'R';
    case 'refused', c = 'x';
    case 'unattempted', c = '.';
    case 'policy_skipped', c = 's';
    case 'unfinished', c = 'u';
    otherwise, c = ' ';
end
end
