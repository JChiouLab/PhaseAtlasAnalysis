function varargout = sa_report(action, varargin)
%SA_REPORT  Per-call records and the run's output files.
%
%   A = SA_REPORT('blank')                    one empty attempt record
%   r = SA_REPORT('blank_refinement_row')     one empty refinement row
%   R = SA_REPORT('blank_refinement')         1x0 refinement table
%   P = SA_REPORT('pairs', st)                identity pairs
%   SA_REPORT('write', ctx)                   attempts.mat/.csv, queue.mat,
%                                             summary.txt (atomic), stamped with
%                                             the attempt count; the root store
%                                             lives in per-node shards roots/
%                                             written by the executor (= the
%                                             two actions below, full render)
%   SA_REPORT('write_state', ctx)             the checkpoint files only:
%                                             attempts.mat, process_time.mat,
%                                             queue.mat (every
%                                             checkpoint.every_calls and at
%                                             the stop; a resume compares
%                                             them with its replay)
%   tf = SA_REPORT('due', ctx)                the human-facing files are due
%                                             (report.files_interval_s seconds
%                                             since ctx.t_reports; checked by
%                                             the executor after each live call)
%   [rep, t] = SA_REPORT('write_reports', ctx [, check])
%                                             the human-facing files only:
%                                             attempts.csv, summary.txt,
%                                             roots_index.csv, rendered with
%                                             the cache ctx.rep (returned
%                                             updated). check = true (finish):
%                                             the cached render is compared
%                                             with a full render from
%                                             scratch; on a difference the
%                                             full one is written and the
%                                             journal and the screen say so;
%                                             t: the next interval's timer
%   [T, rep] = SA_REPORT('render', ctx, rep, wall)
%                                             the three texts (T.attempts_csv,
%                                             T.summary, T.roots_index); rep =
%                                             [] renders from scratch; wall is
%                                             the summary's wall seconds
%       Render cache (memory only, rebuilt from the state after a resume):
%       each attempt row's csv line, plain-words reason and failing-check
%       list (rows never change once appended; an 'interrupted' row is
%       rendered again every time, because its reason names the call that
%       re-ran it, which may come after a later resume), each node's
%       roots_index lines (rendered again when the node's store gained a
%       discovery: next_discovery_number), concatenated in lattice order
%       through the key -> node map, and the running maximum of the
%       per-root amplification (a root's representative never changes).
%       Each output line is the same function of the same state as a full
%       render.
%   SA_REPORT('final', ctx)                   identity_pairs.csv
%   txt = SA_REPORT('summary_text', ctx)
%   R = SA_REPORT('roots', ctx)               every root of every shard
%   C = SA_REPORT('route_counts', ctx)        route membership counts
%   SA_REPORT('pilot', run_dir)               pilot report (sa_pilot_report;
%                                             zero solves, pure function of the
%                                             run folder)
%   t = SA_REPORT('roots_index', ctx)         roots_index.csv text (written
%                                             with every 'write')
%   t = SA_REPORT('multiplicity_note', S)     '' when
%                                             S.targets.policy is
%                                             'all_seedable'; otherwise the
%                                             MULTIPLICITY NOT TESTED line
%                                             (summary.txt, pilot report);
%                                             all_seedable with
%                                             targets.thin -> MULTIPLICITY
%                                             PARTIALLY TESTED; per_seed ->
%                                             MULTIPLICITY NOT TESTED BY COLD
%   summary.txt under per_seed adds the seed count and the seedable
%   (family, node) pairs not called cold.
%   summary.txt: with a finite cold.mesh_cap, the line "cold calls cut
%   by the cap (NMax = mesh_reject = c): n of m"; the search block prints
%   the route defaults (rescue rungs, ring, root_pass.enabled), the
%   no_basis count and the per-root pass rounds.
%   summary.txt says "mesh proposals recorded
%   for n of m calls", and the cut-by-cap line says whether its count is
%   exact (every cold call's proposals recorded) or a lower bound.
%   r = SA_REPORT('plain_reason', A)          the plain-words reason of one
%                                             attempt row ('' for a stored
%                                             root); derived at write time,
%                                             never stored in a row
%   T = SA_REPORT('plain_table')              the check / outcome -> words table
%   SA_REPORT('refuse_if_compacted', run_dir, tool)
%       error SA:REPORT:E_COMPACTED ("compacted on <date>: n files removed")
%       when internal/compaction.mat exists: a tool that enumerates attempt
%       files or reads records by id refuses rather than miscounts
%   X = SA_REPORT('cutset', run_dir, c)
%   X = SA_REPORT('cutset', run_dir, c, kinds, o)
%       Exact cut recount, zero solves, read only on the run folder
%       (outputs go to o.out_dir, default run_dir). kinds default {'cold'}
%       (refinement / continuation / rescue / pass only when asked). run_dir
%       is the run folder (its state folder is resolved; o.out_dir defaults
%       to the state folder); refused on a compacted folder. Each
%       non-interrupted call of those kinds is classed by the exact rule
%       for the hypothetical bound pair (c, c), with P
%       its largest proposal, m its returned mesh and nmax0 its own NMax:
%         c <= nmax0: cut iff P > c or (returned and m >= c), else not_cut;
%         c > nmax0:  cut if P > c; else, when the call did not stop on a
%                     proposal above nmax0 (over_nmax false): cut iff
%                     returned and m >= c, else not_cut; otherwise
%                     unknown_above_nmax;
%         P not recorded (NaN): unknown_not_recorded.
%       Interrupted records are excluded from every count. Writes
%       cutset_<c>.csv (one row per call), cutset_<c>_summary.txt (the
%       one-line summary, the histogram of P, and the cap-loss curve at
%       c = 3,601 / 5,000 / 10,000: cut calls and, among them, stored
%       accepted returns, each marked exact or not). Calls without a trace
%       are also classed by the cap-cut prediction rule (certain: NMax-bound or
%       returned mesh >= c; possible: a no-return with wall_s above
%       o.cut_time_s, default 0), labelled as a lower-bound prediction,
%       never presented as exact. o.write (default true).
%
%   Outcomes: accepted (valid node return, entered the store), invalid
%   (returned, failed the boundary), no_return (nothing returned),
%   refined_valid (refinement call that passed), substep_valid /
%   substep_invalid (route records at non-lattice coordinates), interrupted.
%   Cold targets never reached are listed as unfinished.

switch action
    case 'blank',                varargout{1} = blank();
    case 'blank_refinement_row', varargout{1} = blank_ref();
    case 'blank_refinement',     varargout{1} = repmat(blank_ref(), 1, 0);
    case 'pairs',                varargout{1} = pairs(varargin{:});
    case 'write',                write_all(varargin{:});
    case 'write_state',          write_state(varargin{:});
    case 'write_reports',        [varargout{1}, varargout{2}] = write_reports(varargin{:});
    case 'due',                  varargout{1} = reports_due(varargin{:});
    case 'render',               [varargout{1}, varargout{2}] = render(varargin{:});
    case 'final',                write_pairs(varargin{:});
    case 'summary_text',         varargout{1} = summary_text_full(varargin{:});
    case 'roots',                varargout{1} = all_roots(varargin{:});
    case 'route_counts',         varargout{1} = route_counts(varargin{:});
    case 'pilot',                [varargout{1:nargout}] = sa_pilot_report(varargin{:});
    case 'roots_index',          varargout{1} = roots_index_full(varargin{:});
    case 'multiplicity_note',    varargout{1} = multiplicity_note(varargin{:});
    case 'cutset',               varargout{1} = cutset(varargin{:});
    case 'plain_reason',         varargout{1} = plain_reason(varargin{:});
    case 'plain_table',          varargout{1} = plain_table();
    case 'refuse_if_compacted',  refuse_if_compacted(varargin{:});
    otherwise, error('SA:REPORT', 'Unknown sa_report action %s.', action);
end
end

% =========================================================================

function A = blank()
A = struct('attempt_id', '', 'kind', '', 'parent_attempt', '', 'target_number', NaN, ...
    'family', '', 'target_key', '', 'index', [], 'coordinate', [], ...
    'seed_left', '', 'seed_right', '', 'outcome', '', 'failing', '', 'flags', '', ...
    'valid', false, 'corrector_class', '', 'corrector_reason', '', ...
    'corrector_nodal_min', NaN, 'convergence_class', '', 'wall_s', NaN, 'reltol', NaN, ...
    'returned_mesh', NaN, 'maxres', NaN, 'positivity_min', NaN, 'positivity_row', '', ...
    'positivity_xi', NaN, 'nodal_min', NaN, 'uniform_grid_min', NaN, ...
    'residual_scaled_sup', NaN, 'residual_scaled_rms', NaN, 'residual_unscaled_sup', NaN, ...
    'seam_concentration', NaN, 'seam_flux', NaN, 's_max', NaN, 'observed_family', '', ...
    'root_id', '', 'store_status', '', 'representative_distance', NaN, 'file', '', ...
    'route', '', 'phase', '', 'lineage', 0, 'edge_id', 0, 'fraction', NaN, ...
    'parent_root', '', 'guess_source', '', 'retry_of', '', 'unit_key', '', ...
    'positivity_sampled_min', NaN, 'residual_bvp4c', NaN, ...
    'nmax', NaN, 'cut_by_cap', false, ...                  % cold-cap fields
    'max_proposed_mesh', NaN);                             % NaN = not recorded
end

function r = blank_ref()
r = struct('family', '', 'parent_attempt', '', 'parent_root', '', 'target_key', '', ...
    'parent_maxres', NaN, 'refine_attempt', '', 'refined_maxres', NaN, 'status', '', ...
    'failing', '', 'scaled_rms', NaN, 'relsup_identity', NaN, 'relsup_all', NaN, ...
    'ratio', NaN, 'within_expectation', false);
end

% =========================================================================

function P = pairs(st)
P = struct('kind', {}, 'a', {}, 'b', {}, 'target_key', {}, 'family_a', {}, ...
    'family_b', {}, 'root_a', {}, 'root_b', {}, 'scaled_rms', {}, ...
    'relsup_identity', {}, 'relsup_all', {}, 'ratio', {}, 'route_a', {}, 'route_b', {});
% Same-target pairs over the accepted node returns of every route
acc = st.attempts(strcmp({st.attempts.outcome}, 'accepted'));
keys = unique({acc.target_key});
for k = 1:numel(keys)
    g = acc(strcmp({acc.target_key}, keys{k}));
    % each attempt file of the node is read once (the same load, so the
    % same sol, as one read per pair member)
    sol = cell(1, numel(g));
    if numel(g) >= 2
        for a = 1:numel(g)
            Z = load(fullfile(st.state_dir, g(a).file), 'sol');
            sol{a} = Z.sol;
        end
    end
    for a = 1:numel(g)
        for b = a+1:numel(g)
            m = sa_identity('pair', st.M, sol{a}, sol{b}, st.K, st.S);
            P(end+1) = struct('kind', 'same_target', 'a', g(a).attempt_id, ...
                'b', g(b).attempt_id, 'target_key', keys{k}, 'family_a', g(a).family, ...
                'family_b', g(b).family, 'root_a', g(a).root_id, 'root_b', g(b).root_id, ...
                'scaled_rms', m.scaled_rms, 'relsup_identity', m.relsup_identity, ...
                'relsup_all', m.relsup_all, 'ratio', m.ratio, 'route_a', g(a).route, ...
                'route_b', g(b).route); %#ok<AGROW>
        end
    end
end
for k = 1:numel(st.refinement)
    r = st.refinement(k);
    if ~strcmp(r.status, 'refined'), continue; end
    P(end+1) = struct('kind', 'refinement', 'a', r.parent_attempt, 'b', r.refine_attempt, ...
        'target_key', r.target_key, 'family_a', r.family, 'family_b', r.family, ...
        'root_a', r.parent_root, 'root_b', '', 'scaled_rms', r.scaled_rms, ...
        'relsup_identity', r.relsup_identity, 'relsup_all', r.relsup_all, ...
        'ratio', r.ratio, 'route_a', 'cold', 'route_b', 'refinement'); %#ok<AGROW>
end
end

% =========================================================================

function write_all(st)
write_state(st);
T = render(st, [], toc(st.t_run));
write_texts(st, T);
end

function write_state(st)
d = st.state_dir;                   % state files
c = numel(st.attempts);
gn = guess_nodes_of(st);
fp = guess_fp_of(st);
sa_write('mat', fullfile(d, 'attempts.mat'), struct('attempts', st.attempts, ...
    'refinement', st.refinement, 'calls', st.calls, 'stop_reason', st.stop_reason, ...
    'targets', st.T, 'attempt_count', c, 'guess_nodes', gn, 'guess_fingerprint', {fp}));
% The process time of this segment (volatile; its own file so the record
% files stay clock-free; the pilot report reads it)
sa_write('mat', fullfile(d, 'process_time.mat'), struct('attempt_count', c, 'process_s', toc(st.t_run)));
sa_write('mat', fullfile(d, 'queue.mat'), struct('attempt_count', c, 'Q', st.Q_at_last, ...
    'ring_calls', st.ring_calls_at_last));
end

function write_texts(st, T)
top = st.S.run_dir;                 % human-facing files
sa_write('text', fullfile(top, 'attempts.csv'), T.attempts_csv);
sa_write('text', fullfile(top, 'summary.txt'), T.summary);
sa_write('text', fullfile(top, 'roots_index.csv'), T.roots_index);
end

function tf = reports_due(st)
% the human-facing files are due: report.files_interval_s seconds since
% st.t_reports (the last write, or the start of the run); the interval is
% the resolved setting (its default lives in sa_settings only)
tf = toc(st.t_reports) >= st.S.report.files_interval_s;
end

function [rep, t_reports] = write_reports(st, check)
% t_reports: a new interval timer, started after the write
if nargin < 2, check = false; end
rep = [];
if isfield(st, 'rep'), rep = st.rep; end
wall = toc(st.t_run);                % one value for both renders
[T, rep] = render(st, rep, wall);
if check
    [F, rep_full] = render(st, [], wall);
    if ~isequal(T, F)
        names = {'attempts.csv', 'summary.txt', 'roots_index.csv'};
        f = fieldnames(T);
        bad = names(cellfun(@(x) ~isequal(T.(x), F.(x)), f.'));
        msg = sprintf(['report self-check: the cached render differs from the full render (%s); the full ' ...
            'render was written'], strjoin(bad, ', '));
        try
            sa_journal('note', st.state_dir, msg);
        catch
        end
        fprintf('WARNING %s\n', msg);
        T = F; rep = rep_full;
    end
end
write_texts(st, T);
t_reports = tic;
end

function t = summary_text_full(st)
T = render(st, [], toc(st.t_run));
t = T.summary;
end

function t = roots_index_full(st)
T = render(st, [], NaN);
t = T.roots_index;
end

% =========================================================================
% rendering (cached)

function rep = blank_cache(n_nodes)
rep = struct('n', 0, 'csv', {cell(1, 0)}, 'reason', {cell(1, 0)}, 'fl', {cell(1, 0)}, ...
    'vol', false(1, 0), 'n_nodes', n_nodes, 'ri', {repmat({cell(1, 0)}, 1, n_nodes)}, ...
    'ri_stamp', zeros(1, n_nodes), 'amp_done', zeros(1, n_nodes), 'amp_max', []);
end

function [T, rep] = render(st, rep, wall)
nn = numel(st.L.key);
A = st.attempts;
if isempty(rep) || ~isstruct(rep) || ~isfield(rep, 'n_nodes') || rep.n_nodes ~= nn || rep.n > numel(A)
    rep = blank_cache(nn);
end
gn = guess_nodes_of(st);
fp = guess_fp_of(st);
w = plain_words();
rt = retried_as(A);
[rep, lines] = render_rows(rep, A, gn, fp, w, rt);
T = struct();
T.attempts_csv = [strjoin([{csv_header()}, lines], newline) newline];
[rep, rlines, n_roots] = render_roots(st, rep, A);
T.summary = summary_text(st, wall, struct('reasons', {rep.reason}, 'fl', {[rep.fl{:}]}, ...
    'n_roots', n_roots, 'amp', rep.amp_max));
T.roots_index = [strjoin([{roots_header()}, rlines], newline) newline];
end

function [rep, lines] = render_rows(rep, A, gn, fp, w, rt)
% csv line, plain reason and failing-check list per row; rows 1..rep.n are
% cached (the interrupted ones rendered again), the new rows appended
cols = csv_cols();
has = isfield(A, cols);
n = numel(A);
n0 = rep.n;
rep.csv(n0+1:n) = {''}; rep.reason(n0+1:n) = {''}; rep.fl(n0+1:n) = {{}}; rep.vol(n0+1:n) = false;
redo = [find(rep.vol(1:n0)), n0+1:n];
for k = redo
    a = A(k);
    r = plain_reason(a, rt, w);
    rep.reason{k} = r;
    rep.csv{k} = csv_line(a, cols, has, gn(k), fp{k}, r);
    if k > n0
        rep.vol(k) = strcmp(a.outcome, 'interrupted');
        if ~isempty(a.failing) && ~strcmp(a.outcome, 'no_return')
            rep.fl{k} = strsplit(a.failing, ';');
        end
    end
end
rep.n = n;
lines = rep.csv(1:n);
end

function [rep, lines, total] = render_roots(st, rep, A)
% roots_index lines per node (lattice order) and the amplification maximum
keys = st.shards.keys;
pos = zeros(1, numel(keys));
for k = 1:numel(keys), pos(k) = node_pos(st, keys{k}); end
[pos, o] = sort(pos);
keys = keys(o);
total = 0;
fs = st.K.field_scale;
for k = 1:numel(keys)
    store = st.shards(keys{k});
    p = pos(k);
    if rep.ri_stamp(p) ~= store.next_discovery_number
        rep.ri{p} = root_lines(st, keys{k}, store, A, p);
        rep.ri_stamp(p) = store.next_discovery_number;
    end
    nr = numel(store.roots);
    for r = rep.amp_done(p)+1:nr
        a = sa_identity('amplification', store.roots(r).representative_profile, fs);
        if isempty(rep.amp_max), rep.amp_max = a; else, rep.amp_max = max(rep.amp_max, a); end
    end
    rep.amp_done(p) = max(rep.amp_done(p), nr);
    total = total + nr;
end
lines = [cell(1, 0), rep.ri{pos}];
end

function p = node_pos(st, key)
if isfield(st, 'node_of') && isKey(st.node_of, key)
    p = st.node_of(key);
else
    p = find(strcmp(st.L.key, key), 1);
end
end

function lines = root_lines(st, key, store, A, p)
% One roots_index row per root of one node: key, coordinate, routes,
% observed shape, maxres, positivity minimum, discoveries, representative.
co = st.L.coordinate(p, :);
lines = cell(1, numel(store.roots));
for k = 1:numel(store.roots)
    r = store.roots(k);
    a = attempt_row(A, r.representative_attempt_id);
    if isempty(a), mr = NaN; pm = NaN; ob = ''; else, mr = a.maxres; pm = a.positivity_min; ob = a.observed_family; end
    lines{k} = sprintf('"%s","%s","%s","%s","%s",%.6g,%.6g,%d,%s', sprintf('%s/%s', key, r.root_id), key, ...
        mat2str(co, 10), strjoin(unique({r.discoveries.route}), ';'), ob, mr, pm, ...
        numel(r.discoveries), r.representative_attempt_id);
end
end

function a = attempt_row(A, id)
% the attempt row of id (ids are A###### = the row number; else a search,
% the last match as the former id -> row map kept)
a = [];
k = sscanf(id, 'A%d');
if isscalar(k) && k >= 1 && k <= numel(A) && k == fix(k) && strcmp(A(k).attempt_id, id)
    a = A(k); return
end
j = find(strcmp({A.attempt_id}, id), 1, 'last');
if ~isempty(j), a = A(j); end
end

function h = roots_header()
h = 'root_id,key,coordinate,routes,observed_family,maxres,positivity_min,n_discoveries,representative_attempt_id';
end

function fp = guess_fp_of(st)
% Report-only guess fingerprint per attempt row ('' when unknown)
fp = repmat({''}, 1, numel(st.attempts));
if isfield(st, 'guess_fp')
    m = min(numel(st.guess_fp), numel(fp));
    fp(1:m) = st.guess_fp(1:m);
    fp(cellfun(@isempty, fp)) = {''};
end
end

function R = all_roots(st)
% Every root of every shard, keys in lattice order; root_id is the wrapper
% name <key>/R######.
R = struct('key', {}, 'root_id', {}, 'root', {});
keys = shard_keys(st);
for k = 1:numel(keys)
    store = st.shards(keys{k});
    for r = 1:numel(store.roots)
        R(end+1) = struct('key', keys{k}, 'root_id', sprintf('%s/%s', keys{k}, ...
            store.roots(r).root_id), 'root', store.roots(r)); %#ok<AGROW>
    end
end
end

function keys = shard_keys(st)
keys = st.shards.keys;
if isempty(keys), return; end
pos = zeros(1, numel(keys));         % a loop, not a function handle over st (the run context)
for k = 1:numel(keys), pos(k) = node_pos(st, keys{k}); end
[~, o] = sort(pos);
keys = keys(o);
end

function C = route_counts(st)
% From root MEMBERSHIP, never a distance threshold.
R = all_roots(st);
C = struct('n_roots', numel(R), 'n_multi_route_roots', 0, 'n_cross_route_pairs', 0, ...
    'routes', {{}});
for k = 1:numel(R)
    rt = {R(k).root.discoveries.route};
    C.routes = union(C.routes, rt);
    if numel(unique(rt)) >= 2, C.n_multi_route_roots = C.n_multi_route_roots + 1; end
    for a = 1:numel(rt)
        for b = a+1:numel(rt)
            C.n_cross_route_pairs = C.n_cross_route_pairs + ~strcmp(rt{a}, rt{b});
        end
    end
end
end

function write_pairs(st, ~)
% The final scheduler state (queue.mat holds the state stamped at the
% last call's start, for the replay comparison; the pilot report needs the
% final edges). Derived, never compared on replay.
if isfield(st, 'Q')
    sa_write('mat', fullfile(st.state_dir, 'queue_final.mat'), struct('attempt_count', ...
        numel(st.attempts), 'Q', st.Q, 'ring_calls', st.ring_calls, 'stop_reason', st.stop_reason));
end
P = st.pairs;
lines = {'kind,a,b,target_key,family_a,family_b,root_a,root_b,scaled_rms_v1,relsup_identity,relsup_all,ratio,in_near_threshold_band,route_a,route_b'};
band = st.S.identity.near_threshold;
for k = 1:numel(P)
    p = P(k);
    % root ids are <key>/R###### and keys hold commas: quoted like target_key
    lines{end+1} = sprintf('%s,%s,%s,"%s",%s,%s,"%s","%s",%.6g,%.6g,%.6g,%.6g,%d,%s,%s', p.kind, p.a, p.b, ...
        p.target_key, p.family_a, p.family_b, p.root_a, p.root_b, p.scaled_rms, ...
        p.relsup_identity, p.relsup_all, p.ratio, ...
        p.scaled_rms >= band(1) && p.scaled_rms < band(2), p.route_a, p.route_b); %#ok<AGROW>
end
sa_write('text', fullfile(st.state_dir, 'identity_pairs.csv'), [strjoin(lines, newline) newline]);
end

function gn = guess_nodes_of(st)
% Guess node count per attempt row, report-only; kept
% beside the rows (never a row field, so older records replay unchanged)
gn = NaN(1, numel(st.attempts));
if isfield(st, 'guess_nodes')
    m = min(numel(st.guess_nodes), numel(gn));
    gn(1:m) = st.guess_nodes(1:m);
end
end

function cols = csv_cols()
cols = {'attempt_id','kind','parent_attempt','family','target_key','seed_left', ...
    'seed_right','outcome','failing','flags','corrector_class','convergence_class', ...
    'wall_s','returned_mesh','maxres','positivity_min','positivity_row','positivity_xi', ...
    'nodal_min','uniform_grid_min','residual_scaled_sup','residual_scaled_rms', ...
    'residual_unscaled_sup','seam_concentration','seam_flux','s_max', ...
    'observed_family','root_id','store_status','representative_distance','file', ...
    'route','phase','lineage','edge_id','fraction','parent_root','retry_of', ...
    'positivity_sampled_min','residual_bvp4c','nmax','cut_by_cap','max_proposed_mesh'};
end

function h = csv_header()
h = [strjoin(csv_cols(), ',') ',guess_fingerprint,guess_nodes,reason'];
end

function line = csv_line(a, cols, has, gnk, fpk, reason)
% one attempts.csv row (has(c): the rows carry column c; legacy rows NaN)
parts = cell(1, numel(cols));
for c = 1:numel(cols)
    if has(c), v = a.(cols{c}); else, v = NaN; end   % legacy rows
    if ischar(v)
        parts{c} = ['"' strrep(v, '"', '""') '"'];
    elseif islogical(v)
        parts{c} = sprintf('%d', v);
    else
        parts{c} = sprintf('%.10g', v);
    end
end
parts{end+1} = fpk;
parts{end+1} = sprintf('%.10g', gnk);
parts{end+1} = ['"' strrep(reason, '"', '""') '"'];
line = strjoin(parts, ',');
end

% =========================================================================

function t = summary_text(st, wall, X)
% X (from render): reasons (plain reason per row), fl (failing checks of
% the returned calls, row order), n_roots, amp (running max amplification)
S = st.S; A = st.attempts;
L = {};
L{end+1} = st.clean_line;
L{end+1} = sprintf('seeded_atlas_run summary  (%s)', S.software_version);
L{end+1} = sprintf('model %s  mode %s  run_dir %s', S.model, S.mode, S.run_dir);
if isfield(S.lattice, 'points_per_axis') && ~isempty(S.lattice.points_per_axis)
    L{end+1} = sprintf('lattice step %s (derived from points_per_axis %s)', mat2str(S.lattice.step, 6), ...
        mat2str(S.lattice.points_per_axis));
end
L{end+1} = sprintf('stop reason: %s   corrector calls: %d of %d   wall %.1f s', ...
    st.stop_reason, st.calls, S.budget.max_calls, wall);
L{end+1} = '';
L{end+1} = 'Seedability census (nodes per family):';
for f = 1:numel(st.C.names)
    L{end+1} = sprintf('  %-12s %6d', st.C.names{f}, st.C.counts(f)); %#ok<AGROW>
end
L{end+1} = sprintf('  nodes in native range %d, outside %d', sum(st.C.in_range), sum(~st.C.in_range));
L{end+1} = '';
L{end+1} = sprintf('Targets: %d declared (policy %s), %d policy-skipped', ...
    numel(st.T.list), S.targets.policy, numel(st.T.skipped));
if any(strcmp(S.targets.policy, {'per_seed', 'per_seed_centre'}))
    L{end+1} = sprintf(['  per_seed: %d library seeds, one cold call each; %d seedable (family, node) pairs ' ...
        'not called cold (settings_effective.mat targets.skipped)'], numel(st.T.list), ...
        sum(startsWith({st.T.skipped.reason}, 'per_seed:')));
end
mnote = multiplicity_note(S);       % report-only
if ~isempty(mnote), L{end+1} = ['  ' mnote]; end
cold = A(strcmp({A.kind}, 'cold'));
L{end+1} = sprintf('%-12s %8s %9s %8s %10s %11s', 'family', 'targets', 'accepted', 'invalid', 'no_return', 'unfinished');
for f = 1:numel(st.T.families)
    fam = st.T.families{f};
    nt = sum(strcmp({st.T.list.family}, fam));
    a = cold(strcmp({cold.family}, fam));
    na = sum(strcmp({a.outcome}, 'accepted')); ni = sum(strcmp({a.outcome}, 'invalid'));
    nn = sum(strcmp({a.outcome}, 'no_return'));
    L{end+1} = sprintf('%-12s %8d %9d %8d %10d %11d', fam, nt, na, ni, nn, nt - numel(a)); %#ok<AGROW>
end
L{end+1} = sprintf('Roots in store: %d (from %d accepted cold returns)', X.n_roots, ...
    sum(strcmp({cold.outcome}, 'accepted')));
if S.routes.continuation || S.routes.rescue
    L{end+1} = sprintf(['  all routes: %d accepted node returns (%d cold, %d continuation, ' ...
        '%d rescue); %d substep records; %d interrupted records'], ...
        sum(strcmp({A.outcome}, 'accepted')), sum(strcmp({A.outcome}, 'accepted') & strcmp({A.kind}, 'cold')), ...
        sum(strcmp({A.outcome}, 'accepted') & strcmp({A.kind}, 'continuation')), ...
        sum(strcmp({A.outcome}, 'accepted') & strcmp({A.kind}, 'rescue')), ...
        sum(startsWith({A.outcome}, 'substep')), sum(strcmp({A.outcome}, 'interrupted')));
end
% Cold calls cut by the cold mesh cap
if isfield(S, 'cold') && isfinite(S.cold.mesh_cap)
    nc = 0; if isfield(A, 'cut_by_cap') && ~isempty(cold), nc = sum([cold.cut_by_cap]); end
    L{end+1} = sprintf(['cold calls cut by the cap (NMax = mesh_reject = %d): %d of %d -- a refusal of the ' ...
        'capped search, not evidence that no root exists%s'], S.cold.mesh_cap, nc, numel(cold), cut_exactness(cold));
end
% Mesh proposals recorded
L{end+1} = sprintf('mesh proposals recorded for %d of %d calls', n_recorded(A), numel(A));
L{end+1} = '';
% failing-check counts (one concatenation of the per-row lists)
fl = X.fl;
if isempty(fl), fl = {}; end
L{end+1} = sprintf('Failing checks (count over returned calls; %d call(s) returned nothing):', ...
    sum(strcmp({A.outcome}, 'no_return')));
if isempty(fl), L{end+1} = '  none'; end
[u, ~, j] = unique(fl);
for k = 1:numel(u)
    L{end+1} = sprintf('  %-28s %d', u{k}, sum(j == k)); %#ok<AGROW>
end
% plain-words reasons (first failing check or outcome); derived, never stored
pr = reshape(X.reasons, 1, []);
pr = pr(~ismember({A.outcome}, {'accepted', 'refined_valid', 'substep_valid'}));
L{end+1} = sprintf('Why calls failed (plain words; %d of %d calls):', numel(pr), numel(A));
if isempty(pr), L{end+1} = '  none'; end
[u, ~, j] = unique(pr);
cnt = accumarray(j(:), 1);
[~, o] = sort(-cnt);
for k = o(:).'
    L{end+1} = sprintf('  %5d  %s', cnt(k), u{k}); %#ok<AGROW>
end
% returns positive at the mesh nodes but failing the sampled positivity rule
rc5 = A(arrayfun(@(a) a.nodal_min > 0 && ~(a.positivity_min > 0), A));
L{end+1} = sprintf('Positive at the mesh nodes but not between them (between-node positivity check): %d', numel(rc5));
for k = 1:numel(rc5)
    L{end+1} = sprintf('  %s %s %s  nodal min %.3g  between-node min %.3g at %s xi=%.6g', ...
        rc5(k).attempt_id, rc5(k).family, rc5(k).target_key, rc5(k).nodal_min, ...
        rc5(k).positivity_min, rc5(k).positivity_row, rc5(k).positivity_xi); %#ok<AGROW>
end
% distributions
acc = cold(strcmp({cold.outcome}, 'accepted'));
L{end+1} = '';
L{end+1} = 'Distributions over returned calls (min / median / max):';
ret = A(~strcmp({A.outcome}, 'no_return') & ~strcmp({A.outcome}, 'interrupted'));
L{end+1} = dist_line('returned mesh points', [ret.returned_mesh]);
L{end+1} = dist_line('seam concentration (returned)', [ret.seam_concentration]);
L{end+1} = dist_line('seam flux (returned)', [ret.seam_flux]);
L{end+1} = dist_line('seam concentration (accepted)', [acc.seam_concentration]);
L{end+1} = dist_line('seam flux (accepted)', [acc.seam_flux]);
L{end+1} = dist_line('scaled residual sup (accepted)', [acc.residual_scaled_sup]);
L{end+1} = sprintf('  residual flags: %d of %d returned calls', ...
    sum(~cellfun(@isempty, {ret.flags})), numel(ret));
% field scale and amplification
L{end+1} = '';
fs = st.K.field_scale;
L{end+1} = sprintf('Field scale (%s): %s', S.identity.metric_id, mat2str(fs.', 4));
if X.n_roots > 0
    L{end+1} = sprintf('Per-field amplification max|value|/field_scale over roots: max %s', ...
        mat2str(X.amp.', 4));
end
% refinement
L{end+1} = '';
L{end+1} = sprintf('Refinement re-solves at tighter tolerances (expectation scaled_rms <= %g):', S.validity.refine_expected_distance);
if isempty(st.refinement), L{end+1} = '  none'; end
for k = 1:numel(st.refinement)
    r = st.refinement(k);
    L{end+1} = sprintf('  %-10s %s -> %s  %s  scaled_rms %.3g  relsup %.3g/%.3g  %s', r.family, ...
        r.parent_attempt, r.refine_attempt, r.status, r.scaled_rms, r.relsup_identity, ...
        r.relsup_all, r.failing); %#ok<AGROW>
end
% identity pairs
if ~isempty(st.pairs)
    band = S.identity.near_threshold;
    v = [st.pairs.scaled_rms];
    nb = sum(v >= band(1) & v < band(2));
    L{end+1} = '';
    L{end+1} = sprintf(['Identity pairs: %d (same-target and refinement); %d in the ' ...
        'near-threshold band [%g, %g).'], numel(v), nb, band(1), band(2));
    % Cross-route identity from root MEMBERSHIP (complete link at
    % similarity_tolerance in the store), never from a distance threshold.
    % Refinement returns are never discoveries, so they never count.
    rc = route_counts(st);
    if nb == 0
        L{end+1} = sprintf(['  near-threshold band empty over %d pair(s). This is necessary ' ...
            'but not sufficient for relying on similarity_tolerance (not yet calibrated).'], numel(v));
    else
        L{end+1} = sprintf('  near-threshold band NOT empty (%d pair(s)).', nb);
    end
    L{end+1} = sprintf(['  cross-route identity by root membership: %d root(s) whose ' ...
        'discoveries span >= 2 routes; %d different-route discovery pair(s) inside one root.'], ...
        rc.n_multi_route_roots, rc.n_cross_route_pairs);
    L{end+1} = '  similarity_tolerance remains PROVISIONAL.';
end
% search report (node classes, edges, lineages, rescue, ring)
if (S.routes.continuation || S.routes.rescue) && isfield(st, 'Q') && ~isempty(st.Q)
    Lc = sa_continuation('report_lines', st);
    if ~isempty(mnote)
        m = ~cellfun(@isempty, regexp(Lc, '^  nodes with >= 2 distinct roots: \d+$', 'once'));
        Lc(m) = strcat(Lc(m), ' (lower bound: multiplicity not tested)');
    end
    L = [L, Lc];
end
% unfinished
L{end+1} = '';
done = [cold.target_number];
unf = setdiff(1:numel(st.T.list), done);
L{end+1} = sprintf('Unfinished targets: %d', numel(unf));
for k = unf
    L{end+1} = sprintf('  %s %s', st.T.list(k).family, st.T.list(k).key); %#ok<AGROW>
end
t = [strjoin(L, newline) newline];
end

function t = multiplicity_note(S)
% Rule: multiplicity is tested only when the cold targets
% are every seedable (family, node) pair (targets.policy 'all_seedable');
% for any other policy the roots found are those connected to the cold-
% seeded nodes, so a node's root count is a lower bound.
t = '';
if strcmp(S.targets.policy, 'all_seedable')
    % The thinning lever tests multiplicity only
    % partially
    if isfield(S.targets, 'thin') && ~isempty(S.targets.thin.families)
        o = S.targets.thin.offset; if isempty(o), o = 0; end
        t = sprintf(['MULTIPLICITY PARTIALLY TESTED: cold targets thinned for %s (stride %d, offset %s); ' ...
            'coexisting roots of those families at thinned nodes are looked for only by continuation, so ' ...
            'nodes with >= 2 roots is a lower bound.'], strjoin(S.targets.thin.families, ', '), ...
            S.targets.thin.stride, mat2str(o));
    end
    return
end
if any(strcmp(S.targets.policy, {'per_seed', 'per_seed_centre'}))
    t = ['MULTIPLICITY NOT TESTED BY COLD: cold targets = per_seed (one cold call per library seed, at the ' ...
        'seedable node nearest its library cell; the anchor node excepted). Other roots at the seed''s other ' ...
        'nodes are looked for only by continuation (capped lineages) and the per-root pass (neighbours of ' ...
        'roots found); rescue only fills nodes with no root. So nodes with >= 2 roots is a lower bound. ' ...
        'Use targets.policy all_seedable for a cold multiplicity census.'];
    return
end
switch S.targets.policy
    case 'central_nine', what = 'the central nine nodes';
    case 'pilot_sample', what = 'the central nine and the pilot sample nodes';
    otherwise, what = sprintf('the %s cold targets', S.targets.policy);
end
t = sprintf(['MULTIPLICITY NOT TESTED: cold targets = %s (not all_seedable). The roots found are the branch ' ...
    'connected to %s only; coexisting roots elsewhere are not looked for, so nodes with >= 2 roots is a ' ...
    'lower bound. Use all_seedable for a multiplicity census.'], S.targets.policy, what);
end

function n = n_recorded(A)
n = 0;
if isfield(A, 'max_proposed_mesh') && ~isempty(A), n = sum(isfinite([A.max_proposed_mesh])); end
end

function t = cut_exactness(cold)
% The cut count is exact when every cold call's proposals are recorded
% (the live label check then holds the label to the exact rule); otherwise
% a lower bound of what the cap changed
c = cold(~strcmp({cold.outcome}, 'interrupted'));
nr = n_recorded(c);
if isempty(c)
    t = '';
elseif nr == numel(c)
    t = '; count exact (mesh proposals recorded for every cold call)';
else
    t = sprintf(['; count a lower bound of the calls the cap changed (mesh proposals not recorded for %d of ' ...
        '%d cold calls)'], numel(c) - nr, numel(c));
end
end

% =========================================================================
% Exact cut recount

function X = cutset(run_dir, c, kinds, o)
if nargin < 3 || isempty(kinds), kinds = {'cold'}; end
if ischar(kinds), kinds = {kinds}; end
if nargin < 4, o = struct(); end
refuse_if_compacted(run_dir, 'sa_report(''cutset'')');
run_dir = sa_layout('state', run_dir);
if ~isfield(o, 'out_dir'), o.out_dir = run_dir; end
if ~isfield(o, 'write'), o.write = true; end
if ~isfield(o, 'cut_time_s'), o.cut_time_s = 0; end
f = dir(fullfile(run_dir, 'attempts', 'A*.mat'));
names = sort({f.name});
names = names(~cellfun(@isempty, regexp(names, '^A\d{6}\.mat$', 'once')));
C = struct('attempt_id', {}, 'kind', {}, 'outcome', {}, 'nmax', {}, 'returned_mesh', {}, 'initial', {}, ...
    'n_proposals', {}, 'P', {}, 'over_nmax', {}, 'cut_by_cap', {}, 'convergence_class', {}, 'wall_s', {}, ...
    'returned', {});
nint = 0;
for k = 1:numel(names)
    file = fullfile(run_dir, 'attempts', names{k});
    v = who('-file', file);
    Z = load(file, v{ismember(v, {'attempt', 'opts', 'mesh_trace'})});
    a = Z.attempt;
    if strcmp(a.outcome, 'interrupted'), nint = nint + 1; continue; end
    if ~any(strcmp(a.kind, kinds)), continue; end
    r = struct('attempt_id', a.attempt_id, 'kind', a.kind, 'outcome', a.outcome, 'nmax', NaN, ...
        'returned_mesh', a.returned_mesh, 'initial', NaN, 'n_proposals', NaN, 'P', NaN, 'over_nmax', false, ...
        'cut_by_cap', isfield(a, 'cut_by_cap') && isequal(a.cut_by_cap, true), ...
        'convergence_class', a.convergence_class, 'wall_s', a.wall_s, 'returned', isfinite(a.returned_mesh));
    if isfield(Z, 'opts') && isstruct(Z.opts) && isfield(Z.opts, 'nmax'), r.nmax = Z.opts.nmax; end
    if isfield(Z, 'mesh_trace') && isstruct(Z.mesh_trace)
        r.initial = Z.mesh_trace.initial; r.n_proposals = Z.mesh_trace.n_proposals;
        r.P = sa_mesh_trace('max_of', Z.mesh_trace); r.over_nmax = Z.mesh_trace.over_nmax;
    end
    C(end+1) = r; %#ok<AGROW>
end
cls = arrayfun(@(r) cut_class(r, c), C, 'UniformOutput', false);
X = struct('c', c, 'kinds', {kinds}, 'calls', C, 'class', {cls}, 'n_interrupted_excluded', nint);
X.n = numel(C);
X.n_cut = sum(strcmp(cls, 'cut')); X.n_not_cut = sum(strcmp(cls, 'not_cut'));
X.n_unknown_not_recorded = sum(strcmp(cls, 'unknown_not_recorded'));
X.n_unknown_above_nmax = sum(strcmp(cls, 'unknown_above_nmax'));
X.exact = X.n_unknown_not_recorded == 0 && X.n_unknown_above_nmax == 0;
X.cut_ids = {C(strcmp(cls, 'cut')).attempt_id};
L = {};
L{end+1} = sprintf(['cutset c = %d (kinds %s; %s): %d calls; cut %d; not_cut %d; unknown_not_recorded %d; ' ...
    'unknown_above_nmax %d; interrupted records excluded %d; %s'], c, strjoin(kinds, ','), run_dir, X.n, X.n_cut, ...
    X.n_not_cut, X.n_unknown_not_recorded, X.n_unknown_above_nmax, nint, ...
    tern(X.exact, 'EXACT', 'NOT exact (unknown classes present)'));
% histogram of P
P = [C.P];
edges = [0 1 500 1000 2000 3601 5000 10001 20001 Inf];
L{end+1} = sprintf('histogram of P (largest proposed mesh) over the %d calls (NaN = not recorded: %d):', ...
    numel(P), sum(isnan(P)));
L{end+1} = sprintf('  P = 0 (no proposal): %d', sum(P == 0));
for b = 2:numel(edges) - 1
    L{end+1} = sprintf('  %6g <= P < %-6g: %d', edges(b), edges(b + 1), sum(P >= edges(b) & P < edges(b + 1))); %#ok<AGROW>
end
% cap-loss curve
L{end+1} = 'cap-loss curve (exact rule; c = NMax = mesh_reject): c, cut calls, of them stored accepted returns, exact';
X.curve = struct('c', {}, 'n_cut', {}, 'n_cut_accepted', {}, 'exact', {}, 'ids', {});
for cc = [3601 5000 10000]
    k2 = arrayfun(@(r) cut_class(r, cc), C, 'UniformOutput', false);
    cut = strcmp(k2, 'cut');
    acc = cut & strcmp({C.outcome}, 'accepted');
    ex = ~any(strcmp(k2, 'unknown_not_recorded') | strcmp(k2, 'unknown_above_nmax'));
    X.curve(end+1) = struct('c', cc, 'n_cut', sum(cut), 'n_cut_accepted', sum(acc), 'exact', ex, ...
        'ids', {{C(cut).attempt_id}});
    L{end+1} = sprintf('  c = %5d: cut %d, of them accepted %d, %s', cc, sum(cut), sum(acc), ...
        tern(ex, 'exact', 'NOT exact')); %#ok<AGROW>
end
% Prediction fallback for calls without a trace
nr = isnan(P);
if any(nr)
    cert = nr & (strcmp({C.convergence_class}, 'UNRESOLVED_RELTOL_NMAX') | [C.returned_mesh] >= c);
    poss = nr & ~cert & strcmp({C.outcome}, 'no_return') & [C.wall_s] > o.cut_time_s;
    X.rc51 = struct('certain', {{C(cert).attempt_id}}, 'possible', {{C(poss).attempt_id}});
    L{end+1} = sprintf(['Prediction (lower bound) for the %d calls without recorded proposals: ' ...
        'certain %d, possible %d (no-return with wall above %g s) -- not exact'], sum(nr), sum(cert), sum(poss), ...
        o.cut_time_s);
else
    X.rc51 = struct('certain', {{}}, 'possible', {{}});
end
X.lines = L;
if o.write
    if exist(o.out_dir, 'dir') ~= 7, mkdir(o.out_dir); end
    rows = {'attempt_id,kind,outcome,nmax,returned_mesh,initial,n_proposals,max_proposed_mesh,over_nmax,class,cut_by_cap_stored,convergence_class,wall_s'};
    for k = 1:numel(C)
        r = C(k);
        rows{end+1} = sprintf('%s,%s,%s,%.10g,%.10g,%.10g,%.10g,%.10g,%d,%s,%d,%s,%.6g', r.attempt_id, r.kind, ...
            r.outcome, r.nmax, r.returned_mesh, r.initial, r.n_proposals, r.P, r.over_nmax, cls{k}, r.cut_by_cap, ...
            r.convergence_class, r.wall_s); %#ok<AGROW>
    end
    sa_write('text', fullfile(o.out_dir, sprintf('cutset_%d.csv', c)), [strjoin(rows, newline) newline]);
    sa_write('text', fullfile(o.out_dir, sprintf('cutset_%d_summary.txt', c)), [strjoin(L, newline) newline]);
end
end

% =========================================================================
% plain-words failure reasons

function T = plain_table()
% check / outcome -> plain words. Record-only checks never fail a call.
T = struct('name', {}, 'words', {}, 'kind', {});
add = @(T, n, w, k) [T, struct('name', n, 'words', w, 'kind', k)];
T = add(T, 'finite', 'the solution has non-finite values', 'check');
T = add(T, 'real_valued', 'the solution has complex values', 'check');
T = add(T, 'convergence_maxres', 'the solver did not reach its tolerance', 'check');
T = add(T, 'convergence_proposed_mesh', 'the solver needed more mesh points than allowed (NMax)', 'check');
T = add(T, 'convergence_warning', 'the solver warned that it did not meet its tolerance', 'check');
T = add(T, 'mesh', 'the returned mesh is at or above the mesh limit', 'check');
T = add(T, 'positivity', 'a concentration is zero or negative somewhere', 'check');
T = add(T, 'positivity_exact', 'a concentration is zero or negative somewhere', 'check');
T = add(T, 'seam_concentration', 'the halves do not match at the seam (value)', 'check');
T = add(T, 'seam_flux', 'the halves do not match at the seam (flux)', 'check');
T = add(T, 'mass_accumulator', 'a pool mass misses its target', 'check');
T = add(T, 'mass_dense', 'a pool mass misses its target', 'check');
T = add(T, 'mass_dense_refinement', 'a pool mass misses its target', 'check');
T = add(T, 'shift_antisymmetry', 'the forcing shifts are inconsistent', 'check');
T = add(T, 'velocity_consistency', 'the velocity estimates are inconsistent', 'check');
T = add(T, 'residual_bvp4c', 'the solver''s own residual estimate is above the tolerance', 'check');
T = add(T, 'positivity_sampled', 'record only: never fails a call', 'record_only');
T = add(T, 'residual', 'record only: never fails a call', 'record_only');
T = add(T, 'no_return:singular_jacobian', 'no solution returned: the solver met a singular Jacobian', 'outcome');
T = add(T, 'no_return:singularity_guard', 'no solution returned: the singularity guard stopped the solver', 'outcome');
T = add(T, 'no_return:mesh_blowup', 'no solution returned: the mesh grew past the limit', 'outcome');
T = add(T, 'no_return:solver_error', 'no solution returned: the solver stopped with an error', 'outcome');
T = add(T, 'cut_by_cap', ['cut by the cold mesh cap: a refusal of the capped search, not evidence that no ' ...
    'root exists'], 'outcome');
T = add(T, 'interrupted', 'the process was killed during this call; it was run again as %s', 'outcome');
T = add(T, 'substep', 'intermediate continuation step, not a lattice node', 'outcome');
T = add(T, 'accepted', '', 'outcome');
T = add(T, 'refined_valid', '', 'outcome');
end

function m = retried_as(A)
% interrupted attempt id -> the id of the call that re-ran its work unit
m = containers.Map('KeyType', 'char', 'ValueType', 'char');
for k = 1:numel(A)
    if isfield(A, 'retry_of') && ~isempty(A(k).retry_of), m(A(k).retry_of) = A(k).attempt_id; end
end
end

function w = plain_words()
% check / outcome name -> words (built once per render, not per row)
T = plain_table();
w = containers.Map({T.name}, {T.words});
end

function r = plain_reason(A, rt, w)
% the plain-words reason of one attempt row: '' for accepted and
% refined_valid returns; the first failing check otherwise. rt (optional):
% retried_as(rows), for the id that re-ran an interrupted call; w
% (optional): plain_words()
if nargin < 2, rt = []; end
if nargin < 3, w = plain_words(); end
o = A.outcome;
r = '';
if any(strcmp(o, {'accepted', 'refined_valid'})), return; end
if strcmp(o, 'interrupted')
    r = sprintf(w('interrupted'), 'the next call of the same work unit');
    if ~isempty(rt) && isKey(rt, A.attempt_id), r = sprintf(w('interrupted'), rt(A.attempt_id)); end
    return
end
pre = '';
if strcmp(o, 'substep_valid'), r = w('substep'); return; end
if startsWith(o, 'substep'), pre = [w('substep') '; ']; end
if isfield(A, 'cut_by_cap') && isequal(A.cut_by_cap, true)
    r = [pre w('cut_by_cap')]; return
end
if strcmp(o, 'no_return') || (startsWith(o, 'substep') && ~isfinite(A.returned_mesh))
    r = [pre no_return_words(A, w)]; return
end
f = strsplit(A.failing, ';');
f = f(~cellfun(@isempty, f));
if isempty(f)
    if startsWith(o, 'substep'), r = w('substep'); else, r = o; end
    return
end
name = f{1};
if startsWith(name, 'mass_') && ~isKey(w, name), name = 'mass_dense'; end
if isKey(w, name), r = [pre w(name)]; else, r = [pre name]; end
end

function t = no_return_words(A, w)
why = lower(A.corrector_reason);
if contains(why, 'singular jacobian')
    t = w('no_return:singular_jacobian');
elseif contains(why, 'singularity_guard')
    t = w('no_return:singularity_guard');
elseif contains(why, 'mesh blow-up')
    t = w('no_return:mesh_blowup');
else
    t = w('no_return:solver_error');
end
end

function refuse_if_compacted(run_dir, tool)
D = sa_layout('state', run_dir);
f = fullfile(D, 'compaction.mat');
if exist(f, 'file') == 2
    Z = load(f, 'compacted_at', 'n_removed');
    error('SA:REPORT:E_COMPACTED', ['%s refused: %s was compacted on %s: %d files removed (the attempt ' ...
        'files it would enumerate are gone; it refuses rather than miscounts).'], tool, run_dir, ...
        Z.compacted_at, Z.n_removed);
end
end

function k = cut_class(r, c)
% the exact rule for the bound pair (c, c)
if isnan(r.P), k = 'unknown_not_recorded'; return; end
over_m = r.returned && r.returned_mesh >= c;
if ~(c > r.nmax)                         % c <= the call's own NMax (or NMax unknown: NaN -> exact form)
    if r.P > c || over_m, k = 'cut'; else, k = 'not_cut'; end
    return
end
if r.P > c, k = 'cut'; return; end
if ~r.over_nmax
    if over_m, k = 'cut'; else, k = 'not_cut'; end
else
    k = 'unknown_above_nmax';
end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function s = dist_line(name, v)
v = v(isfinite(v));
if isempty(v)
    s = sprintf('  %-32s (none)', name);
else
    s = sprintf('  %-32s %.3g / %.3g / %.3g  (n=%d)', name, min(v), median(v), max(v), numel(v));
end
end
