function varargout = sa_search(action, varargin)
%SA_SEARCH  Scheduler and executor.
%
%   ctx = SA_SEARCH('run', ctx)
%       Passes: cold -> import (milestone 10: one call per planned row of
%       ctx.import, sa_import; S.solver, no cold cap, no retry; then the
%       source-edge knowledge ctx.ik) -> refinement ->
%       [continuation -> rescue] up to search.max_passes (repeated only while
%       rescue added roots), per phase: 'interior' then 'ring' when
%       search.ring_last (the ring = lattice nodes outside the library's
%       native range), else one phase 'all'. Routes choose only the guess and
%       the target coordinate: every return, substeps included, crosses
%       sa_validity('evaluate') with the run's solver settings.
%   [ctx, A, R] = SA_SEARCH('call', ctx, unit, tg, guess, solver)
%       THE executor. Budget check first (budget_calls, budget_wall,
%       budget_ring); then either REPLAY (resume: the next stored attempt must
%       carry the same work unit; its verdict is re-derived; its attempt
%       row and the store must replay bit-identically) or a LIVE call: journal
%       'start', one M.solve, the attempt file (atomic), the node's shard
%       (an accepted call only), journal 'end', a checkpoint every
%       checkpoint.every_calls (the checkpoint files attempts.mat,
%       queue.mat, process_time.mat, which a resume compares with its replay
%       and never starts from; the disk floor) and,
%       when report.files_interval_s has passed (checked after each live
%       call, no timer object), the human-facing files attempts.csv,
%       summary.txt and roots_index.csv from the render cache.
%   u = SA_SEARCH('unit')            blank work unit (fixed field order)
%   k = SA_SEARCH('unit_key', u)     its deterministic text form
%   s = SA_SEARCH('solver_for', S, kind)
%   o = SA_SEARCH('opts_for', ctx, solver)   options the module builds (zero
%                                    solves: a probe corrector returns nothing)
%   ctx = SA_SEARCH('checkpoint', ctx)
%   ctx = SA_SEARCH('finish', ctx)   final derived files (pairs, the state,
%                                    the human-facing files with the
%                                    cached-vs-full render self-check)
%   sol = SA_SEARCH('get_sol', ctx, attempt_id)   a stored call's sol (the
%                                    refinement parent, continuation parents,
%                                    the light-shard refill)
%   Memory: no sol is kept in memory between calls (the refinement loads
%   its parent by attempt id, as a resumed run always did), and the shards
%   are kept light (sa_record_call, flush_shards). Root counts per node:
%   ctx.nroots (checked against the shards at every checkpoint).
%   g = SA_SEARCH('free_gib', dir)
%
%   Determinism: the queue is a pure function of settings and completed
%   outcomes; nothing here reads a clock except the budget_wall check, which
%   applies to live calls only.
%   Every call also records a report-only guess fingerprint
%   (sa_modes 'fingerprint': guess x, y bits, coordinate, solver options)
%   beside the attempt rows, like guess_nodes (never a row field or key).
%   The cold pass calls with solver_for(S,
%   'cold') (cold.mesh_cap: NMax = mesh_reject = cap; continuation, rescue
%   and the pass keep S.solver; refinement as before); per phase and search
%   pass the order is lineages (expand) -> the per-root pass (sa_root_pass,
%   when root_pass.enabled) -> rescue -> closure while rescue added a node;
%   the journal 'end' line of a cold call cut by the cap carries the suffix
%   ' cut_by_cap'.
%   Screen output (sa_progress; display only, never a record or a key): a
%   stage line when each stage starts and a status line every
%   report.progress_interval_s seconds, checked after every call (live or
%   replayed). The live call runs inside evalc: the corrector, its warning
%   state and its lastwarn capture run exactly as without it; only the text
%   it would print is captured and passed to sa_progress('solver_text'),
%   which leaves out bvp4c's mesh-limit warning unless
%   report.solver_warnings is true and prints any other text unchanged.

switch action
    case 'run',        varargout{1} = run(varargin{:});
    case 'call',       [varargout{1}, varargout{2}, varargout{3}] = do_call(varargin{:});
    case 'unit',       varargout{1} = blank_unit();
    case 'unit_key',   varargout{1} = unit_key(varargin{:});
    case 'solver_for', varargout{1} = solver_for(varargin{:});
    case 'opts_for',   varargout{1} = opts_for(varargin{:});
    case 'checkpoint', varargout{1} = checkpoint(varargin{:});
    case 'finish',     varargout{1} = finish(varargin{:});
    case 'free_gib',   varargout{1} = free_gib(varargin{:});
    case 'get_sol',    varargout{1} = get_sol(varargin{:});
    case 'end_text',   varargout{1} = end_text(varargin{:});
    otherwise, error('SA:SEARCH', 'Unknown sa_search action %s.', action);
end
end

% =========================================================================
% orchestration

function ctx = run(ctx)
S = ctx.S;
% disk floor at the start (checked again at every checkpoint)
if ~ctx.replay.active && free_gib(S.run_dir) < S.checkpoint.min_free_gib
    ctx = set_stop(ctx, 'disk_floor');
end
% ---- cold pass (budget check, guess, one call per target) ----------------
if ~ctx.stop && S.routes.cold
    stage(ctx, 'cold search (one call per cold target)', 'cold');
    for t = 1:numel(ctx.T.list)
        tg = ctx.T.list(t);
        u = cold_unit(tg, t);
        [ok, why] = budget_ok(ctx, u);
        if ~ok, ctx = set_stop(ctx, why); break; end
        [guess, gwhy] = ctx.M.pair_guess(ctx.lib, ctx.theta, tg.left_type, tg.cell_left, ...
            tg.right_type, tg.cell_right, S.model_options);
        if isempty(guess)
            error('SA:RUN', 'Seedable target %s %s produced no guess: %s', tg.family, tg.key, gwhy);
        end
        ctx = do_call(ctx, u, tg, guess, solver_for(S, 'cold'));   % cold mesh cap
        if ctx.stop, break; end
    end
end
% ---- import (one call per planned source root; milestone 10) -------------
if ~ctx.stop && sa_import('active', ctx)
    stage(ctx, 'import (one call per source root)', 'import');
    rows = ctx.import.rows(strcmp({ctx.import.rows.status}, 'planned'));
    for t = 1:numel(rows)
        [u, tg] = sa_import('unit', ctx, rows(t));
        [ok, why] = budget_ok(ctx, u);
        if ~ok, ctx = set_stop(ctx, why); break; end
        guess = sa_import('guess', ctx, rows(t));      % SHA-256 checked, live and on replay
        ctx = do_call(ctx, u, tg, guess, S.solver);     % no cold cap, one call, no retry
        if ctx.stop, break; end
    end
end
if isfield(ctx, 'import')
    ctx.ik = sa_import('knowledge', ctx);              % source-edge knowledge ('none': inactive)
end
% ---- refinement (budget refusals are recorded, not a stop) ---------------
ctx.refinement = sa_report('blank_refinement');
if S.validity.refine_enabled
    stage(ctx, 'refinement (one re-solve per family at tighter tolerances)', 'refinement');
    ctx = sa_refinement(ctx);           % parents' sols loaded by attempt id
end
% ---- continuation / rescue passes per phase ------------------------------
if ~ctx.stop && (S.routes.continuation || S.routes.rescue)
    ctx.Q = sa_continuation('init', ctx);
    if S.search.ring_last
        phases = {'interior'};
        if any(ctx.Q.box & ctx.Q.ring), phases{end+1} = 'ring'; end
    else
        phases = {'all'};
    end
    for p = 1:numel(phases)
        ctx.Q.phase = phases{p};
        ctx.Q.phases_run{end+1} = phases{p};
        for pass = 1:S.search.max_passes
            ctx.Q.pass = pass;
            where = sprintf('%s phase, search pass %d', phase_name(phases{p}), pass);
            if S.routes.continuation
                stage(ctx, sprintf('continuation (%s)', where), 'continuation');
                ctx = sa_continuation('expand', ctx);
                if ctx.stop, break; end
                if isfield(S, 'root_pass') && S.root_pass.enabled
                    stage(ctx, sprintf('per-root pass (%s)', where), 'root_pass');
                    ctx = sa_root_pass('run', ctx);          % per-root pass
                    if ctx.stop, break; end
                end
            end
            added = 0;
            if S.routes.rescue
                stage(ctx, sprintf('rescue (%s)', where), 'rescue');
                [ctx, added] = sa_rescue('run', ctx);
                if ctx.stop, break; end
            end
            if added == 0, break; end
        end
        if ctx.stop, break; end
    end
end
if ~ctx.stop
    ctx.stop_reason = 'completed';
end
end

function u = cold_unit(tg, t)
u = blank_unit();
u.route = 'cold'; u.kind = 'cold'; u.phase = 'cold';
u.target_number = t; u.family = tg.family; u.key = tg.key; u.index = tg.index;
u.fraction = 1; u.coordinate = tg.coordinate;
u.guess_source = 'library';
u.seed_left = sprintf('%s(%d,%d)', tg.left_type, tg.cell_left);
u.seed_right = sprintf('%s(%d,%d)', tg.right_type, tg.cell_right);
end

% =========================================================================
% work units

function u = blank_unit()
u = struct('route', '', 'kind', '', 'phase', '', 'pass', 0, 'target_number', NaN, ...
    'family', '', 'key', '', 'index', [], 'fraction', 1, 'coordinate', [], ...
    'lineage', 0, 'edge_id', 0, 'parent_attempt', '', 'parent_root', '', ...
    'guess_source', '', 'seed_left', '', 'seed_right', '', 'retry_of', '');
end

function k = unit_key(u)
% Everything that defines the work, never retry_of (a retry is the same unit).
k = sprintf('%s|%s|%s|p%d|t%s|%s|%s|i%s|f%s|c%s|L%d|E%d|%s|%s|%s|%s|%s', u.route, u.kind, ...
    u.phase, u.pass, num2str(u.target_number, 17), u.family, u.key, mat2str(u.index, 17), ...
    num2str(u.fraction, 17), mat2str(u.coordinate, 17), u.lineage, u.edge_id, ...
    u.parent_attempt, u.parent_root, u.guess_source, u.seed_left, u.seed_right);
k = strrep(k, newline, ' ');
end

function s = solver_for(S, kind)
% The solver options of a unit kind. Refinement: refinement
% tolerances, nmax = mesh_reject. A cold call under a finite
% cold.mesh_cap gets nmax = mesh_reject = cap; a legacy settings struct
% without the cold group, or cap Inf, leaves S.solver unchanged (so stored
% older options re-derive as they ran). Every other kind (continuation,
% rescue, the per-root pass) keeps S.solver.
s = S.solver;
if strcmp(kind, 'refinement')
    s.reltol = S.validity.refine_reltol;
    s.abstol = S.validity.refine_abstol;
    s.nmax = S.solver.mesh_reject;
elseif strcmp(kind, 'cold') && isfield(S, 'cold') && isstruct(S.cold) && isfield(S.cold, 'mesh_cap') && ...
        isfinite(S.cold.mesh_cap)
    s.nmax = S.cold.mesh_cap;
    s.mesh_reject = S.cold.mesh_cap;
end
end

function o = opts_for(ctx, solver)
R = ctx.M.solve(ctx.theta, [], [], solver, ctx.S.model_options, @probe_corrector);
o = R.opts;
end

function [sol, info] = probe_corrector(~, ~, ~, ~)
sol = [];
info = struct('class', 'NO_RETURN', 'warning_id', '', 'reason', 'options probe (no solve)');
end

% =========================================================================
% budget and stops

function [ok, why] = budget_ok(ctx, u)
ok = true; why = '';
S = ctx.S;
if ctx.calls >= S.budget.max_calls
    ok = false; why = 'budget_calls';
elseif ~ctx.replay.active && toc(ctx.t_run) >= min(S.budget.max_wall_s, S.budget.wall_s)
    ok = false; why = 'budget_wall';
elseif strcmp(u.phase, 'ring') && ctx.ring_calls >= S.budget.ring_calls
    ok = false; why = 'budget_ring';
end
end

function ctx = set_stop(ctx, why)
ctx.stop = true;
ctx.stop_reason = why;
if ~ctx.replay.active
    sa_journal('note', ctx.state_dir, sprintf('stop %s at attempt count %d', why, numel(ctx.attempts)));
end
end

function integrity(ctx, why)
try
    sa_journal('note', ctx.state_dir, ['integrity ' why]);
catch
end
error('SA:INTEGRITY', 'Integrity stop at attempt count %d: %s', numel(ctx.attempts), why);
end

% =========================================================================
% the executor

function [ctx, A, R] = do_call(ctx, u, tg, guess, solver)
A = []; R = [];
[ok, why] = budget_ok(ctx, u);
if ~ok, ctx = set_stop(ctx, why); return; end
ctx.Q_at_last = ctx.Q;
ctx.ring_calls_at_last = ctx.ring_calls;
target = ctx.M.pool_targets(ctx.theta, u.coordinate).solver_target;
ukey = unit_key(u);
gn = guess_count(guess);            % report-only, never in the row or key
fp = sa_modes('fingerprint', guess, u.coordinate, solver);   % report-only
rec = [];
if ctx.replay.active
    [ctx, rec] = replay_next(ctx, ukey, gn, fp);
    if isempty(rec)
        % replay exhausted on this request: re-check the budget after the
        % interrupted records it consumed (they count toward max_calls)
        [ok, why] = budget_ok(ctx, u);
        if ~ok, ctx = set_stop(ctx, why); return; end
    end
end
if ~isempty(rec)
    % ---- REPLAY: zero solves, verdict re-derivation -----------------
    if ~isequal(rec.target, target)
        integrity(ctx, sprintf('%s: stored solver target differs from the scheduler''s', rec.attempt.attempt_id));
    end
    [V, info] = sa_validity('rederive', ctx.M, rec, ctx.theta, ctx.S, opts_for(ctx, solver), ctx.prov);
    if ~info.integrity_ok
        integrity(ctx, sprintf('%s: %s', rec.attempt.attempt_id, info.reason));
    end
    if ~info.same_version
        integrity(ctx, sprintf('%s was written under another software/MATLAB version or thread mode', ...
            rec.attempt.attempt_id));
    end
    R = struct('sol', rec.sol, 'opts', rec.opts, 'warning_id', rec.warning_id, ...
        'warning_msg', rec.warning_msg, 'corrector_class', rec.attempt.corrector_class, ...
        'corrector_reason', rec.attempt.corrector_reason, 'wall_s', rec.attempt.wall_s, ...
        'corrector_nodal_min', rec.attempt.corrector_nodal_min);
    % The stored mesh trace is a trusted
    % input (like stats and the warning id: it cannot be recomputed without a
    % solve); P is RECOMPUTED from its proposals, never copied from the stored
    % row, so the row isequaln check below verifies the stored field, and
    % sa_record_call re-runs the trace consistency checks.
    if isfield(rec, 'mesh_trace') && isstruct(rec.mesh_trace)
        R.mesh_trace = rec.mesh_trace;
    else
        R.mesh_trace = sa_mesh_trace('blank');
    end
    R.max_proposed_mesh = sa_mesh_trace('max_of', R.mesh_trace);
    u.retry_of = rec.unit.retry_of;
    [ctx, A] = sa_record_call(ctx, R, u, tg, target, V, false);
    ctx.guess_nodes(numel(ctx.attempts)) = gn;
    ctx = set_fp(ctx, fp);
    if ~isequaln(A, rec.attempt)
        integrity(ctx, sprintf('%s: attempt row differs from its replay', A.attempt_id));
    end
    ctx = flush_shards(ctx, false);
    ctx = compare_checkpoints(ctx);
    if ~ctx.replay.active
        ctx = end_of_replay(ctx);   % the last stored attempt has been processed
    end
else
    % ---- LIVE --------------------------------------------------------------
    n = numel(ctx.attempts) + 1;
    id = sprintf('A%06d', n);
    if ~isempty(ctx.pending_retry)
        if strcmp(ctx.pending_retry.key, ukey), u.retry_of = ctx.pending_retry.id; end
        ctx.pending_retry = [];
    end
    sa_journal('start', ctx.state_dir, id, ukey);
    R = quiet_solve(ctx, guess, target, solver);
    [ctx, A] = sa_record_call(ctx, R, u, tg, target, [], true);
    ctx.guess_nodes(numel(ctx.attempts)) = gn;
    ctx = set_fp(ctx, fp);
    ctx = flush_shards(ctx, true);
    sa_journal('end', ctx.state_dir, id, end_text(A));
    if mod(n, ctx.S.checkpoint.every_calls) == 0
        ctx = checkpoint(ctx);
    end
    if reports_due(ctx)
        ctx = write_reports(ctx);
    end
end
if strcmp(u.phase, 'ring'), ctx.ring_calls = ctx.ring_calls + 1; end
if sa_progress('due'), sa_progress('status', progress_info(ctx, sa_progress('kind'))); end
end

function R = quiet_solve(ctx, guess, target, solver)
% One module solve inside evalc: the call, its warning state and lastwarn
% are unchanged; the captured text goes to the screen through sa_progress.
R = []; %#ok<NASGU>
txt = evalc('R = ctx.M.solve(ctx.theta, guess, target, solver, ctx.S.model_options, ctx.corrector);');
sa_progress('solver_text', txt, report_leaf(ctx.S, 'solver_warnings', false));
end

% =========================================================================
% progress lines (screen only)

function stage(ctx, name, kind)
sa_progress('stage', name, progress_info(ctx, kind));
end

function I = progress_info(ctx, kind)
I = struct('kind', kind, 'calls', ctx.calls, 'max_calls', ctx.S.budget.max_calls, ...
    'roots', count_roots(ctx), 'within', '');
switch kind
    case 'cold'
        I.within = sprintf('%d of %d cold targets called', count_kind(ctx, 'cold'), numel(ctx.T.list));
    case 'refinement'
        I.within = sprintf('%d refinement call(s)', count_kind(ctx, 'refinement'));
    case 'import'
        I.within = sprintf('%d of %d planned imports called', count_kind(ctx, 'import'), ...
            sum(strcmp({ctx.import.rows.status}, 'planned')));
    case {'continuation', 'root_pass', 'rescue'}
        if isstruct(ctx.Q) && isfield(ctx.Q, 'box')
            [r, b] = reached(ctx);
            I.within = sprintf('box nodes with a root %d of %d; edges %d', r, b, numel(ctx.Q.edges));
            if strcmp(kind, 'rescue'), I.within = sprintf('%s; nodes rescue-tried %d', I.within, ...
                    sum(ctx.Q.rescue_done)); end
        end
end
if isfield(ctx, 'replay') && isstruct(ctx.replay) && ctx.replay.active
    I.replay = sprintf('replaying stored calls (no solves): %d of %d', ctx.replay.pos - 1, ...
        numel(ctx.replay.files));
end
end

function n = count_kind(ctx, kind)
n = 0;
if ~isempty(ctx.attempts), n = sum(strcmp({ctx.attempts.kind}, kind)); end
end

function n = count_roots(ctx)
if isfield(ctx, 'nroots'), n = sum(ctx.nroots); return; end
n = 0;
k = ctx.shards.keys;
for q = 1:numel(k), n = n + numel(ctx.shards(k{q}).roots); end
end

function [r, b] = reached(ctx)
nodes = find(ctx.Q.box);
b = numel(nodes);
if isfield(ctx, 'nroots'), r = sum(ctx.nroots(nodes) > 0); return; end
r = 0;
for q = reshape(nodes, 1, [])
    key = ctx.L.key{q};
    if isKey(ctx.shards, key) && ~isempty(ctx.shards(key).roots), r = r + 1; end
end
end

function t = phase_name(p)
switch p
    case 'interior', t = 'interior';
    case 'ring',     t = 'ring (outside the library range)';
    otherwise,       t = 'whole lattice';
end
end

function v = report_leaf(S, name, d)
v = d;
if isfield(S, 'report') && isstruct(S.report) && isfield(S.report, name), v = S.report.(name); end
end

function ctx = flush_shards(ctx, write)
% write = false: replay; shards are compared, never written (a shard is
% written by a live call, after its node changes, and by repair_shards at
% the end of a complete replay). The store is full
% here (sa_record_call refilled it before the add); it is written or
% compared full and then kept light in memory (sa_store_io 'light').
keys = unique(ctx.dirty, 'stable');
ctx.dirty = {};
c = numel(ctx.attempts);
for k = 1:numel(keys)
    store = ctx.shards(keys{k});
    f = fullfile(ctx.state_dir, 'roots', sa_store_io('shard_file', keys{k}));
    if write
        sa_store_io('save_shard', store, f, keys{k}, c);
    elseif isKey(ctx.replay.shard_stamp, keys{k}) && ctx.replay.shard_stamp(keys{k}) == c
        Z = sa_store_io('read_shard', f);       % the decoded store
        if ~isequaln(Z.store, store)
            integrity(ctx, sprintf('stored shard %s (count %d) differs from the replay', keys{k}, c));
        end
        ctx.replay.shard_matched(keys{k}) = true;
    end
    ctx.shards(keys{k}) = sa_store_io('light', store);
end
end

function ctx = compare_checkpoints(ctx)
c = numel(ctx.attempts);
if ctx.replay.queue_stamp == c
    Z = load(fullfile(ctx.state_dir, 'queue.mat'));
    if ~isequaln(Z.Q, ctx.Q_at_last) || ~isequaln(Z.ring_calls, ctx.ring_calls_at_last)
        integrity(ctx, sprintf('queue.mat (count %d) differs from the replay state', c));
    end
    ctx.replay.queue_matched = true;
end
if ctx.replay.attempts_stamp == c
    Z = load(fullfile(ctx.state_dir, 'attempts.mat'));
    if ~isequaln(Z.attempts, ctx.attempts)
        integrity(ctx, sprintf('attempts.mat (count %d) differs from the replay state', c));
    end
    ctx.replay.attempts_matched = true;
end
end

function ctx = set_fp(ctx, fp)
if isfield(ctx, 'guess_fp'), ctx.guess_fp{numel(ctx.attempts)} = fp; end
end

function [ctx, rec] = replay_next(ctx, ukey, gn, fp)
rec = [];
P = ctx.replay;
while P.pos <= numel(P.files)
    try
        r = load(P.files{P.pos});
    catch ME
        integrity(ctx, sprintf('attempt file %s cannot be read: %s', P.files{P.pos}, ME.message));
    end
    P.pos = P.pos + 1;
    if isfield(r, 'ukey'), rk = r.ukey; else, rk = unit_key(r.unit); end
    if ~strcmp(rk, ukey)
        ctx.replay = P;
        integrity(ctx, sprintf(['stored attempt %s carries work unit "%s" but the scheduler ' ...
            'requested "%s"'], r.attempt.attempt_id, rk, ukey));
    end
    if strcmp(r.attempt.outcome, 'interrupted')
        if numel(ctx.attempts) + 1 ~= sscanf(r.attempt.attempt_id, 'A%d')
            integrity(ctx, 'attempt numbering gap at an interrupted record');
        end
        ctx.attempts(end+1) = r.attempt;
        ctx.guess_nodes(numel(ctx.attempts)) = gn;   % same work unit, same guess
        ctx = set_fp(ctx, fp);
        ctx.calls = numel(ctx.attempts);
        ctx.pending_retry = struct('id', r.attempt.attempt_id, 'key', ukey);
        continue
    end
    if numel(ctx.attempts) + 1 ~= sscanf(r.attempt.attempt_id, 'A%d')
        ctx.replay = P;
        integrity(ctx, sprintf('attempt numbering gap at %s', r.attempt.attempt_id));
    end
    ctx.pending_retry = [];
    rec = r;
    break
end
if P.pos > numel(P.files)
    P.active = false;
end
ctx.replay = P;
if ~P.active && isempty(rec)
    ctx = end_of_replay(ctx);       % exhausted on interrupted records only
end
end

function ctx = end_of_replay(ctx)
% Every stored checkpoint must have been met at its own attempt count; then
% a shard file left behind its node's last discovery is rewritten
% (repair_shards).
P = ctx.replay;
if P.queue_stamp > 0 && ~P.queue_matched
    integrity(ctx, sprintf('queue.mat stamped %d was never met by the replay', P.queue_stamp));
end
if P.attempts_stamp > 0 && ~P.attempts_matched
    integrity(ctx, sprintf('attempts.mat stamped %d was never met by the replay', P.attempts_stamp));
end
keys = P.shard_stamp.keys;
for k = 1:numel(keys)
    if ~isKey(P.shard_matched, keys{k})
        integrity(ctx, sprintf('shard %s stamped %d was never met by the replay', keys{k}, ...
            P.shard_stamp(keys{k})));
    end
end
ctx = repair_shards(ctx);
sa_journal('note', ctx.state_dir, sprintf('replay complete at attempt count %d; live from here', ...
    numel(ctx.attempts)));
end

function ctx = repair_shards(ctx)
% A kill after an accepted call's attempt file is written and before its
% node's shard is flushed leaves that shard file one discovery behind (or
% missing, for a node's first discovery); the replay rebuilt the node's
% store from the attempt files, but it compares shards only at their own
% stamps and never writes them. Here, after every stored shard has been met
% at its stamp (the checks above), each node whose shard file is missing
% or stamped below n_last, the attempt number of its last discovery, gets
% the file an uninterrupted run wrote at that call: the replayed store
% (which has no discovery after n_last) made full from the attempt files
% and saved stamped n_last. Two journal notes per file: an intent note
% before the save and a 'rewritten' note after it, so a stop inside the
% repair leaves an intent note without its 'rewritten' note (the next
% resume repairs again if the file is still behind, else finds nothing to
% do; a repeated save writes the same file). The in-memory store stays
% light. On a folder without such a gap nothing is written or noted. A shard
% stamped above its node's last discovery, or a shard of a node without a
% discovery, is never met by the replay and has already stopped the run as
% an integrity error above; it is never rewritten.
P = ctx.replay;
keys = ctx.shards.keys;
for k = 1:numel(keys)
    store = ctx.shards(keys{k});
    n_last = 0;
    for r = 1:numel(store.roots)
        ids = {store.roots(r).discoveries.attempt_id};
        for q = 1:numel(ids), n_last = max(n_last, sscanf(ids{q}, 'A%d')); end
    end
    if n_last == 0, continue; end
    if isKey(P.shard_stamp, keys{k})
        s = P.shard_stamp(keys{k});
        if s >= n_last, continue; end
        was = sprintf('stamped %d', s);
    else
        was = 'missing';
    end
    sd = struct('state_dir', ctx.state_dir);     % the handle captures the state folder only, not ctx
    fs = sa_store_io('full', store, ctx.M, ctx.S, ctx.K, @(id) get_sol(sd, id));
    f = fullfile(ctx.state_dir, 'roots', sa_store_io('shard_file', keys{k}));
    sa_journal('note', ctx.state_dir, sprintf(['shard %s repair starts: stamp %d at attempt count %d; ' ...
        'the file is %s, behind its last discovery A%06d (consistent with a stop between that call''s ' ...
        'attempt file and its shard write)'], keys{k}, n_last, numel(ctx.attempts), was, n_last));
    sa_store_io('save_shard', fs, f, keys{k}, n_last);
    sa_journal('note', ctx.state_dir, sprintf('shard %s rewritten stamped %d at attempt count %d (was %s)', ...
        keys{k}, n_last, numel(ctx.attempts), was));
end
end

% =========================================================================
% checkpoints and final outputs

function ctx = checkpoint(ctx)
% The checkpoint files (attempts.mat, process_time.mat, queue.mat; a
% resume compares them with its replay and never starts from them) every
% checkpoint.every_calls live calls, and the disk floor. The human-facing
% files (attempts.csv, summary.txt, roots_index.csv) are written on their
% own interval (reports_due) and at finish.
if ctx.replay.active, return; end
check_root_counts(ctx);
sa_report('write_state', ctx);
if free_gib(ctx.S.run_dir) < ctx.S.checkpoint.min_free_gib
    ctx = set_stop(ctx, 'disk_floor');
end
end

function tf = reports_due(ctx)
% Checked after every live call (like sa_progress('due'); no timer object):
% sa_report('due') holds the interval clock (report.files_interval_s since
% the files were last written, or since the run started). Its answer only
% decides when the human-facing files are rewritten; it never reaches the
% queue, a record or a key.
tf = ~ctx.replay.active && isfield(ctx, 't_reports') && sa_report('due', ctx);
end

function ctx = write_reports(ctx)
% attempts.csv, summary.txt, roots_index.csv from the render cache (ctx.rep)
[ctx.rep, ctx.t_reports] = sa_report('write_reports', ctx);
end

function check_root_counts(ctx)
% The node root-count vector must equal the shards' root counts (it steers
% the scheduler's scans; a difference would change the order of work)
if ~isfield(ctx, 'nroots'), return; end
n = zeros(size(ctx.nroots));
k = ctx.shards.keys;
for q = 1:numel(k), n(ctx.node_of(k{q})) = numel(ctx.shards(k{q}).roots); end
if ~isequal(n, ctx.nroots)
    integrity(ctx, 'the node root-count vector differs from the shards');
end
end

function ctx = finish(ctx)
if ctx.replay.active
    % stopped before the stored attempts were all replayed: keep every
    % stored derived file (a later resume replays the rest)
    if ~ctx.stop
        integrity(ctx, 'the scheduler finished with stored attempts left unreplayed');
    end
    return
end
check_root_counts(ctx);
ctx.pairs = sa_report('pairs', ctx);
sa_report('write_state', ctx);
% the human-facing files, with the self-check: the cached render must equal
% a full render from scratch, else the full one is written and the
% journal and the screen say so
[ctx.rep, ctx.t_reports] = sa_report('write_reports', ctx, true);
sa_report('final', ctx);
end

function g = free_gib(d)
g = NaN;
try
    g = double(java.io.File(d).getUsableSpace()) / 2^30;
catch
end
if ~(isfinite(g) && g > 0)
    [st, out] = system(sprintf('df -k "%s" | tail -1', d));
    if st == 0
        p = regexp(strtrim(out), '\s+', 'split');
        g = str2double(p{4}) / 2^20;
    end
end
end

function t = end_text(A)
% The journal 'end' text; a cold call cut by the cap gets ' cut_by_cap'
% (sa_resume's reconcile writes the same from the stored row)
t = A.outcome;
if isfield(A, 'cut_by_cap') && isequal(A.cut_by_cap, true), t = [t ' cut_by_cap']; end
end

function n = guess_count(guess)
% Node count of the guess the corrector receives
% (report-only leaf; NaN when the guess carries no mesh).
n = NaN;
if isstruct(guess) && isfield(guess, 'x'), n = numel(guess.x); end
end

function sol = get_sol(ctx, id)
% reads ctx.state_dir only (callers may pass struct('state_dir', d))
Z = load(fullfile(ctx.state_dir, 'attempts', [id '.mat']), 'sol');
sol = Z.sol;
end
