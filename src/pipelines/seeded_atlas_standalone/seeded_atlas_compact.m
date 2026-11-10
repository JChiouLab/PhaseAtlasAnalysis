function out = seeded_atlas_compact(run_dir, opts)
%SEEDED_ATLAS_COMPACT  Optional clean-up of a completed or cleanly stopped run
%   (default: never).
%
%   out = SEEDED_ATLAS_COMPACT(run_dir)
%   out = SEEDED_ATLAS_COMPACT(run_dir, opts)
%
%   1. Refused unless the run's status (sa_resume 'status' against its own
%      stored settings) is completed and not yet compacted; a flat (software
%      0.5) folder is evidence and is refused. A cleanly STOPPED run (not
%      completed: queue_final.mat present, its attempt count equal to the
%      attempt files, no journal start without an end, no integrity stop;
%      stop reasons such as disk_floor, budget_calls, budget_wall) is
%      compacted only with opts.stopped_ok = true, and a line then says
%      that it can no longer be resumed; without the opt-in the refusal
%      says how to opt in. A killed run (no final queue, a count mismatch
%      or a call started without an end) stays refused even with it.
%      compaction.mat / .txt record compacted_from 'stopped' and the stop
%      reason for such a run (fields present only then). A compacted
%      stopped run remains an import source (sa_import reads
%      compaction.mat for its attempt count).
%   2. The verified root load seeded_atlas_roots(run_dir, struct('verify',
%      true)) runs FIRST; any failure stops the compaction with every file
%      in place.
%   3. Only then: internal/compaction.mat / compaction.txt list every file
%      to be removed (name, bytes, SHA-256) and the attempt count before
%      compaction, and the files are removed: the attempt files of
%      non-discovery calls (invalid, no-return, substeps, refinement,
%      interrupted), queue.mat and process_time.mat. Kept: everything
%      verified loading and traces need (settings_effective.mat,
%      attempts.mat, roots/, every discovery's attempt file), plus
%      settings_user.mat, journal.txt, queue_final.mat, resume_*.mat,
%      identity_pairs.csv, pilot/*.csv, every top-level file, plots/ and
%      traces/.
%   4. The verified load runs again and must equal the first (isequaln).
%   Lost afterwards: resume (SA:RESUME:E_COMPACTED), re-checks of failed
%   calls, the exact cut recount, rebuilding the pilot report, equivalence
%   comparisons, inspecting failed iterates; those tools refuse on
%   internal/compaction.mat. seeded_atlas_roots, seeded_atlas_plot,
%   seeded_atlas_trace and seeded_atlas_compare keep working.
%   opts.stopped_ok: true allows a cleanly stopped run (above); default false.
%   opts.fail_load (test hook): an error thrown in place of the first
%   verified load.
%   out: n_removed, bytes_removed, attempt_count_before, n_roots, manifest.

if nargin < 2, opts = struct(); end
if strcmp(sa_layout('layout', run_dir), 'flat')
    error('SA:COMPACT:E_LEGACY', '%s is a flat (software 0.5 or earlier) folder: evidence, never compacted.', ...
        run_dir);
end
D = sa_layout('state', run_dir);
if exist(fullfile(D, 'compaction.mat'), 'file') == 2
    error('SA:COMPACT:E_DONE', '%s is already compacted.', run_dir);
end
if exist(fullfile(D, 'settings_effective.mat'), 'file') ~= 2
    error('SA:COMPACT:E_NO_RUN', '%s holds no run.', run_dir);
end
E = load(fullfile(D, 'settings_effective.mat'), 'settings');
% the stored settings, with the focus-mode leaves a folder made before them
% lacks at their neutral values (sa_resume compares filled structs)
Z = sa_resume('status', run_dir, sa_settings('fill_focus_leaves', E.settings), struct());
stopped = '';                                     % the stop reason of a cleanly stopped run
if ~strcmp(Z.state, 'completed')
    [clean, why] = clean_stop(D, Z);
    if ~clean
        error('SA:COMPACT:E_NOT_COMPLETED', ['compaction refused: %s is %s (%s); only a completed run, or a cleanly ' ...
            'stopped one with opts.stopped_ok = true, is compacted; this one is not cleanly stopped (%s): resume it ' ...
            'to a stop first.'], run_dir, Z.state, Z.message, why);
    end
    if ~(isfield(opts, 'stopped_ok') && isequal(opts.stopped_ok, true))
        error('SA:COMPACT:E_NOT_COMPLETED', ['compaction refused: %s stopped cleanly (%s) but is not completed (%s). ' ...
            'To compact it anyway, call seeded_atlas_compact(run_dir, struct(''stopped_ok'', true)); the run can ' ...
            'then no longer be resumed (its files of calls that found no root are removed).'], run_dir, ...
            Z.stop_reason, Z.state);
    end
    stopped = Z.stop_reason;
    fprintf(['seeded_atlas_compact: %s stopped (%s) without completing; after this compaction it can NO LONGER ' ...
        'be resumed (opts.stopped_ok).\n'], run_dir, stopped);
end
% ---- 2. verified load first -----------------------------------------------------
if isfield(opts, 'fail_load') && ~isempty(opts.fail_load)
    error('SA:COMPACT:E_LOAD', 'compaction stopped before any removal: %s', opts.fail_load);
end
R0 = seeded_atlas_roots(run_dir, struct('verify', true));
% ---- 3. the manifest, then the removal ------------------------------------------------
keep = containers.Map('KeyType', 'char', 'ValueType', 'logical');
sh = sa_store_io('shard_files', D);              % skips sa_write temporary files
for k = 1:numel(sh)
    S = sa_store_io('read_shard', fullfile(sh(k).folder, sh(k).name));
    for r = 1:numel(S.store.roots)
        keep(S.store.roots(r).representative_attempt_id) = true;
        for q = 1:numel(S.store.roots(r).discoveries)
            keep(S.store.roots(r).discoveries(q).attempt_id) = true;
        end
    end
end
fa = dir(fullfile(D, 'attempts', 'A*.mat'));
fa = fa(~cellfun(@isempty, regexp({fa.name}, '^A\d{6}\.mat$', 'once')));
n_before = numel(fa);
rm = {};
for k = 1:numel(fa)
    id = fa(k).name(1:end-4);
    if ~isKey(keep, id), rm{end+1} = fullfile('attempts', fa(k).name); end %#ok<AGROW>
end
for f = {'queue.mat', 'process_time.mat'}
    if exist(fullfile(D, f{1}), 'file') == 2, rm{end+1} = f{1}; end %#ok<AGROW>
end
files = struct('file', {}, 'bytes', {}, 'sha256', {});
for k = 1:numel(rm)
    p = fullfile(D, rm{k});
    i = dir(p);
    files(end+1) = struct('file', rm{k}, 'bytes', i.bytes, 'sha256', sa_sha256(p)); %#ok<AGROW>
end
M = struct('compacted_at', datestr(now, 'yyyy-mm-ddTHH:MM:SS'), 'n_removed', numel(files), ... %#ok<DATST,TNOW1>
    'bytes_removed', sum([files.bytes]), 'attempt_count_before', n_before, 'n_discovery_files', keep.Count, ...
    'removed', files, 'software_version', sa_settings('version'), 'done', false);
if ~isempty(stopped)
    M.compacted_from = 'stopped';                 % not completed (opts.stopped_ok)
    M.stop_reason = stopped;
end
sa_write('mat', fullfile(D, 'compaction.mat'), M);
sa_write('text', fullfile(D, 'compaction.txt'), manifest_text(M, run_dir));
for k = 1:numel(rm)
    delete(fullfile(D, rm{k}));
end
M.done = true;
sa_write('mat', fullfile(D, 'compaction.mat'), M);
% ---- 4. the verified load again --------------------------------------------------------
R1 = seeded_atlas_roots(run_dir, struct('verify', true));
if ~isequaln(R0, R1)
    error('SA:COMPACT:E_AFTER', ['compaction of %s: the verified root load differs after the removal ' ...
        '(the removed files are listed with SHA-256 in compaction.txt).'], run_dir);
end
out = struct('n_removed', M.n_removed, 'bytes_removed', M.bytes_removed, ...
    'attempt_count_before', n_before, 'n_roots', numel(R1), 'manifest', fullfile(D, 'compaction.txt'));
fprintf('compacted %s: %d files removed (%.1f MB); %d roots verified before and after\n', run_dir, ...
    M.n_removed, M.bytes_removed / 2^20, numel(R1));
end

function [ok, why] = clean_stop(D, Z)
% A cleanly stopped run: not refused (an integrity stop, a settings change),
% queue_final.mat present with the attempt-file count, every journal start
% ended (a killed run fails one of these)
ok = false; why = '';
if ~any(strcmp(Z.state, {'resumable', 'budget_exhausted'}))
    why = sprintf('status %s', Z.state); return
end
qf = fullfile(D, 'queue_final.mat');
if exist(qf, 'file') ~= 2
    why = 'no queue_final.mat: killed or still running'; return
end
Q = load(qf, 'attempt_count');
n = Z.attempt_count;
if Q.attempt_count ~= n
    why = sprintf('queue_final.mat counts %d calls but %d attempt files exist: killed or still running', ...
        Q.attempt_count, n); return
end
J = sa_journal('read', D);
open_ = setdiff({J.starts.id}, {J.ends.id});
if ~isempty(open_)
    why = sprintf('the journal has a call started without an end (%s): killed or still running', open_{1}); return
end
ok = true;
end

function t = manifest_text(M, run_dir)
L = {sprintf('compaction of %s on %s (%s)', run_dir, M.compacted_at, M.software_version), ...
    sprintf(['attempt files before: %d; discovery files kept: %d; files removed: %d (%d bytes); the run can no ' ...
    'longer be resumed'], M.attempt_count_before, M.n_discovery_files, M.n_removed, M.bytes_removed), ...
    'file,bytes,sha256'};
if isfield(M, 'compacted_from')
    L = [L(1:2), {sprintf(['compacted from a STOPPED run (stop reason %s), not a completed one ' ...
        '(opts.stopped_ok): it can no longer be resumed'], M.stop_reason)}, L(3)];
end
for k = 1:numel(M.removed)
    L{end+1} = sprintf('%s,%d,%s', M.removed(k).file, M.removed(k).bytes, M.removed(k).sha256); %#ok<AGROW>
end
t = [strjoin(L, newline) newline];
end
