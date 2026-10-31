function out = seeded_atlas_compact(run_dir, opts)
%SEEDED_ATLAS_COMPACT  Optional clean-up of a completed run (default: never).
%
%   out = SEEDED_ATLAS_COMPACT(run_dir)
%   out = SEEDED_ATLAS_COMPACT(run_dir, opts)
%
%   1. Refused unless the run's status (sa_resume 'status' against its own
%      stored settings) is completed and not yet compacted; a flat (software
%      0.5) folder is evidence and is refused.
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
Z = sa_resume('status', run_dir, E.settings, struct());
if ~strcmp(Z.state, 'completed')
    error('SA:COMPACT:E_NOT_COMPLETED', 'compaction refused: %s is %s (%s); only a completed run is compacted.', ...
        run_dir, Z.state, Z.message);
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

function t = manifest_text(M, run_dir)
L = {sprintf('compaction of %s on %s (%s)', run_dir, M.compacted_at, M.software_version), ...
    sprintf(['attempt files before: %d; discovery files kept: %d; files removed: %d (%d bytes); the run can no ' ...
    'longer be resumed'], M.attempt_count_before, M.n_discovery_files, M.n_removed, M.bytes_removed), ...
    'file,bytes,sha256'};
for k = 1:numel(M.removed)
    L{end+1} = sprintf('%s,%d,%s', M.removed(k).file, M.removed(k).bytes, M.removed(k).sha256); %#ok<AGROW>
end
t = [strjoin(L, newline) newline];
end
