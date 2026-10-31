function R = sag_compare_runs(mode, dirA, dirB)
%SAG_COMPARE_RUNS  Milestone-2C regression comparator (packet 2C section 1a).
%
%   R = SAG_COMPARE_RUNS('leaves', dirA, dirB)
%       RC-6 leaf list, for a baseline against the stored pre-repair M2 smoke:
%       every attempts/A####.mat: sol.x, sol.y, sol.yp, sol.stats, opts,
%       warning_id, and the value and pass flag of every check (matched by
%       name; the check-name sets must agree); the roots in roots.mat
%       (profiles and memberships); identity_pairs.csv; summary.txt with the
%       volatile run_dir and wall entries masked. Record STRUCTURE is not
%       compared (RC-6).
%   R = SAG_COMPARE_RUNS('whole', dirA, dirB)
%       Every leaf of every output file (settings_effective.mat, attempts.mat,
%       attempts/*.mat, roots.mat, attempts.csv, identity_pairs.csv,
%       summary.txt), skipping ONLY the frozen volatile list below.
%   L = SAG_COMPARE_RUNS('volatile')
%       The frozen volatile list (frozen from code reading before any
%       comparison was made):
%         run_dir and file paths  settings.run_dir and any text leaf that
%                                 contains the run folder (masked to <RUN>)
%         started_at              settings_effective.started_at
%         resolved_at             settings.resolved_at (sa_settings)
%         wall_s                  per-call corrector wall time (attempt
%                                 records, attempts.csv column) and the
%                                 summary's "wall %.1f s"
%         path_entries, pwd       settings_effective.runtime
%         software_version        only if bumped (it is compared here)
%
%   R.n_diff, R.diffs (cellstr 'path: detail'), R.n_compared, R.volatile_seen.
%   Read only: nothing is written.
%   Milestone 3 (additive): a run folder without roots.mat is read from its
%   per-node shards roots/*.mat (roots concatenated in shard-file order), and
%   'whole' also walks queue.mat and every shard file when present.

switch mode
    case 'volatile', R = volatile(); return
    case 'leaves',   R = leaves(dirA, dirB);
    case 'whole',    R = whole(dirA, dirB);
    otherwise, error('SAG:COMPARE', 'Unknown mode %s.', mode);
end
R.pass = R.n_diff == 0;
end

% =========================================================================

function L = volatile()
L = {'run_dir', 'started_at', 'resolved_at', 'wall_s', 'path_entries', 'pwd'};
end

function R = blank()
R = struct('n_diff', 0, 'diffs', {{}}, 'n_compared', 0, 'volatile_seen', {{}}, 'pass', false);
end

function R = note(R, path, detail)
R.n_diff = R.n_diff + 1;
R.diffs{end+1} = sprintf('%s: %s', path, detail);
end

% =========================================================================
% RC-6 leaf list

function R = leaves(dA, dB)
R = blank();
fa = dir(fullfile(sa_layout('state', dA), 'attempts', 'A*.mat'));
fb = dir(fullfile(sa_layout('state', dB), 'attempts', 'A*.mat'));
na = sort({fa.name}); nb = sort({fb.name});
if ~isequal(na, nb)
    R = note(R, 'attempts/', sprintf('file sets differ (%d vs %d)', numel(na), numel(nb)));
end
for k = 1:numel(na)
    if ~any(strcmp(na{k}, nb)), continue; end
    A = load(fullfile(sa_layout('state', dA), 'attempts', na{k}));
    B = load(fullfile(sa_layout('state', dB), 'attempts', na{k}));
    p = ['attempts/' na{k}];
    for f = {'x', 'y', 'yp', 'stats'}
        [R] = cmp_leaf(R, [p ':sol.' f{1}], getsol(A.sol, f{1}), getsol(B.sol, f{1}));
    end
    R = cmp_leaf(R, [p ':opts'], A.opts, B.opts);
    R = cmp_leaf(R, [p ':warning_id'], A.warning_id, B.warning_id);
    ca = A.validity.checks; cb = B.validity.checks;
    if ~isequal(sort({ca.name}), sort({cb.name}))
        R = note(R, [p ':check names'], sprintf('{%s} vs {%s}', strjoin({ca.name}, ','), ...
            strjoin({cb.name}, ',')));
    end
    for q = 1:numel(ca)
        h = find(strcmp(ca(q).name, {cb.name}));
        if numel(h) ~= 1, continue; end
        R = cmp_leaf(R, sprintf('%s:check %s.value', p, ca(q).name), ca(q).value, cb(h).value);
        R = cmp_leaf(R, sprintf('%s:check %s.pass', p, ca(q).name), ca(q).pass, cb(h).pass);
    end
    R = cmp_leaf(R, [p ':validity.valid'], A.validity.valid, B.validity.valid);
end
% roots: profiles and memberships
ra = run_roots(dA); rb = run_roots(dB);
R = cmp_leaf(R, 'roots.mat:n_roots', numel(ra), numel(rb));
for r = 1:min(numel(ra), numel(rb))
    p = sprintf('roots.mat:roots(%d)', r);
    for f = {'root_id', 'target_key', 'coordinate', 'comparison_mesh', ...
             'representative_profile', 'representative_attempt_id'}
        R = cmp_leaf(R, [p '.' f{1}], ra(r).(f{1}), rb(r).(f{1}));
    end
    da = ra(r).discoveries; db = rb(r).discoveries;
    R = cmp_leaf(R, [p '.n_discoveries'], numel(da), numel(db));
    for q = 1:min(numel(da), numel(db))
        for f = {'discovery_number', 'attempt_id', 'coordinate', 'seed_family', ...
                 'observed_family', 'route', 'returned_mesh', 'returned_profile', ...
                 'comparison_profile', 'member_root_id', 'representative_distance', 'tolerance'}
            R = cmp_leaf(R, sprintf('%s.discoveries(%d).%s', p, q, f{1}), da(q).(f{1}), db(q).(f{1}));
        end
    end
end
R = cmp_leaf(R, 'identity_pairs.csv', fileread(fullfile(sa_layout('state', dA), 'identity_pairs.csv')), ...
    fileread(fullfile(sa_layout('state', dB), 'identity_pairs.csv')));
% summary.txt COUNTS (RC-6): every line, masked, except the identity-pair
% wording sentences under "Identity pairs:" (rewritten by the D5/F2 repair
% ba5aee0 after the stored smoke); the "Identity pairs: N ...; M in the band"
% count line itself is kept.
R = cmp_leaf(R, 'summary.txt (counts, masked)', ...
    summary_counts(mask_text(fileread(fullfile(dA, 'summary.txt')), dA)), ...
    summary_counts(mask_text(fileread(fullfile(dB, 'summary.txt')), dB)));
end

function r = run_roots(d)
if exist(fullfile(d, 'roots.mat'), 'file') == 2
    S = load(fullfile(d, 'roots.mat')); r = S.store.roots; return
end
f = dir(fullfile(sa_layout('state', d), 'roots', '*.mat'));
r = [];
for k = 1:numel(f)
    S = load(fullfile(f(k).folder, f(k).name));
    if isfield(S, 'format') && strcmp(S.format, 'sa_store_io/shard-2')   % AM-9 decoded
        S = sa_store_io('read_shard', fullfile(f(k).folder, f(k).name));
    end
    r = [r, S.store.roots]; %#ok<AGROW>
end
end

function t = summary_counts(t)
L = strsplit(t, newline);
keep = true(size(L));
k = find(startsWith(L, 'Identity pairs:'), 1);
if ~isempty(k)
    q = k + 1;
    while q <= numel(L) && startsWith(L{q}, '  ')
        keep(q) = false; q = q + 1;
    end
end
t = strjoin(L(keep), newline);
end

function v = getsol(s, f)
if isstruct(s) && isfield(s, f), v = s.(f); else, v = '<absent>'; end
end

function R = cmp_leaf(R, path, a, b)
R.n_compared = R.n_compared + 1;
if ~isequaln(a, b) || ~strcmp(class(a), class(b))
    R = note(R, path, describe(a, b));
end
end

function s = describe(a, b)
if isnumeric(a) && isnumeric(b) && isequal(size(a), size(b)) && ~isempty(a)
    s = sprintf('max |diff| %.3g', max(abs(double(a(:)) - double(b(:)))));
elseif ischar(a) && ischar(b)
    s = sprintf('text differs (%d vs %d chars)', numel(a), numel(b));
else
    s = sprintf('%s %s vs %s %s', class(a), mat2str(size(a)), class(b), mat2str(size(b)));
end
end

function t = mask_text(t, d)
t = strrep(t, d, '<RUN>');
t = regexprep(t, 'wall [0-9.]+ s', 'wall <W> s');
end

% =========================================================================
% whole-record comparison

function R = whole(dA, dB)
R = blank();
mats = {'settings_effective.mat', 'attempts.mat'};
for f = {'roots.mat', 'queue.mat'}
    if exist(fullfile(sa_layout('state', dA), f{1}), 'file') == 2 || exist(fullfile(sa_layout('state', dB), f{1}), 'file') == 2
        mats{end+1} = f{1}; %#ok<AGROW>
    end
end
ra = dir(fullfile(sa_layout('state', dA), 'roots', '*.mat')); rb = dir(fullfile(sa_layout('state', dB), 'roots', '*.mat'));
if ~isequal(sort({ra.name}), sort({rb.name}))
    R = note(R, 'roots/', 'shard file sets differ');
end
for k = 1:numel(ra), mats{end+1} = fullfile('roots', ra(k).name); end %#ok<AGROW>
fa = dir(fullfile(sa_layout('state', dA), 'attempts', 'A*.mat'));
fb = dir(fullfile(sa_layout('state', dB), 'attempts', 'A*.mat'));
if ~isequal(sort({fa.name}), sort({fb.name}))
    R = note(R, 'attempts/', 'file sets differ');
end
for k = 1:numel(fa), mats{end+1} = fullfile('attempts', fa(k).name); end %#ok<AGROW>
ctx = struct('dA', dA, 'dB', dB, 'vol', {volatile()});
for k = 1:numel(mats)
    A = load(fullfile(sa_layout('state', dA), mats{k})); B = load(fullfile(sa_layout('state', dB), mats{k}));
    R = walk(R, A, B, mats{k}, ctx);
end
R = cmp_leaf(R, 'identity_pairs.csv', mask_text(fileread(fullfile(sa_layout('state', dA), 'identity_pairs.csv')), dA), ...
    mask_text(fileread(fullfile(sa_layout('state', dB), 'identity_pairs.csv')), dB));
R = cmp_leaf(R, 'summary.txt', mask_text(fileread(fullfile(dA, 'summary.txt')), dA), ...
    mask_text(fileread(fullfile(dB, 'summary.txt')), dB));
R = cmp_leaf(R, 'attempts.csv (wall_s column dropped)', ...
    drop_col(fileread(fullfile(dA, 'attempts.csv')), 'wall_s'), ...
    drop_col(fileread(fullfile(dB, 'attempts.csv')), 'wall_s'));
end

function t = drop_col(t, name)
L = strsplit(strtrim(t), newline);
hdr = strsplit(L{1}, ',');
c = find(strcmp(hdr, name));
for k = 1:numel(L)
    parts = regexp(L{k}, ',(?=(?:[^"]*"[^"]*")*[^"]*$)', 'split');
    parts(c) = [];
    L{k} = strjoin(parts, ',');
end
t = strjoin(L, newline);
end

function R = walk(R, a, b, path, ctx)
last = regexp(path, '[^.:(]+$', 'match', 'once');
if any(strcmp(last, ctx.vol))
    R.volatile_seen{end+1} = path;
    return
end
if ~strcmp(class(a), class(b))
    R = note(R, path, sprintf('class %s vs %s', class(a), class(b))); return
end
if isstruct(a)
    if ~isequal(size(a), size(b))
        R = note(R, path, sprintf('size %s vs %s', mat2str(size(a)), mat2str(size(b)))); return
    end
    fa = fieldnames(a); fb = fieldnames(b);
    if ~isequal(sort(fa), sort(fb))
        R = note(R, path, sprintf('fields differ: {%s} vs {%s}', strjoin(setdiff(fa, fb)', ','), ...
            strjoin(setdiff(fb, fa)', ',')));
    end
    common = intersect(fa, fb, 'stable');
    for e = 1:numel(a)
        if numel(a) == 1, pe = path; else, pe = sprintf('%s(%d)', path, e); end
        for f = 1:numel(common)
            R = walk(R, a(e).(common{f}), b(e).(common{f}), [pe '.' common{f}], ctx);
        end
    end
elseif iscell(a)
    if ~isequal(size(a), size(b))
        R = note(R, path, sprintf('size %s vs %s', mat2str(size(a)), mat2str(size(b)))); return
    end
    for e = 1:numel(a)
        R = walk(R, a{e}, b{e}, sprintf('%s{%d}', path, e), ctx);
    end
elseif ischar(a)
    R = cmp_leaf(R, path, strrep(a, ctx.dA, '<RUN>'), strrep(b, ctx.dB, '<RUN>'));
elseif isa(a, 'function_handle')
    R = cmp_leaf(R, path, func2str(a), func2str(b));
else
    R = cmp_leaf(R, path, a, b);
end
end
