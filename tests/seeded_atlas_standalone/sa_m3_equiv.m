function D = sa_m3_equiv(dU, dR, opts)
%SA_M3_EQUIV  Strict resume equivalence (milestone-3 packet section 2e).
%
%   D = SA_M3_EQUIV(dirU, dirR) compares an uninterrupted run U with a run R
%   that was killed and resumed: the same completed work units in the same
%   order, bit-identical iterates (sol) and validity records, the same
%   attempt rows and the same store (per-node shards: memberships,
%   distances, per-node root ids). Permitted differences: 'interrupted'
%   records, attempt ids shifted by them (every A###### token in R is mapped
%   to U's numbering before comparing), resume_<n>.mat, and volatile leaves
%   (wall_s; times).
%   D = SA_M3_EQUIV(dirU, dirR, struct('prefix_only', true)) compares only
%   R's completed units against the first ones of U (a run R that stopped
%   earlier).
%   D.pass, D.n_units, D.n_interrupted, D.diffs, D.added (M5: attempt-row
%   fields present on one side only, e.g. nmax / cut_by_cap against an M4
%   run; never a difference).
%   M6 (packet v0.2 section 2.3): the row field max_proposed_mesh and the
%   attempt-file field mesh_trace are added leaves against an older run;
%   when BOTH sides carry mesh_trace it is compared too (isequaln without
%   the per-process call token): D.n_trace_compared.

if nargin < 3, opts = struct(); end
sa_report('refuse_if_compacted', dU, 'sa_m3_equiv');   % M7a RC-7.6: compaction-blind, refuses
sa_report('refuse_if_compacted', dR, 'sa_m3_equiv');
prefix = isfield(opts, 'prefix_only') && opts.prefix_only;
D = struct('pass', false, 'n_units', 0, 'n_interrupted', 0, 'diffs', {{}}, 'added', {{}}, 'n_trace_compared', 0);
[U, ~] = records(dU);
[R, nint] = records(dR);
D.n_interrupted = nint;
if prefix
    if numel(R) > numel(U), D.diffs{end+1} = 'R has more units than U'; end
    U = U(1:min(numel(U), numel(R)));
end
if numel(U) ~= numel(R)
    D.diffs{end+1} = sprintf('completed units %d vs %d', numel(U), numel(R));
end
% id map R -> U numbering
m = containers.Map('KeyType', 'char', 'ValueType', 'char');
for k = 1:numel(R), m(R(k).attempt.attempt_id) = sprintf('A%06d', k); end
n = min(numel(U), numel(R));
D.n_units = n;
for k = 1:n
    a = U(k); b = R(k);
    if ~strcmp(a.ukey, mapstr(b.ukey, m))
        D.diffs{end+1} = sprintf('unit %d differs: %s | %s', k, a.ukey, b.ukey);
        continue
    end
    if ~isequaln(a.sol, b.sol), D.diffs{end+1} = sprintf('unit %d: sol differs', k); end
    if ~isequaln(a.validity, b.validity), D.diffs{end+1} = sprintf('unit %d: validity differs', k); end
    if ~isequaln(a.opts, b.opts), D.diffs{end+1} = sprintf('unit %d: opts differ', k); end
    if isstruct(a.mesh_trace) && isstruct(b.mesh_trace)
        D.n_trace_compared = D.n_trace_compared + 1;
        if ~isequaln(notoken(a.mesh_trace), notoken(b.mesh_trace))
            D.diffs{end+1} = sprintf('unit %d: mesh trace differs', k);
        end
    elseif isstruct(a.mesh_trace) ~= isstruct(b.mesh_trace)
        D.added = union(D.added, {'mesh_trace'});
    end
    A = clean_row(a.attempt, []); B = clean_row(b.attempt, m);
    % M5 (packet v0.2 section 6): row fields present on one side only (nmax,
    % cut_by_cap against an M2-M4 run) are added leaves, not differences
    only = setxor(fieldnames(A), fieldnames(B));
    if ~isempty(only)
        D.added = union(D.added, only);
        A = rmfield(A, intersect(only, fieldnames(A))); B = rmfield(B, intersect(only, fieldnames(B)));
        B = orderfields(B, A);
    end
    if ~isequaln(A, B)
        f = fieldnames(A); bad = f(cellfun(@(x) ~isequaln(A.(x), B.(x)), f));
        D.diffs{end+1} = sprintf('unit %d: attempt row differs in %s', k, strjoin(bad.', ','));
    end
end
if ~prefix
    fu = dir(fullfile(sa_layout('state', dU), 'roots', '*.mat')); fr = dir(fullfile(sa_layout('state', dR), 'roots', '*.mat'));
    if ~isequal(sort({fu.name}), sort({fr.name}))
        D.diffs{end+1} = 'shard file sets differ';
    else
        for k = 1:numel(fu)
            SU = sa_store_io('read_shard', fullfile(fu(k).folder, fu(k).name));   % AM-9 decoded
            SR = sa_store_io('read_shard', fullfile(sa_layout('state', dR), 'roots', fu(k).name));
            if ~isequaln(SU.store, map_store(SR.store, m))
                D.diffs{end+1} = sprintf('shard %s differs', fu(k).name);
            end
        end
    end
end
D.pass = isempty(D.diffs);
end

function [L, nint] = records(d)
f = dir(fullfile(sa_layout('state', d), 'attempts', 'A*.mat'));
names = sort({f.name});
names = names(~cellfun(@isempty, regexp(names, '^A\d{6}\.mat$', 'once')));
L = struct('attempt', {}, 'sol', {}, 'validity', {}, 'opts', {}, 'ukey', {}, 'mesh_trace', {});
nint = 0;
for k = 1:numel(names)
    r = load(fullfile(sa_layout('state', d), 'attempts', names{k}));
    if strcmp(r.attempt.outcome, 'interrupted'), nint = nint + 1; continue; end
    tr = [];
    if isfield(r, 'mesh_trace'), tr = r.mesh_trace; end
    L(end+1) = struct('attempt', r.attempt, 'sol', r.sol, 'validity', r.validity, ...
        'opts', r.opts, 'ukey', sa_search('unit_key', r.unit), 'mesh_trace', tr); %#ok<AGROW>
end
end

function A = clean_row(A, m)
A.wall_s = 0; A.retry_of = '';
if ~isempty(m)
    for f = {'attempt_id', 'file', 'parent_attempt', 'guess_source', 'unit_key', 'root_id'}
        A.(f{1}) = mapstr(A.(f{1}), m);
    end
end
end

function s = mapstr(s, m)
[tok, sp] = regexp(s, 'A\d{6}', 'match', 'split');
for k = 1:numel(tok)
    if isKey(m, tok{k}), tok{k} = m(tok{k}); end
end
c = [sp; [tok, {''}]];
s = [c{:}];
end

function st = map_store(st, m)
for r = 1:numel(st.roots)
    st.roots(r).representative_attempt_id = mapstr(st.roots(r).representative_attempt_id, m);
    for q = 1:numel(st.roots(r).discoveries)
        st.roots(r).discoveries(q).attempt_id = mapstr(st.roots(r).discoveries(q).attempt_id, m);
    end
end
end

function t = notoken(t)
if isfield(t, 'token'), t = rmfield(t, 'token'); end
end
