function rows = sa_t48_am10(E)
%SA_T48_AM10  Milestone-4 T48 (packet v0.3 AM-10, RC-4.2): R3 resumes the
%   refused edge. Stub with the M3 refusal pattern: every call within 0.4 of
%   a lattice step of the wall node q* = [2 0] (coordinate [1.3 0]) is
%   refused, so every continuation, R1 and R2 attempt at q* fails and R3
%   runs. Asserts: R3 starts at the resumed edge's last valid fraction with
%   step s_fail/2 from that substep's stored iterate, is refused at its third
%   refusal, records resumed_edge / start_fraction / start_step, keeps
%   h = max_halvings + 2; no R3 call repeats an earlier call (guess
%   fingerprint = guess bits + coordinate + options) and none coincides with
%   an R1 call; the bound (<= 10 calls at h = 4) by exhaustive enumeration
%   of the step rule (validated against the realised edges); the fallback
%   (no refused cardinal edge into q: the v0.4 rule from the parent); a kill
%   inside R3, resumed, equals the uninterrupted run.
T = 'T48 AM-10 R3 resumes the refused edge';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't48');
qs = [1.3 0];
fail = @(dd, b, ex) wall(dd, qs);
box = [3 1];
U = fullfile(d, 'U');
r = sa_m3_run(U, box, 'fail', fail);
Q = r.Q; A = r.attempts;
Z = load(fullfile(sa_layout('state', U), 'attempts.mat'));
fp = Z.guess_fingerprint;
q = find(all(Q.index == [2 0], 2));
e3 = Q.edges(strcmp({Q.edges.kind}, 'rescue_deep') & [Q.edges.q] == q);
cand = Q.edges([Q.edges.q] == q & strcmp({Q.edges.kind}, 'continuation') & strcmp({Q.edges.status}, 'refused'));
ok = numel(e3) == 1;
det = '';
if ok
    er = Q.edges(e3.resumed_edge);
    fr = max([cand.fraction_valid]);
    last = A(strcmp({A.attempt_id}, er.attempts{end}));
    sfail = last.fraction - er.fraction_valid;
    a3 = A(ismember({A.attempt_id}, e3.attempts));
    sv = A(ismember({A.attempt_id}, er.attempts) & [A.valid] & [A.fraction] == er.fraction_valid);
    ok = ok && strcmp(er.kind, 'continuation') && er.fraction_valid == fr && e3.start_fraction == fr && ...
        e3.start_step == sfail / 2 && e3.h == 6 && strcmp(e3.status, 'refused') && e3.refusals == 3 && ...
        abs(a3(1).fraction - (fr + sfail / 2)) < 1e-15 && strcmp(a3(1).parent_attempt, sv(end).attempt_id) && ...
        e3.calls <= 10;
    det = sprintf(['resumed edge %d (fraction %.6g, s_fail %.6g); R3 calls %d at fractions %s, refusals %d, ' ...
        'status %s'], er.id, fr, sfail, e3.calls, mat2str([a3.fraction], 6), e3.refusals, e3.status);
end
rows = sa_check(rows, T, ['R3 resumes the refused cardinal edge with the largest last valid fraction: pos = fraction_valid, ' ...
    's = s_fail/2, guess = that substep''s iterate; refused at the third refusal; record fields; h kept 6'], ok, det);
% ---- no repeated call ---------------------------------------------------------
i3 = find(strcmp({A.route}, 'rescue_deep'));
i1 = find(strcmp({A.route}, 'rescue_diagonal'));
rep = 0; rep1 = 0;
for k = i3
    rep = rep + any(strcmp(fp(1:k-1), fp{k}));
    rep1 = rep1 + any(strcmp(fp(i1), fp{k}));
end
allfp = fp(~cellfun(@isempty, fp));
rows = sa_check(rows, T, 'no R3 call repeats an earlier call (guess, coordinate, options) and none coincides with an R1 call', ...
    ~isempty(i3) && rep == 0 && rep1 == 0 && numel(allfp) == numel(A), ...
    sprintf('%d R3 calls; repeats %d; R1 coincidences %d; run-wide repeated fingerprints %d', numel(i3), rep, ...
    rep1, numel(allfp) - numel(unique(allfp))));
% ---- the bound by exhaustive enumeration --------------------------------------------
agree = true;
for e = [cand, e3]
    a = A(ismember({A.attempt_id}, e.attempts));
    pat = arrayfun(@(x) x.valid, a);
    if strcmp(e.kind, 'rescue_deep')
        [n, fr] = sim(e.start_fraction, e.start_step, 2, pat);
    else
        [n, fr] = sim(0, 1, 4, pat);
    end
    agree = agree && n == numel(a) && max(abs(fr - [a.fraction])) == 0;
end
ends = edge_ends(4);                      % every refused v0.4 end state (pos, s_fail), h = 4
worst = 0; nstates = 0;
for k = 1:size(ends, 1)
    pos = ends(k, 1); s = ends(k, 2) / 2;
    worst = max(worst, worst_calls(pos, s, 2));
    nstates = nstates + 1;
end
rows = sa_check(rows, T, ['bound: over every refused h = 4 end state and every outcome sequence R3 makes <= 10 calls ' ...
    '(the step rule simulator reproduces the realised edges)'], agree && worst <= 10, ...
    sprintf('%d end states; worst %d calls', nstates, worst));
% ---- fallback: no refused cardinal edge into q -------------------------------------------
o = struct('targets', struct('list', [1 0; 3 0]), 'routes', struct('continuation', false), ...
    'rescue', struct('rungs', {{'deep'}}));
rF = sa_m3_run(fullfile(d, 'fallback'), box, 'fail', fail, 'settings', o);
eF = rF.Q.edges(strcmp({rF.Q.edges.kind}, 'rescue_deep') & [rF.Q.edges.q] == q);
aF = rF.attempts(ismember({rF.attempts.attempt_id}, eF.attempts));
rows = sa_check(rows, T, 'fallback without a refused edge: v0.4 rule (from the parent at the full step, h = 6, <= 13 calls)', ...
    numel(eF) == 1 && eF.resumed_edge == 0 && eF.start_fraction == 0 && eF.start_step == 1 && ...
    aF(1).fraction == 1 && eF.calls <= 13 && strcmp(eF.status, 'refused') && eF.refusals == 7, ...
    sprintf('calls %d refusals %d', eF.calls, eF.refusals));
% ---- kill inside R3, resume -----------------------------------------------------------------
k = i3(2);
K = fullfile(d, 'K');
try
    sa_m3_run(K, box, 'fail', fail, 'kill_at', k);
catch ME
    assert(strcmp(ME.identifier, 'SA:TEST:KILL'), ME.message);
end
rK = sa_m3_run(K, box, 'fail', fail, 'settings', struct('resume', true));
D = sa_m3_equiv(U, K);
rows = sa_check(rows, T, 'kill inside R3 (its second call), resume == uninterrupted', ...
    D.pass && isequal(rK.replay_interrupted, k), strjoin(D.diffs, ' | '));
rmpath(E.stub);
end

function c = wall(dd, qs)
c = 'S';
if max(abs(dd - qs)) < 0.26 - 1e-12, c = 'R'; end
end

function [n, fr] = sim(pos, s, h, pat)
% the adaptive step rule (sa_continuation run_edge) on an outcome pattern
n = 0; r = 0; fr = [];
while true
    f = pos + s; n = n + 1; fr(end+1) = f; %#ok<AGROW>
    if n > numel(pat), ok = false; else, ok = pat(n); end
    if ok
        pos = f;
        if pos == 1, return; end
        s = min(2 * s, 1 - pos);
    else
        s = s / 2; r = r + 1;
        if r > h, return; end
    end
end
end

function E = edge_ends(h)
% every (pos, s_fail) at which a v0.4 adaptive edge (from pos 0, s 1) is refused
E = zeros(0, 2);
stack = {[0 1 0]};
while ~isempty(stack)
    st = stack{end}; stack(end) = [];
    pos = st(1); s = st(2); r = st(3);
    f = pos + s;
    % success
    if f ~= 1
        stack{end+1} = [f min(2 * s, 1 - f) r]; %#ok<AGROW>
    end
    % refusal
    if r + 1 > h
        E(end+1, :) = [pos s]; %#ok<AGROW>
    else
        stack{end+1} = [pos s / 2 r + 1]; %#ok<AGROW>
    end
end
E = unique(E, 'rows');
end

function w = worst_calls(pos, s, h)
% longest outcome sequence of the step rule from (pos, s) with h further
% refusals allowed (refused at refusal h + 1)
w = rec(pos, s, 0);
    function n = rec(p, s1, r)
        f = p + s1;
        n2 = 0;
        if r + 1 <= h, n2 = rec(p, s1 / 2, r + 1); end
        n1 = 0;
        if f ~= 1, n1 = rec(f, min(2 * s1, 1 - f), r); end
        n = 1 + max(n1, n2);
    end
end
