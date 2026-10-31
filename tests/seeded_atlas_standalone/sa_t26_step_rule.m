function rows = sa_t26_step_rule(E)
%SA_T26_STEP_RULE  Milestone-3 T26 (RC-1 with the review addendum): the
%   adaptive rule refuses an edge when refusals exceed h, within 2h+1 calls;
%   adversarial alternating refusal/success -> REFUSAL at the cap (9 calls at
%   h = 4, 13 at h = 6); R^4 S^5 completes in 9; myxo_restart <= 31 calls
%   (attained by failing every level at its last substep). The stub's
%   scripted sequence starts with the anchor's cold call; the first edge is
%   anchor -> k:-1,0 (box [1 0], rescue off).
T = 'T26 step rule';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't26');
off = struct('routes', struct('rescue', false));
e = first_edge(fullfile(d, 'alt4'), ['S' repmat('RS', 1, 8)], off);
rows = sa_check(rows, T, 'h = 4, alternating R,S,R,S,...: REFUSED at r = h+1 = 5 after 2h+1 = 9 calls', ...
    strcmp(e.status, 'refused') && e.calls == 9 && e.refusals == 5 && e.fraction_valid == 15/16, desc(e));
e = first_edge(fullfile(d, 'r4s5'), ['S' 'RRRRSSSSS'], off);
rows = sa_check(rows, T, 'h = 4, R^4 S^5: completes in 9 calls', strcmp(e.status, 'completed') && e.calls == 9, desc(e));
o6 = off; o6.continuation = struct('max_halvings', 6);
e = first_edge(fullfile(d, 'alt6'), ['S' repmat('RS', 1, 10)], o6);
rows = sa_check(rows, T, 'h = 6, alternating: REFUSED at r = 7 after 2h+1 = 13 calls', ...
    strcmp(e.status, 'refused') && e.calls == 13 && e.refusals == 7, desc(e));
e = first_edge(fullfile(d, 'allR'), ['S' repmat('R', 1, 9)], off);
rows = sa_check(rows, T, 'h = 4, all refusals: refused after h+1 = 5 calls, last valid fraction 0', ...
    strcmp(e.status, 'refused') && e.calls == 5 && e.fraction_valid == 0, desc(e));
rs = RandStream('mt19937ar', 'Seed', 26);
ok = true; mx = 0; n = 0;
for t = 1:12
    seq = ['S' char('R' + ('S' - 'R') * (rand(rs, 1, 12) < 0.5))];
    e = first_edge(fullfile(d, sprintf('rnd%02d', t)), seq, off);
    mx = max(mx, e.calls); n = n + 1;
    ok = ok && e.calls <= 9 && (strcmp(e.status, 'completed') || e.refusals == 5);
end
rows = sa_check(rows, T, '12 random R/S sequences: every edge <= 9 calls; refused exactly at r = 5', ok, ...
    sprintf('max calls %d over %d edges', mx, n));
om = off; om.continuation = struct('step_rule', 'myxo_restart');
seq = ['S' 'R' 'SR' 'SSSR' 'SSSSSSSR' repmat('S', 1, 16)];
e = first_edge(fullfile(d, 'myxo31'), seq, om);
rows = sa_check(rows, T, 'myxo_restart: levels fail at their last substep -> completes at level 4 in 1+2+4+8+16 = 31 calls', ...
    strcmp(e.status, 'completed') && e.calls == 31, desc(e));
e = first_edge(fullfile(d, 'myxoR'), ['S' repmat('R', 1, 6)], om);
rows = sa_check(rows, T, 'myxo_restart: all refusals -> refused after h+1 = 5 calls (<= 31)', ...
    strcmp(e.status, 'refused') && e.calls == 5, desc(e));
rmpath(E.stub);
end

function e = first_edge(run_dir, seq, o)
res = sa_m3_run(run_dir, [1 0], 'seq', seq, 'settings', o);
e = res.Q.edges(1);
assert(e.p == res.Q.anchor_node);
end

function s = desc(e)
s = sprintf('status %s, calls %d, refusals %d, last valid fraction %g', e.status, e.calls, e.refusals, e.fraction_valid);
end
