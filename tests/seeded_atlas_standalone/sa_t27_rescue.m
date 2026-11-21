function rows = sa_t27_rescue(E)
%SA_T27_RESCUE  Rescue ladder. Wall nodes (k1 = 2, |k2| <=
%   1) refuse every non-exact guess at their own coordinate, so continuation
%   edges and the R1 diagonal edge are refused; only the R2 secant guess
%   2*deval(P1, x1) - deval(P2, x1) (exact for linear targets) succeeds, from
%   collinear roots of one lineage. Caps: <= 9 calls for R1, 1 per R2
%   direction, rescue.max_calls_per_node per node (a cap of 5 cuts R1).
T = 'T27 rescue';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't27');
res = sa_m3_run(fullfile(d, 'a'), [3 1], 'fail', @(dd, b, ex) sa_m3_wall(dd, ex));
Q = res.Q; L = res.L;
wall = find(Q.box & abs(L.index(:, 1) - 2) == 0 & abs(L.index(:, 2)) <= 1).';
rr = Q.rescue_result(wall).';
rows = sa_check(rows, T, 'the 3 wall nodes are rescued by the secant rung', ...
    all(strcmp(rr, 'rescued_secant')), strjoin(rr, ','));
E1 = Q.edges(strcmp({Q.edges.kind}, 'rescue_diagonal'));
E2 = Q.edges(strcmp({Q.edges.kind}, 'rescue_secant'));
rows = sa_check(rows, T, 'R1 diagonal edges refused within 2h+1 = 9 calls each; R2 one call each', ...
    ~isempty(E1) && all(strcmp({E1.status}, 'refused')) && all([E1.calls] <= 9) && all([E2.calls] == 1), ...
    sprintf('R1 calls %s; R2 calls %s', mat2str([E1.calls]), mat2str([E2.calls])));
A = res.attempts;
sa = A(strcmp({A.route}, 'rescue_secant') & strcmp({A.outcome}, 'accepted'));
col = true;
for k = 1:numel(sa)
    t = regexp(sa(k).seed_left, 'secant k:(-?\d+),(-?\d+) k:(-?\d+),(-?\d+)', 'tokens', 'once');
    p1 = str2double(t(1:2)); p2 = str2double(t(3:4)); q = sa(k).index;
    col = col && isequal(q - p1, p1 - p2) && sum(abs(q - p1)) == 1;
end
rows = sa_check(rows, T, 'every accepted secant used P1 = q - e and P2 = q - 2e (collinear, one lattice step)', ...
    col && numel(sa) == 3, sprintf('%d secant roots', numel(sa)));
th = sa_m3_stub('params');
s1 = sa_m3_stub('sol', th, [2 2 2 2], 'A'); s2 = sa_m3_stub('sol', th, [1.5 2.5 2 2], 'A');
g = sa_rescue('secant_guess', s1, s2);
rows = sa_check(rows, T, 'secant guess = 2*deval(P1, x1) - deval(P2, x1) on P1''s mesh, all states', ...
    isequal(g.x, s1.x) && isequal(g.y, 2 * deval(s1, s1.x) - deval(s2, s1.x)) && size(g.y, 1) == 8);
rc = Q.rescue_calls(wall);
rows = sa_check(rows, T, 'rescue calls per node <= max_calls_per_node (16)', all(rc <= 16), mat2str(rc.'));
res2 = sa_m3_run(fullfile(d, 'cap5'), [3 1], 'fail', @(dd, b, ex) sa_m3_wall(dd, ex), ...
    'settings', struct('rescue', struct('max_calls_per_node', 5)));
Q2 = res2.Q;
rr2 = Q2.rescue_result(wall).';
E12 = Q2.edges(strcmp({Q2.edges.kind}, 'rescue_diagonal'));
rows = sa_check(rows, T, 'cap 5: R1 cut at 5 calls (rescue_cap), no later rung, result refused_cap', ...
    all(strcmp(rr2, 'refused_cap')) && all(strcmp({E12.status}, 'rescue_cap')) && all([E12.calls] == 5) && ...
    all(Q2.rescue_calls(wall) == 5) && ~any(strcmp({Q2.edges.kind}, 'rescue_secant')), ...
    sprintf('%s; R1 calls %s', strjoin(rr2, ','), mat2str([E12.calls])));
rmpath(E.stub);
end
