function rows = sa_t64_cutset(E)
%SA_T64_CUTSET  The exact cut recount sa_report('cutset'). Zero solves; planted attempt
%   files in a scratch folder.
%   At c = 10,000 with the call's own NMax 20,000: P > c -> cut; m >= c ->
%   cut; P = c with m < c -> not_cut; P = c = m -> cut; N0 > c with every
%   proposal <= c and m < c -> not_cut; a no-return with P > c -> cut, with
%   P <= c -> not_cut; P NaN -> unknown_not_recorded. Above the call's own
%   NMax (3,000): a mesh-bound (over_nmax) call with P > c -> cut; over_nmax
%   with P <= c -> unknown_above_nmax; a non-mesh-bound call -> not_cut. An
%   interrupted row is excluded; refinement and continuation are excluded
%   unless asked; a folder without traces prints the predicted classes labelled
%   as a lower bound and is never exact; the cap-loss
%   curve at 3,601 / 5,000 / 10,000; the csv and summary files. When the recall
%   fixture fixtures/m6/recall_A000125.mat exists (it is not shipped; the check
%   is skipped without it): A000125, a cut the prediction missed, is cut at 10,000 by its recorded P.
T = 'T64 exact recount';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
d = sa_m3_testdir(E, 't64');
% id, kind, outcome, nmax0, N0, proposals, m (NaN = no return), over, expected class at c = 10,000
P = {
    'A000001', 'cold', 'invalid',     20000, 300,   [5000 12000],  3000,  false, 'cut'                   % P > c
    'A000002', 'cold', 'invalid',     20000, 12000, [9000],        12000, false, 'cut'                   % m >= c
    'A000003', 'cold', 'accepted',    20000, 300,   [10000 6000],  6000,  false, 'not_cut'               % P = c, m < c
    'A000004', 'cold', 'invalid',     20000, 300,   [10000],       10000, false, 'cut'                   % P = c = m
    'A000005', 'cold', 'accepted',    20000, 12000, [8000 6000],   6000,  false, 'not_cut'               % N0 > c
    'A000006', 'cold', 'no_return',   20000, 300,   [11000],       NaN,   false, 'cut'                   % no-return, P > c
    'A000007', 'cold', 'no_return',   20000, 300,   [9000],        NaN,   false, 'not_cut'               % no-return, P <= c
    'A000008', 'cold', 'invalid',     20000, NaN,   [],            4000,  false, 'unknown_not_recorded'  % no trace
    'A000009', 'cold', 'invalid',     3000,  300,   [2000 12000],  2000,  true,  'cut'                   % over_nmax, P > c > NMax
    'A000010', 'cold', 'invalid',     3000,  300,   [2000 4000],   2000,  true,  'unknown_above_nmax'    % over_nmax, P <= c
    'A000011', 'cold', 'accepted',    3000,  300,   [2000 2500],   2500,  false, 'not_cut'               % not mesh-bound, c > NMax
    'A000012', 'interrupted', 'interrupted', 20000, NaN, [],       NaN,   false, ''                      % excluded
    'A000013', 'refinement', 'refined_valid', 20000, 300, [15000 9000], 9000, false, 'cut'               % only when asked
    'A000014', 'continuation', 'accepted', 20000, 300, [12000 8000], 8000, false, 'cut'                  % only when asked
    'A000015', 'cold', 'accepted',    20000, 300,   [4000 3000],   3000,  false, 'not_cut'
    };
F = fullfile(d, 'planted');
plant(F, P);
X = sa_report('cutset', F, 10000, {}, struct('out_dir', fullfile(d, 'out')));
cold = ~ismember(P(:, 2), {'interrupted', 'refinement', 'continuation'});
want = P(cold, 9).';
got = X.class;
ids = {X.calls.attempt_id};
rows = sa_check(rows, T, ['exact rule at c = 10,000 (NMax 20,000): P > c, m >= c, P = c = m, a no-return with P > c -> ' ...
    'cut; P = c with m < c, N0 > c with proposals <= c and m < c, a no-return with P <= c -> not_cut; P NaN -> ' ...
    'unknown_not_recorded; above the own NMax 3,000: over_nmax with P > c -> cut, over_nmax with P <= c -> ' ...
    'unknown_above_nmax, not mesh-bound -> not_cut'], isequal(ids, P(cold, 1).') && isequal(got, want), ...
    sprintf('%s', strjoin(strcat(ids, ':', got), ' ')));
rows = sa_check(rows, T, ['an interrupted row is excluded (counted apart); refinement and continuation excluded by ' ...
    'default; the count is NOT exact with unknown classes present'], X.n_interrupted_excluded == 1 && ...
    X.n == sum(cold) && ~X.exact && X.n_cut == 5 && X.n_unknown_not_recorded == 1 && X.n_unknown_above_nmax == 1, ...
    X.lines{1});
Xa = sa_report('cutset', F, 10000, {'cold', 'refinement', 'continuation'}, struct('write', false));
ka = strcmp({Xa.calls.attempt_id}, 'A000013') | strcmp({Xa.calls.attempt_id}, 'A000014');
rows = sa_check(rows, T, 'refinement and continuation calls are classed when asked (both cut at 10,000)', ...
    Xa.n == sum(cold) + 2 && all(strcmp(Xa.class(ka), 'cut')), Xa.lines{1});
% exact folder (only exact classes) and the cap-loss curve
G = fullfile(d, 'exact');
Pe = P([1:7 11 15], :);
plant(G, Pe);
Xe = sa_report('cutset', G, 10000, {}, struct('out_dir', fullfile(d, 'out_exact')));
cv = Xe.curve;
% at 3,601: 1 (P 12000), 2 (m 12000), 3 (P 10000), 4, 5 (N0 12000 is not tested; P 8000 > 3601), 6, 7 (P 9000), 15 (P 4000)
% cut; 11 (NMax 3,000 < 3,601: not mesh-bound, P 2500 <= c, m 2500 < c) not_cut
okc = Xe.exact && isequal([cv.c], [3601 5000 10000]) && isequal([cv.n_cut], [8 7 4]) && ...
    isequal([cv.n_cut_accepted], [3 2 0]) && all([cv.exact]);
sm = fileread(fullfile(d, 'out_exact', 'cutset_10000_summary.txt'));
cs = strsplit(strtrim(fileread(fullfile(d, 'out_exact', 'cutset_10000.csv'))), newline);
rows = sa_check(rows, T, ['cap-loss curve at 3,601 / 5,000 / 10,000 (cut calls; of them stored accepted returns; exact): ' ...
    '8/3, 7/2, 4/0, all exact; cutset_10000.csv one row per call; the summary file holds the one-line summary, the ' ...
    'histogram of P and the curve'], okc && numel(cs) == size(Pe, 1) + 1 && contains(sm, 'EXACT') && ...
    contains(sm, 'histogram of P') && contains(sm, 'cap-loss curve') && contains(cs{1}, 'max_proposed_mesh,over_nmax,class'), ...
    sprintf('n_cut %s accepted %s exact %d', mat2str([cv.n_cut]), mat2str([cv.n_cut_accepted]), Xe.exact));
% a folder without traces: the prediction fallback, labelled a lower bound
H = fullfile(d, 'legacy');
Pl = P([1 2 6 7], :);
plant(H, Pl, true);
Xl = sa_report('cutset', H, 10000, {}, struct('out_dir', fullfile(d, 'out_legacy'), 'cut_time_s', 5));
rows = sa_check(rows, T, ['a folder without traces: every call unknown_not_recorded, NOT exact, and the call classes ' ...
    'printed as "RC-5.1 prediction (lower bound, K-37)" (certain: NMax-bound or returned mesh >= c; possible: ' ...
    'no-returns above the cut time; A000001, which proposed 12,000 but returned 3,000, is missed: K-37)'], ~Xl.exact && Xl.n_unknown_not_recorded == 4 && ...
    any(contains(Xl.lines, 'Prediction (lower bound) for the')) && isequal(Xl.rc51.certain, {'A000002'}) && ...
    isequal(Xl.rc51.possible, {'A000006', 'A000007'}), strjoin(Xl.lines(end), ''));
% the recall fixture (A000125)
fx = fullfile(E.tests, 'fixtures', 'm6', 'recall_A000125.mat');
if exist(fx, 'file') == 2
    Z = load(fx);
    Q = {'A000125', 'cold', Z.outcome, Z.nmax, Z.mesh_trace.initial, Z.mesh_trace.proposals, Z.returned_mesh, ...
        Z.mesh_trace.over_nmax, ''};
    K = fullfile(d, 'a125');
    plant(K, Q, false, Z.mesh_trace);
    Xk = sa_report('cutset', K, 10000, {}, struct('write', false));
    rows = sa_check(rows, T, ['fixture: the B1 recall of A000125 (uncapped, NMax 20,000; V2a''s unpredicted cut): cut at ' ...
        '10,000 by its recorded P (> 10,000) although it returned below the cap'], ...
        strcmp(Xk.class{1}, 'cut') && Z.mesh_trace.max_proposed > 10000 && Z.returned_mesh < 10000, ...
        sprintf('P %g, returned %g, proposals %s', Z.mesh_trace.max_proposed, Z.returned_mesh, mat2str(Z.mesh_trace.proposals)));
end
end

% =========================================================================

function plant(F, P, notrace, trace)
if nargin < 3, notrace = false; end
mkdir(fullfile(sa_layout('state', F), 'attempts'));
for k = 1:size(P, 1)
    id = P{k, 1};
    a = sa_report('blank');
    a.attempt_id = id; a.kind = P{k, 2}; a.outcome = P{k, 3};
    m = P{k, 7};
    a.returned_mesh = m; a.wall_s = 10;
    a.convergence_class = 'CONVERGED';
    if isfinite(m) && ~isempty(P{k, 6}) && P{k, 6}(end) > P{k, 4}, a.convergence_class = 'UNRESOLVED_RELTOL_NMAX'; end
    if strcmp(a.outcome, 'no_return'), a.convergence_class = 'NO_RETURN'; end
    Z = struct('attempt', a, 'opts', struct('nmax', P{k, 4}, 'mesh_reject', P{k, 4}));
    if ~notrace && ~strcmp(a.outcome, 'interrupted')
        if nargin >= 4
            Z.mesh_trace = trace;
        else
            n0 = P{k, 5}; pr = P{k, 6};
            tr = sa_mesh_trace('blank');
            tr.initial = n0; tr.proposals = pr; tr.n_proposals = numel(pr);
            if ~isempty(pr), tr.max_proposed = max(pr); elseif isfinite(n0), tr.max_proposed = 0; end
            tr.over_nmax = P{k, 8};
            Z.mesh_trace = tr;
        end
        a.max_proposed_mesh = Z.mesh_trace.max_proposed;
        Z.attempt = a;
    end
    sa_write('mat', fullfile(sa_layout('state', F), 'attempts', [id '.mat']), Z);
end
end
