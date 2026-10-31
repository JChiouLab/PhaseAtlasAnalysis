function ctx = sa_refinement(ctx, ~)
%SA_REFINEMENT  Refinement check, declared before any solve.
%
%   Sample: in each selected family, the accepted cold return with the largest
%   stats.maxres (ties: earliest attempt); fewer if a family has none.
%   Re-solve from that root's own returned profile with RelTol =
%   validity.refine_reltol, AbsTol = validity.refine_abstol and
%   NMax = mesh_reject = solver.mesh_reject. The refined return must pass the
%   whole validity boundary; the report gives the identity distance
%   (scaled_rms_v1) and the relative-sup distances, and flags distances above
%   validity.refine_expected_distance. A refined return that fails the
%   boundary is reported as "refinement not achieved", never dropped; a sample
%   the budget does not reach is "not attempted (<budget reason>)".
%   Refinement calls count against budget.max_calls; refined returns never
%   enter the root store. Every call goes through the executor
%   sa_search('call') (journal, replay on resume).
%   The parent's sol is loaded from its attempt file by attempt id
%   (sa_search 'get_sol'), live and on resume alike (a MAT round trip is
%   lossless); no cold sol is kept in memory. A second argument (the former
%   in-memory sols cache) is accepted and ignored.

S = ctx.S; M = ctx.M;
cold = ctx.attempts(strcmp({ctx.attempts.kind}, 'cold') & strcmp({ctx.attempts.outcome}, 'accepted'));
rows = ctx.refinement;
for f = 1:numel(ctx.T.families)
    fam = ctx.T.families{f};
    pool = cold(strcmp({cold.family}, fam));
    if isempty(pool), continue; end
    [~, k] = max([pool.maxres]);           % max returns the first of ties
    A = pool(k);
    r = sa_report('blank_refinement_row');
    r.family = fam; r.parent_attempt = A.attempt_id; r.parent_root = A.root_id;
    r.target_key = A.target_key; r.parent_maxres = A.maxres;
    tg = ctx.T.list(A.target_number);
    u = sa_search('unit');
    u.route = 'refinement'; u.kind = 'refinement'; u.phase = 'refinement';
    u.target_number = A.target_number; u.family = tg.family; u.key = tg.key;
    u.index = tg.index; u.fraction = 1; u.coordinate = tg.coordinate;
    u.parent_attempt = A.attempt_id; u.parent_root = A.root_id;
    u.guess_source = ['attempt ' A.attempt_id];
    u.seed_left = A.seed_left; u.seed_right = A.seed_right;
    if ctx.stop
        ok = false; why = ctx.stop_reason;
    else
        [ok, why] = budget_check(ctx, u);
    end
    if ~ok
        r.status = sprintf('not attempted (%s)', why);
        rows(end+1) = r; %#ok<AGROW>
        continue
    end
    sol0 = sa_search('get_sol', ctx, A.attempt_id);
    guess = struct('x', sol0.x, 'y', sol0.y);
    [ctx, A2, R] = sa_search('call', ctx, u, tg, guess, sa_search('solver_for', S, 'refinement'));
    if isempty(A2)                          % the executor refused (stop)
        r.status = sprintf('not attempted (%s)', ctx.stop_reason);
        rows(end+1) = r; %#ok<AGROW>
        ctx.stop = false;                   % a refinement refusal is not a stop
        continue
    end
    r.refine_attempt = A2.attempt_id;
    r.refined_maxres = A2.maxres;
    if A2.valid
        P = sa_identity('pair', M, sol0, R.sol, ctx.K, S);
        r.status = 'refined';
        r.scaled_rms = P.scaled_rms;
        r.relsup_identity = P.relsup_identity;
        r.relsup_all = P.relsup_all;
        r.ratio = P.ratio;
        r.within_expectation = P.scaled_rms <= S.validity.refine_expected_distance;
    else
        r.status = 'refinement not achieved';
        r.failing = A2.failing;
    end
    rows(end+1) = r; %#ok<AGROW>
    ctx.refinement = rows;
end
ctx.refinement = rows;
end

function [ok, why] = budget_check(ctx, u) %#ok<INUSD>
ok = true; why = '';
S = ctx.S;
if ctx.calls >= S.budget.max_calls
    ok = false; why = 'budget_calls';
elseif ~ctx.replay.active && toc(ctx.t_run) >= min(S.budget.max_wall_s, S.budget.wall_s)
    ok = false; why = 'budget_wall';
end
end
