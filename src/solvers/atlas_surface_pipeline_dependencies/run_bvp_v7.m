% run_bvp_v7  --  shared V7 per-point BVP core (extracted 2026-07-09, Packet 1B).
% Single source of truth for run_hybrid_sweeper_pdefb.m AND build_mpss_branch.m.
% Extracted verbatim from run_hybrid_sweeper_pdefb.m; do not fork. Equations, BCs,
% tolerances, and mass conventions are unchanged.
function [bvp_sol, success] = run_bvp_v7(guess_sol, L, par, D, nG, nP)
    options = bvpset('RelTol', 1e-6, 'NMax', 10000);
    bvp_sol = []; success = false;
    ws = warning('off', 'all');
    cleanup = onCleanup(@() warning(ws));
    try
        bvp_sol = bvp4c(@(x,y) ode_sys_v7(x,y,par,D), ...
                        @(ya,yb) bc_aug_v7(ya,yb,nG,nP), ...
                        guess_sol, options);
        if length(bvp_sol.x) >= 3000, bvp_sol = []; return; end
        if any(isnan(bvp_sol.y(:))) || any(isinf(bvp_sol.y(:))), bvp_sol = []; return; end
        if min(bvp_sol.y([1,3,5,7,9],:), [], 'all') < -1e-3, bvp_sol = []; return; end
        success = true;
    catch
        bvp_sol = []; success = false;
    end
end
