% get_hss_local  --  shared V7 per-point BVP core (extracted 2026-07-09, Packet 1B).
% Single source of truth for run_hybrid_sweeper_pdefb.m AND build_mpss_branch.m.
% Extracted verbatim from run_hybrid_sweeper_pdefb.m; do not fork. Equations, BCs,
% tolerances, and mass conventions are unchanged.
function [u, q, v, w, pv] = get_hss_local(TG, TP, L, par)
    opts = optimset('Display','off');
    [vals, ~, exitflag] = fsolve(@(y) hss_eqs(y, par, L, TG, TP), [0.5, 0.5], opts);
    if exitflag <= 0
        [vals, ~] = fsolve(@(y) hss_eqs(y, par, L, TG, TP), [0.1, 0.1], opts);
    end
    u  = vals(1);
    pv = vals(2);
    q  = (par.c/par.d)*u*pv;
    w  = (par.d/par.e)*q;
    v  = (par.b*u + par.d*q) / (par.a*u^2 + par.f);
end
