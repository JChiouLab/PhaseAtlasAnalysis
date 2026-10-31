% hss_eqs  --  shared V7 per-point BVP core (extracted 2026-07-09, Packet 1B).
% Single source of truth for run_hybrid_sweeper_pdefb.m AND build_mpss_branch.m.
% Extracted verbatim from run_hybrid_sweeper_pdefb.m; do not fork. Equations, BCs,
% tolerances, and mass conventions are unchanged.
function F = hss_eqs(y, par, L, TG, TP)
    u = y(1); pv = y(2);
    q = (par.c/par.d)*u*pv;
    w = (par.d/par.e)*q;
    v = (par.b*u + par.d*q) / (par.a*u^2 + par.f);
    F = [L*(u + q + v + w) - TG;
         L*(pv + q) - TP];
end
