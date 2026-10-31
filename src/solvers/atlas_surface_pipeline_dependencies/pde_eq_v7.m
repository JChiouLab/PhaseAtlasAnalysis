% pde_eq_v7  --  shared V7 per-point BVP core (extracted 2026-07-09, Packet 1B).
% Single source of truth for run_hybrid_sweeper_pdefb.m AND build_mpss_branch.m.
% Extracted verbatim from run_hybrid_sweeper_pdefb.m; do not fork. Equations, BCs,
% tolerances, and mass conventions are unchanged.
function [c,f,s] = pde_eq_v7(~, ~, u, dudx, par, D)
    c = ones(5,1);
    f = D .* dudx;
    Ru = (par.a*u(1)^2 + par.f)*u(3) - par.b*u(1) - par.c*u(1)*u(5);
    Rq = par.c*u(1)*u(5) - par.d*u(2);
    Rv = -(par.a*u(1)^2 + par.f)*u(3) + par.b*u(1) + par.e*u(4);
    Rw = par.d*u(2) - par.e*u(4);
    Rp = -par.c*u(1)*u(5) + par.d*u(2);
    s = [Ru; Rq; Rv; Rw; Rp];
end
