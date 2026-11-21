% pde_eq_v7  --  pdepe equations of the fission yeast model (u, q, v, w, p) for the half-cell PDE.
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
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
