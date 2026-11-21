% hss_eqs  --  algebraic equations for the uniform state at given total GEF (TG) and PAK (TP) amounts.
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
function F = hss_eqs(y, par, L, TG, TP)
    u = y(1); pv = y(2);
    q = (par.c/par.d)*u*pv;
    w = (par.d/par.e)*q;
    v = (par.b*u + par.d*q) / (par.a*u^2 + par.f);
    F = [L*(u + q + v + w) - TG;
         L*(pv + q) - TP];
end
