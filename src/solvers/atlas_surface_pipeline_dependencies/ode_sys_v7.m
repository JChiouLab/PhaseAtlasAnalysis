% ode_sys_v7  --  shared V7 per-point BVP core (extracted 2026-07-09, Packet 1B).
% Single source of truth for run_hybrid_sweeper_pdefb.m AND build_mpss_branch.m.
% Extracted verbatim from run_hybrid_sweeper_pdefb.m; do not fork. Equations, BCs,
% tolerances, and mass conventions are unchanged.
function dydx = ode_sys_v7(~, y, par, D)
    u = max(0,y(1)); q = max(0,y(3)); v = max(0,y(5)); w = max(0,y(7)); pv = max(0,y(9));
    s1 = y(11); s2 = y(12);
    Ru = (par.a*u^2 + par.f)*v - par.b*u - par.c*u*pv;
    Rq = par.c*u*pv - par.d*q;
    Rv = -(par.a*u^2 + par.f)*v + par.b*u + par.e*w;
    Rw = par.d*q - par.e*w;
    Rp = -par.c*u*pv + par.d*q;
    dydx = [y(2);  -(Ru + s1)/D(1);
            y(4);  -Rq/D(2);
            y(6);  -Rv/D(3);
            y(8);  -Rw/D(4);
            y(10); -(Rp + s2)/D(5);
            0; 0;
            u + q + v + w;
            pv + q];
end
