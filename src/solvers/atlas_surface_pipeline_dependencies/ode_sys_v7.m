% ode_sys_v7  --  right-hand side of the isolated half-cell BVP (fission yeast model).
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
% State: u,u',q,q',v,v',w,w',p,p' (y1..y10), the two uniform sources s1, s2 (y11, y12) and the two
% mass accumulators (y13, y14). a = kgg, b = dg, c = kpg, d = dgp, e = k*, f = kg of Table 2.
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
