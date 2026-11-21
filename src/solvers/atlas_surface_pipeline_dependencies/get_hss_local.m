% get_hss_local  --  uniform (no-peak) state of the fission yeast model for given total GEF and PAK amounts.
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
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
