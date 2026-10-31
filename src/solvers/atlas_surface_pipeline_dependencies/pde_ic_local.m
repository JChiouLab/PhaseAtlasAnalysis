% pde_ic_local  --  shared V7 per-point BVP core (extracted 2026-07-09, Packet 1B).
% Single source of truth for run_hybrid_sweeper_pdefb.m AND build_mpss_branch.m.
% Extracted verbatim from run_hybrid_sweeper_pdefb.m; do not fork. Equations, BCs,
% tolerances, and mass conventions are unchanged.
function u0 = pde_ic_local(xq, y_arr, x_arr)
    u0 = zeros(5,1);
    for k = 1:5
        u0(k) = interp1(x_arr, y_arr(k,:), xq, 'pchip');
    end
end
