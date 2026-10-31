% pde_bc_v7  --  shared V7 per-point BVP core (extracted 2026-07-09, Packet 1B).
% Single source of truth for run_hybrid_sweeper_pdefb.m AND build_mpss_branch.m.
% Extracted verbatim from run_hybrid_sweeper_pdefb.m; do not fork. Equations, BCs,
% tolerances, and mass conventions are unchanged.
function [pl,ql,pr,qr] = pde_bc_v7(~,~,~,~,~)
    pl = zeros(5,1); ql = ones(5,1);
    pr = zeros(5,1); qr = ones(5,1);
end
