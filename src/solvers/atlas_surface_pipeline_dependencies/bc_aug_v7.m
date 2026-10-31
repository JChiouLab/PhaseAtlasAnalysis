% bc_aug_v7  --  shared V7 per-point BVP core (extracted 2026-07-09, Packet 1B).
% Single source of truth for run_hybrid_sweeper_pdefb.m AND build_mpss_branch.m.
% Extracted verbatim from run_hybrid_sweeper_pdefb.m; do not fork. Equations, BCs,
% tolerances, and mass conventions are unchanged.
function res = bc_aug_v7(ya, yb, nG, nP)
    res = [ya(2); yb(2); ya(4); yb(4); ya(6); yb(6); ya(8); yb(8); ya(10); yb(10); ...
           ya(13); yb(13) - nG; ya(14); yb(14) - nP];
end
