% bc_aug_v7  --  boundary conditions of the isolated half-cell BVP: no flux at both ends, and the
% two mass accumulators equal the prescribed GEF and PAK amounts nG, nP.
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
function res = bc_aug_v7(ya, yb, nG, nP)
    res = [ya(2); yb(2); ya(4); yb(4); ya(6); yb(6); ya(8); yb(8); ya(10); yb(10); ...
           ya(13); yb(13) - nG; ya(14); yb(14) - nP];
end
