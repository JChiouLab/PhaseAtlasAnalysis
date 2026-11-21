% pde_bc_v7  --  no-flux boundary conditions of the pdepe half-cell PDE (fission yeast model).
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
function [pl,ql,pr,qr] = pde_bc_v7(~,~,~,~,~)
    pl = zeros(5,1); ql = ones(5,1);
    pr = zeros(5,1); qr = ones(5,1);
end
