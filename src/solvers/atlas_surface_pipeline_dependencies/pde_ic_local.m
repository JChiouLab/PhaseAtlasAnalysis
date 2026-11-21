% pde_ic_local  --  initial condition for pdepe: interpolates a stored profile onto the mesh.
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
function u0 = pde_ic_local(xq, y_arr, x_arr)
    u0 = zeros(5,1);
    for k = 1:5
        u0(k) = interp1(x_arr, y_arr(k,:), xq, 'pchip');
    end
end
