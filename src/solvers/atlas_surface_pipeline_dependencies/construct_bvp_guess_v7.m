% construct_bvp_guess_v7  --  builds a bvp4c starting guess (with sources and accumulators) from a profile.
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
function guess_sol = construct_bvp_guess_v7(sol_in)
    % No integral rescale: s1, s2 absorb the residual mass.
    x = sol_in.x; N = length(x);
    y = zeros(14, N);
    u = sol_in.y(1,:); q = sol_in.y(3,:);
    v = sol_in.y(5,:); w = sol_in.y(7,:); pv = sol_in.y(9,:);
    y(1,:) = u; y(3,:) = q; y(5,:) = v; y(7,:) = w; y(9,:) = pv;
    % Derivative w.r.t. x on the actual (possibly non-uniform) mesh.
    y(2,:)  = gradient(u , x);
    y(4,:)  = gradient(q , x);
    y(6,:)  = gradient(v , x);
    y(8,:)  = gradient(w , x);
    y(10,:) = gradient(pv, x);
    y(11,:) = 0; y(12,:) = 0;
    y(13,:) = cumtrapz(x, u + q + v + w);
    y(14,:) = cumtrapz(x, pv + q);
    guess_sol.x = x;
    guess_sol.y = y;
end
