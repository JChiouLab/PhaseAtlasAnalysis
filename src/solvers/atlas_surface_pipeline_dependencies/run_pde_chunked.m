% run_pde_chunked  --  chunked pdepe relaxation to steady state (fission yeast model; pdepe, not central differences + ode15s).
% Shared by run_hybrid_sweeper_pdefb.m and build_mpss_branch.m; keep a single copy.
function sol_out = run_pde_chunked(y0, x, L, par, D, T_end, T_chunk, max_chunks)
    % PDE relaxation: chunked integration with mass-conserving
    % cos-perturbation IC and early stopping when max|du/dt| chunk-to-chunk is
    % below 1e-6.  T_chunk and max_chunks are chosen so T_end <= T_chunk*max_chunks.
    sol_out = [];
    if T_end <= 0, return; end
    if T_chunk > T_end, T_chunk = T_end; end
    if T_chunk <= 0, T_chunk = T_end / 2; end
    if max_chunks < 1, max_chunks = 1; end

    % Broadcast scalar HSS initial condition to full spatial grid
    N = length(x);
    if size(y0, 2) == 1
        y0 = repmat(y0, 1, N);
    end

    % Apply mass-conserving cos perturbation (zero mean over [0,L])
    u_ss = y0(1,:); p_ss = y0(5,:);
    perturb = 0.5 * u_ss .* cos(pi * x / L);
    y0(1,:) = u_ss + perturb;
    y0(5,:) = p_ss - perturb;

    curr_y  = y0;
    converged = false;
    ws = warning('off', 'all');
    cleanup = onCleanup(@() warning(ws));
    try
        for k = 1:max_chunks
            tspan = [0, T_chunk/2, T_chunk];
            sol = pdepe(0, ...
                @(xq,t,uu,dudx) pde_eq_v7(xq,t,uu,dudx,par,D), ...
                @(xq) pde_ic_local(xq, curr_y, x), ...
                @pde_bc_v7, x, tspan);
            frame_end = squeeze(sol(end,:,:))';
            if max(abs(frame_end(1,:) - curr_y(1,:))) < 1e-6
                curr_y = frame_end;
                converged = true;
                break;
            end
            curr_y = frame_end;
        end
    catch
        return;
    end
    if ~converged
        % Still return whatever we have - downstream BVP will determine usefulness.
    end

    sol_out.x = x;
    yf = zeros(14, length(x));
    yf(1,:) = curr_y(1,:); yf(3,:) = curr_y(2,:); yf(5,:) = curr_y(3,:);
    yf(7,:) = curr_y(4,:); yf(9,:) = curr_y(5,:);
    sol_out.y = yf;
end
