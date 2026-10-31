function J = gnx_jacobian(u, par)
%GNX_JACOBIAN  Analytic Jacobian of the Guzzo_noXA source vector.
%
%   J = GNX_JACOBIAN(u, par)
%
% u is a 4 x 1 state [Rm; Rc; Bm; Bc] (or 4 x N; then J is 4 x 4 x N).
% J(i,j) = d s_i / d u_j for s = [sR; -sR; sB; -sB] (gnx_reaction). Derived:
%   dsR/dRm = -dR,  dsR/dRc = kR + kRB*Bm,  dsR/dBm = kRB*Rc,  dsR/dBc = 0
%   dsB/dRm = -dBA*Bm^2/(1+(Bm/KA)^2),  dsB/dRc = 0,
%   dsB/dBm = kBB*Bc - dB/(1+Bm/K)^2 - dBA*Rm*2*Bm/(1+(Bm/KA)^2)^2,
%   dsB/dBc = kB + kBB*Bm.
% The structural zeros dsR/dBc and dsB/dRc are exact zeros (packet §5.1 Z3).

N = size(u, 2);
J = zeros(4, 4, N);
for n = 1:N
    Rm = u(1, n); Rc = u(2, n); Bm = u(3, n); Bc = u(4, n);
    q = 1 + (Bm / par.KA)^2;
    rR = [-par.dR, par.kR + par.kRB * Bm, par.kRB * Rc, 0];
    rB = [-par.dBA * Bm^2 / q, 0, ...
          par.kBB * Bc - par.dB / (1 + Bm / par.K)^2 - par.dBA * Rm * 2 * Bm / q^2, ...
          par.kB + par.kBB * Bm];
    J(:, :, n) = [rR; -rR; rB; -rB];
end
end
