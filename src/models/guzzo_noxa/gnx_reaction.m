function [sR, sB, absR, absB] = gnx_reaction(u, par)
%GNX_REACTION  M. xanthus reaction terms (Eqs. 3-4), unclipped.
%
%   [sR, sB] = GNX_REACTION(u, par)
%   [sR, sB, absR, absB] = GNX_REACTION(u, par)
%
% u is 4 x N with rows Rm, Rc, Bm, Bc (paper rm, rc, bm, bc); par carries kR kRB kB kBB
% dR dB dBA K KA (already scaled, gnx_model('parameters')). Returns 1 x N rows
%   sR = (kR + kRB*Bm)*Rc - dR*Rm
%   sB = (kB + kBB*Bm)*Bc - dB*Bm/(1 + Bm/K) - dBA*Aeff*Bm^2,
%   Aeff = Rm/(1 + (Bm/KA)^2).
% The source vector is [sR; -sR; sB; -sB]. Concentrations are never clamped;
% positivity is an acceptance check, not a change to the equations.
% absR, absB: sum of the absolute values of the terms (scale for residual
% tolerances).

Rm = u(1, :); Rc = u(2, :); Bm = u(3, :); Bc = u(4, :);
aeff = Rm ./ (1 + (Bm / par.KA).^2);
sR = (par.kR + par.kRB .* Bm) .* Rc - par.dR .* Rm;
sB = (par.kB + par.kBB .* Bm) .* Bc ...
    - par.dB .* Bm ./ (1 + Bm / par.K) ...
    - par.dBA .* aeff .* Bm.^2;
if nargout > 2
    absR = abs(par.kR .* Rc) + abs(par.kRB .* Bm .* Rc) + abs(par.dR .* Rm);
    absB = abs(par.kB .* Bc) + abs(par.kBB .* Bm .* Bc) + ...
        abs(par.dB .* Bm ./ (1 + Bm / par.K)) + abs(par.dBA .* aeff .* Bm.^2);
end
end
