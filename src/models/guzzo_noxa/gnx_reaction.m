function [sR, sB, absR, absB] = gnx_reaction(u, par)
%GNX_REACTION  Guzzo_noXA reaction terms, UNCLIPPED (model page; Guzzo_noXA.m:219-233).
%
%   [sR, sB] = GNX_REACTION(u, par)
%   [sR, sB, absR, absB] = GNX_REACTION(u, par)
%
% u is 4 x N with rows Rm, Rc, Bm, Bc (field order); par carries kR kRB kB kBB
% dR dB dBA K KA (already scaled, gnx_model('parameters')). Returns 1 x N rows
%   sR = (kR + kRB*Bm)*Rc - dR*Rm
%   sB = (kB + kBB*Bm)*Bc - dB*Bm/(1 + Bm/K) - dBA*Aeff*Bm^2,
%   Aeff = Rm/(1 + (Bm/KA)^2).
% The source vector is [sR; -sR; sB; -sB]. Concentrations are never clamped;
% positivity is an acceptance check, not a change to the equations. Expression
% order follows bvp10_noxa.m:16-20 so positive states give the same values.
% absR, absB: sum of the absolute values of the terms (scale for residual
% tolerances, packet §4e, §5.1 Z1).

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
