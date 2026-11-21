function solinit = gnx_guess(mode, varargin)
%GNX_GUESS  Starting guesses for the 12-state M. xanthus (Guzzo_noXA) half-cell BVP.
%
%   solinit = GNX_GUESS('profile', x, g4)
%       g4 is 4 x n (rows Rm Rc Bm Bc) on the increasing mesh x. The convention of the reference M. xanthus code:
%       derivative rows by gradient(g, x), shifts
%       zero, accumulators cumsum(g_R)*dx and cumsum(g_B)*dx with
%       dx = x(2) - x(1). On a non-uniform mesh (only the orientation
%       re-solve uses one) the accumulators are cumtrapz instead (cumsum*dx assumes a uniform mesh;
%       the boundary conditions fix the accumulators anyway).
%   solinit = GNX_GUESS('reverse', sol)
%       Orientation re-solve guess: on the reversed mesh
%       ell - fliplr(sol.x), concentration rows reversed, derivative rows
%       reversed and negated, shifts reversed (constants), accumulators
%       rebuilt as above from the reversed concentrations.
%   solinit = GNX_GUESS('constructed', x, nR, nB, amp_frac, w, par)
%       Constructed half-peak guess of the reference M. xanthus code, as a
%       'profile' guess.

switch mode
    case 'profile'
        solinit = from_profile(varargin{:});
    case 'reverse'
        sol = varargin{1};
        x = sol.x(:).';
        ell = x(end);
        xr = ell - fliplr(x);
        xr(1) = 0; xr(end) = ell;          % exact ends
        y = zeros(12, numel(x));
        y([1 3 5 7], :) = fliplr(sol.y([1 3 5 7], :));
        y([2 4 6 8], :) = -fliplr(sol.y([2 4 6 8], :));
        y([9 10], :) = fliplr(sol.y([9 10], :));
        y([11 12], :) = accumulators(xr, y([1 3 5 7], :));
        solinit = struct('x', xr, 'y', y);
    case 'constructed'
        [x, nR, nB, amp_frac, w, par] = varargin{:};
        g4 = constructed(x, nR, nB, amp_frac, w, par);
        solinit = from_profile(x, g4);
    otherwise
        error('GNX:GUESS', 'Unknown gnx_guess mode %s.', mode);
end
end

function solinit = from_profile(x, g4)
x = x(:).';
y = zeros(12, numel(x));
y(1, :) = g4(1, :); y(3, :) = g4(2, :);
y(5, :) = g4(3, :); y(7, :) = g4(4, :);
y(2, :) = gradient(g4(1, :), x);
y(4, :) = gradient(g4(2, :), x);
y(6, :) = gradient(g4(3, :), x);
y(8, :) = gradient(g4(4, :), x);
y([11 12], :) = accumulators(x, g4);
solinit = struct('x', x, 'y', y);
end

function A = accumulators(x, g4)
dx = x(2) - x(1);
gR = g4(1, :) + g4(2, :);
gB = g4(3, :) + g4(4, :);
if max(abs(diff(x) - dx)) <= 1e-12 * x(end)
    A = [cumsum(gR) * dx; cumsum(gB) * dx];        % accumulators
else
    A = [cumtrapz(x, gR); cumtrapz(x, gB)];
end
end

function g4 = constructed(xg, nR, nB, amp_frac, w, par)
% constructed half-peak guess (Lh = xg(end))
Lh = xg(end);
bump = exp(-(xg / w).^2);
bump = bump / trapz(xg, bump);
Bm = amp_frac * nB * bump;
Bc = (1 - amp_frac) * nB / Lh * ones(size(xg));
g = par.kR + par.kRB .* Bm;
shape = g ./ (g + par.dR);
Rm = shape .* (nR / Lh);
Rc = max((nR - trapz(xg, Rm)) / Lh, 0.01) * ones(size(xg));
g4 = [Rm; Rc; Bm; Bc];
end
