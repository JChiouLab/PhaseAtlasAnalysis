function [m, x_at, k_at] = gnx_hermite_min(x, y, yp)
%GNX_HERMITE_MIN  Exact minimum of bvp4c's piecewise cubic Hermite interpolant.
%
%   [m, x_at, k_at] = GNX_HERMITE_MIN(x, y, yp)
%   x (1 x n increasing), y and yp (1 x n) node values and slopes of one row.
%   On each interval [x_k, x_k+1] deval (ntrp3h) evaluates
%     p(s) = y_k + h*d*s^3 + h*c*s^2 + h*yp_k*s,  s in [0, 1],
%     slope = (y_k+1 - y_k)/h, c = 3*slope - 2*yp_k - yp_k+1,
%     d = yp_k + yp_k+1 - 2*slope.
%   Its minimum is the smaller of the end values and p at the real roots in
%   (0, 1) of p'(s)/h-scaled = 3d s^2 + 2c s + yp_k (closed form, stable
%   quadratic formula). Returns the global minimum m, its abscissa x_at and
%   the interval index k_at (packet §4b, R2).

x = x(:).'; y = y(:).'; yp = yp(:).';
h = diff(x);
y0 = y(1:end-1); y1 = y(2:end); m0 = yp(1:end-1); m1 = yp(2:end);
slope = (y1 - y0) ./ h;
c = 3 * slope - 2 * m0 - m1;
d = m0 + m1 - 2 * slope;
A = 3 * d; B = 2 * c; C = m0;
disc = B.^2 - 4 * A .* C;
sq = sqrt(max(disc, 0));
sg = sign(B); sg(sg == 0) = 1;
q = -0.5 * (B + sg .* sq);
s1 = q ./ A;
s2 = C ./ q;
s1(disc < 0) = NaN; s2(disc < 0) = NaN;
p = @(s) y0 + h .* d .* s.^3 + h .* c .* s.^2 + h .* m0 .* s;
v1 = p(s1); v2 = p(s2);
v1(~(s1 > 0 & s1 < 1) | ~isfinite(v1)) = Inf;
v2(~(s2 > 0 & s2 < 1) | ~isfinite(v2)) = Inf;
cand = [y0; y1; v1; v2];                  % 4 x (n-1)
xs = [x(1:end-1); x(2:end); x(1:end-1) + s1 .* h; x(1:end-1) + s2 .* h];
[mk, r] = min(cand, [], 1);
[m, k_at] = min(mk);
x_at = xs(r(k_at), k_at);
if any(~isfinite(y)) || any(~isfinite(yp))
    m = NaN; x_at = NaN; k_at = NaN;      % fail closed on non-finite input
end
end
