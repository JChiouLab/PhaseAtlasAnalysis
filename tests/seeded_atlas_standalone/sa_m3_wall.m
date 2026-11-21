function c = sa_m3_wall(dd, ex)
%SA_M3_WALL  Stub rule: every call AT a node with k1 = 2 and |k2| <= 1 is
%   refused unless the guess is exact (only the secant rescue guess is); calls at
%   substep coordinates succeed. Forces refused edges and the rescue ladder.
c = 'S';
if abs(dd(1) - 1.3) < 1e-9 && abs(dd(2)) < 0.7 && ~ex, c = 'R'; end
end
