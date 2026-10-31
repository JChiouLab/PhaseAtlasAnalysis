function out = gnx_bvp_ode(x, y, par)
%GNX_BVP_ODE  12-state steady BVP of Guzzo_noXA on [0, ell] (packet §4a).
%
%   dydx = GNX_BVP_ODE(x, y, par)
%       State rows (gnx_model('descriptor').state_layout):
%       1 Rm  2 Rm'  3 Rc  4 Rc'  5 Bm  6 Bm'  7 Bc  8 Bc'  9 sigmaR  10 sigmaB
%       11 mR  12 mB, with
%       Rm'' = -(sR + sigmaR)/DM,  Rc'' = sR/DC,  Bm'' = -(sB + sigmaB)/DM,
%       Bc'' = sB/DC,  sigma' = 0,  mR' = Rm + Rc,  mB' = Bm + Bc.
%       Reactions are unclipped (gnx_reaction). This is bvp10_noxa.m:12-26
%       without the clamp, with the two shifts carried as constant states.
%   Singularity guard (packet §4c): if any evaluated Bm <= -K/2 the whole
%   right side is NaN and a guard counter is incremented.
%   n = GNX_BVP_ODE('guard_count') ; GNX_BVP_ODE('guard_reset')

persistent n_guard
if isempty(n_guard), n_guard = 0; end
if ischar(x)
    switch x
        case 'guard_count', out = n_guard;
        case 'guard_reset', n_guard = 0; out = 0;
        otherwise, error('GNX:ODE', 'Unknown action %s.', x);
    end
    return;
end
Bm = y(5, :);
if any(Bm <= -0.5 * par.K)
    n_guard = n_guard + 1;
    out = nan(12, size(y, 2));
    return;
end
[exR, exB] = gnx_reaction(y([1 3 5 7], :), par);
z = zeros(1, size(y, 2));
out = [y(2, :); -(exR + y(9, :)) / par.DM;
       y(4, :); -(-exR) / par.DC;
       y(6, :); -(exB + y(10, :)) / par.DM;
       y(8, :); -(-exB) / par.DC;
       z; z;
       y(1, :) + y(3, :); y(5, :) + y(7, :)];
end
