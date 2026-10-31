function varargout = sag_test_fixture(action, varargin)
%SAG_TEST_FIXTURE  Synthetic 24-row candidates and run settings for the 2C tests.
%
%   p = SAG_TEST_FIXTURE('params', f)        owner-setting parameters at signal f
%   s = SAG_TEST_FIXTURE('settings', run_dir) commissioning settings (packet 2C)
%   [sol, tgt] = SAG_TEST_FIXTURE('poly', theta, slope, n)
%       A Hermite-exact synthetic candidate: every concentration is a cubic
%       c(x) = a + b q(x) + e r(x) on [0, Lh] (q = 3(x/Lh)^2 - 2(x/Lh)^3,
%       r = x^3/(3 Lh^2)), so c'(0) = 0 (pole no flux) and c'(Lh) = e. The
%       right half has e_R = -e_L (seam flux match) and a b_R making the seam
%       continuous; shifts sigma = -J/Lh per half (so the section 2c
%       identities hold); accumulators are the exact quartic integrals and the
%       targets their end values. slope = 0 gives zero seam slopes and zero
%       shifts. yp holds the exact derivatives (conc: c', deriv: c'', acc:
%       density), so bvp4c's interpolant is the exact cubic. It solves no
%       reaction: only the gate structure is exercised.
%   sol = SAG_TEST_FIXTURE('uniform', theta, u, n)
%       The exact uniform 24-row solution at the 4-vector u (both halves).
%   [Y, YP] = SAG_TEST_FIXTURE('eval', theta, slope, x)
%       The polynomial fixture's exact rows and derivatives at any x.
%   R = SAG_TEST_FIXTURE('R', sol)            wrap as a corrector return

switch action
    case 'params'
        varargout{1} = struct('L', 6, 'r_speed', 0.35, 'dba_factor', varargin{1}, ...
            'rho_R', 1, 'rho_B', 1, 'res', 25, 'dn_half_range', 1.8);
    case 'settings'
        s = struct('model', 'guzzo_noxa', 'run_dir', varargin{1});
        s.lattice = struct('step', 0.15, 'anchor', [0 0], 'margin', 1e-4);
        s.solver = struct('reltol', 1e-8, 'abstol', 1e-10, 'nmax', 20000, 'mesh_reject', 20000);
        s.validity = struct('refine_reltol', 1e-10, 'refine_abstol', 1e-12, ...
            'refine_expected_distance', 1e-4);
        varargout{1} = s;
    case 'poly'
        [varargout{1}, varargout{2}] = poly(varargin{:});
    case 'uniform'
        varargout{1} = uniform(varargin{:});
    case 'eval'
        [varargout{1}, varargout{2}] = eval_poly(varargin{:});
    case 'R'
        varargout{1} = struct('sol', varargin{1}, 'opts', struct('reltol', 1e-8, ...
            'abstol', 1e-10, 'nmax', 20000, 'mesh_reject', 20000), ...
            'warning_id', '', 'warning_msg', '');
    otherwise
        error('SAG:FIXTURE', 'Unknown action %s.', action);
end
end

function P = coeffs(theta, slope)
% columns: Rm Rc Bm Bc ; rows a b e, left then right
aL = [0.8 1e-3 0.5 0.4];  bL = [0.3 5e-4 0.2 0.05];  eL = slope * [0.02 1e-5 0.05 0.002];
aR = [0.9 1.2e-3 0.6 0.42];
eR = -eL;
cLs = aL + bL + eL * theta.Lh / 3;                      % c_L(Lh)
bR = cLs - aR - eR * theta.Lh / 3;
P = struct('aL', aL, 'bL', bL, 'eL', eL, 'aR', aR, 'bR', bR, 'eR', eR);
end

function [Y, YP] = eval_poly(theta, slope, x)
P = coeffs(theta, slope);
Lh = theta.Lh; t = x / Lh;
q = 3 * t.^2 - 2 * t.^3;  dq = 6 * (t - t.^2) / Lh;  ddq = 6 * (1 - 2 * t) / Lh^2;
r = x.^3 / (3 * Lh^2);    dr = x.^2 / Lh^2;          ddr = 2 * x / Lh^2;
Iq = x.^3 / Lh^2 - x.^4 / (2 * Lh^3);  Ir = x.^4 / (12 * Lh^2);
n = numel(x);
Y = zeros(24, n); YP = zeros(24, n);
dens = zeros(4, n);
for h = 1:2
    if h == 1, a = P.aL; b = P.bL; e = P.eL; off = 0; else, a = P.aR; b = P.bR; e = P.eR; off = 8; end
    I = zeros(4, n);
    for k = 1:4
        c = a(k) + b(k) * q + e(k) * r;
        dc = b(k) * dq + e(k) * dr;
        ddc = b(k) * ddq + e(k) * ddr;
        Y(off + 2*k - 1, :) = c;  YP(off + 2*k - 1, :) = dc;
        Y(off + 2*k, :) = dc;     YP(off + 2*k, :) = ddc;
        I(k, :) = a(k) * x + b(k) * Iq + e(k) * Ir;
    end
    dens(2*h - 1, :) = Y(off + 1, :) + Y(off + 3, :);   % R pool
    dens(2*h, :)     = Y(off + 5, :) + Y(off + 7, :);   % B pool
    acc = [I(1, :) + I(2, :); I(3, :) + I(4, :)];
    Y(20 + 2*h - 1, :) = acc(1, :);  Y(20 + 2*h, :) = acc(2, :);
end
YP(21, :) = dens(1, :); YP(22, :) = dens(2, :); YP(23, :) = dens(3, :); YP(24, :) = dens(4, :);
% shifts: sigma = -J/Lh, J = DM*Rm'(Lh) + DC*Rc'(Lh) per half pool
J = [theta.DM * P.eL(1) + theta.DC * P.eL(2), theta.DM * P.eL(3) + theta.DC * P.eL(4), ...
     theta.DM * P.eR(1) + theta.DC * P.eR(2), theta.DM * P.eR(3) + theta.DC * P.eR(4)];
sig = -J / Lh;                                          % RL BL RR BR
Y(17:20, :) = repmat(sig(:), 1, n);
end

function [sol, tgt] = poly(theta, slope, n)
if nargin < 3, n = 401; end
x = linspace(0, theta.Lh, n);
[Y, YP] = eval_poly(theta, slope, x);
sol = struct('solver', 'bvp4c', 'x', x, 'y', Y, 'yp', YP, ...
    'stats', struct('nmeshpoints', n, 'maxres', 1e-12, 'nODEevals', 0, 'nBCevals', 0));
tgt = struct('nRL', Y(21, end), 'nRR', Y(23, end), 'nBL', Y(22, end), 'nBR', Y(24, end));
end

function sol = uniform(theta, u, n)
if nargin < 4, n = 51; end
x = linspace(0, theta.Lh, n);
Y = zeros(24, n);
Y([1 3 5 7], :) = repmat(u(:), 1, n); Y([9 11 13 15], :) = repmat(u(:), 1, n);
Y(21, :) = (u(1) + u(2)) * x; Y(22, :) = (u(3) + u(4)) * x;
Y(23, :) = Y(21, :); Y(24, :) = Y(22, :);
F = sag_coupled24('ode', x, Y, theta);
sol = struct('solver', 'bvp4c', 'x', x, 'y', Y, 'yp', F, ...
    'stats', struct('nmeshpoints', n, 'maxres', 1e-12, 'nODEevals', 0, 'nBCevals', 0));
end
