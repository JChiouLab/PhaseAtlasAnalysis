function rows = sag_t20b_linear_stub(E) %#ok<INUSD>
%SAG_T20B_LINEAR_STUB  Milestone-3 T20' (2C F-5): the noxa shift/velocity
%   identities on a SOLVED linear-stub BVP. The reaction is replaced by the
%   linear exchange s = k1*Rm - k2*Rc (per half pool), for which the coupled
%   24-state BVP (sag_coupled24 structure: DM Rm'' = -(s + sigma),
%   DC Rc'' = s, pole no-flux, seam continuity and flux matching,
%   accumulators) has the closed-form solution
%     u = DM Rm + DC Rc = u0 - sigma xi^2/2,  s = A cos(kappa xi) - k1 sigma/(DM kappa^2),
%     kappa^2 = k1/DM + k2/DC; right half sigma_R = -sigma_L, A_R = -A_L,
%     u0_R = u0_L - sigma_L Lh^2, A_L cos(kappa Lh) = k1 sigma_L/(DM kappa^2).
%   sigma is NOT built from J (unlike the 2C fixture, F-5): J + Lh sigma = 0
%   holds only because the profile solves the equations. Checked: the stub
%   ODE and the module's boundary conditions hold at round-off; the module
%   diagnostics then give antisymmetry 0 and |J + Lh sigma| at round-off, and
%   myxo's reading (4/L) J = -2 sigma.
T = 'T20'' linear-stub BVP';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
M = sa_model('guzzo_noxa');
th = M.normalize(sag_test_fixture('params', 1));
DM = th.DM; DC = th.DC; Lh = th.Lh;
k1 = 0.7; k2 = 0.3;
kap = sqrt(k1 / DM + k2 / DC);
den = k1 * DC + k2 * DM;
x = linspace(0, Lh, 401);
pools = struct('sig', {0.05, -0.02}, 'u0', {30, 25});   % R and B pools
y = zeros(24, numel(x)); yp = y;
for p = 1:2
    sL = pools(p).sig; sR = -sL;
    AL = k1 * sL / (DM * kap^2 * cos(kap * Lh)); AR = -AL;
    halves = {struct('sig', sL, 'A', AL, 'u0', pools(p).u0), ...
              struct('sig', sR, 'A', AR, 'u0', pools(p).u0 - sL * Lh^2)};
    for h = 1:2
        H = halves{h};
        u = H.u0 - H.sig * x.^2 / 2; du = -H.sig * x; ddu = -H.sig * ones(size(x));
        s = H.A * cos(kap * x) - k1 * H.sig / (DM * kap^2);
        ds = -kap * H.A * sin(kap * x); dds = -kap^2 * H.A * cos(kap * x);
        Rm = (DC * s + k2 * u) / den; Rc = (k1 * u - DM * s) / den;
        dRm = (DC * ds + k2 * du) / den; dRc = (k1 * du - DM * ds) / den;
        ddRm = (DC * dds + k2 * ddu) / den; ddRc = (k1 * ddu - DM * dds) / den;
        base = 8 * (h - 1) + 4 * (p - 1);           % rows Rm Rm' Rc Rc' of this half pool
        y(base + (1:4), :) = [Rm; dRm; Rc; dRc];
        yp(base + (1:4), :) = [dRm; ddRm; dRc; ddRc];
        y(16 + 2 * (h - 1) + p, :) = H.sig;           % sigma_RL sigma_BL sigma_RR sigma_BR
        acc = 20 + 2 * (h - 1) + p;                   % I_RL I_BL I_RR I_BR
        Iex = @(xx) arrayfun(@(t) integral(@(z) mass(z, H, k1, k2, DM, DC, kap, den), 0, t, ...
            'AbsTol', 1e-14, 'RelTol', 1e-13), xx);
        y(acc, :) = Iex(x);
        yp(acc, :) = Rm + Rc;
        res_ode(h, p) = max(abs([DM * ddRm + (k1 * Rm - k2 * Rc) + H.sig, DC * ddRc - (k1 * Rm - k2 * Rc)])); %#ok<AGROW>
    end
end
sol = struct('solver', 'bvp4c', 'x', x, 'y', y, 'yp', yp);
tg = struct('nRL', y(21, end), 'nBL', y(22, end), 'nRR', y(23, end), 'nBR', y(24, end));
r = sag_coupled24('bc', y(:, 1), y(:, end), th, tg);
rows = sa_check(rows, T, 'solved: linear-stub ODE residual and the module''s 24 boundary residuals at round-off', ...
    max(res_ode(:)) <= 1e-11 && max(abs(r)) <= 1e-9, sprintf('ODE %.3g; BC %.3g', max(res_ode(:)), max(abs(r))));
D = sag_coupled24('diagnostics', sol, th, tg);
rows = sa_check(rows, T, 'diagnostics on the solved BVP: antisymmetry 0, |J + Lh sigma| at round-off, readings equal', ...
    D.antisym_max <= 1e-15 && D.vel_discrepancy_max <= 1e-10 * max(1, D.J_max) && D.s_max > 1e-3 && ...
    max(abs(D.vel_from_slope - D.vel_from_shift)) <= 1e-10, ...
    sprintf('s_max %.3g antisym %.3g |J+Lh s| %.3g J %s', D.s_max, D.antisym_max, D.vel_discrepancy_max, mat2str(D.J, 6)));
rows = sa_check(rows, T, 'myxo reading (4/L) J_L = -2 sigma_L', max(abs((4 / th.L) * D.J(1:2) + 2 * D.sigma(1:2))) <= 1e-10, '');
end

function m = mass(z, H, k1, k2, DM, DC, kap, den)
u = H.u0 - H.sig * z.^2 / 2;
s = H.A * cos(kap * z) - k1 * H.sig / (DM * kap^2);
m = ((DC * s + k2 * u) + (k1 * u - DM * s)) / den;
end
