function rows = sag_t20_identities(E)
%SAG_T20_IDENTITIES  Packet 2C section 2c identities.
%   (a) Algebraic, on the real right-hand side: DM*F(Rm'') + DC*F(Rc'') =
%       -sigma for each half pool, on random states (the reaction cancels).
%   (b) On a synthetic Hermite-exact candidate with nonzero seam slopes and
%       shifts sigma = -J/Lh (sag_test_fixture 'poly'), the diagnostics give
%       antisymmetry 0, velocity identity at round-off, the two velocity
%       readings equal, and myxo's reading (4/L) J = -2 sigma.
%   Deviation from the packet wording ("linear-reaction stub"): no exact
%   solution of a linear stub BVP is constructed; (a) is the identity the
%   integration rests on, (b) exercises the diagnostics that measure it.
T = 'T20 identities';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
M = sa_model('guzzo_noxa');
th = M.normalize(sag_test_fixture('params', 1));
rs = RandStream('mt19937ar', 'Seed', 20);
Y = 0.01 + 2 * rand(rs, 24, 2000); Y(17:20, :) = randn(rs, 4, 2000);
F = sag_coupled24('ode', 0, Y, th);
res = [th.DM * F(2, :) + th.DC * F(4, :) + Y(17, :);  th.DM * F(6, :) + th.DC * F(8, :) + Y(18, :);
       th.DM * F(10, :) + th.DC * F(12, :) + Y(19, :); th.DM * F(14, :) + th.DC * F(16, :) + Y(20, :)];
scale = abs(Y(17:20, :)) + 1;
rows = sa_check(rows, T, '(a) DM*F(Rm'''') + DC*F(Rc'''') + sigma = 0 per half pool on 2000 random states', ...
    max(abs(res(:)) ./ scale(:)) <= 1e-12, sprintf('max %.3g', max(abs(res(:)))));
[sol, tg] = sag_test_fixture('poly', th, 1);
D = sag_coupled24('diagnostics', sol, th, tg);
rows = sa_check(rows, T, '(b) synthetic: sigma_L = -sigma_R exactly, |J + Lh sigma| at round-off, readings equal', ...
    D.antisym_max <= 1e-15 && D.vel_discrepancy_max <= 1e-14 * max(1, D.J_max) && ...
    max(abs(D.vel_from_slope - D.vel_from_shift)) <= 1e-14 && D.s_max > 1e-3, ...
    sprintf('s_max %.3g antisym %.3g vel %.3g slope %s shift %s', D.s_max, D.antisym_max, ...
    D.vel_discrepancy_max, mat2str(D.vel_from_slope, 6), mat2str(D.vel_from_shift, 6)));
rows = sa_check(rows, T, '(b) myxo reading (4/L) J_L equals -2 sigma_L', ...
    max(abs((4 / th.L) * D.J(1:2) + 2 * D.sigma(1:2))) <= 1e-14, '');
end
