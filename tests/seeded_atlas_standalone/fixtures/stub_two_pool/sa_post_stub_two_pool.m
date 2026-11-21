function P = sa_post_stub_two_pool()
%SA_POST_STUB_TWO_POOL  Post module of the stub model, used by the stub-model tests.
%   Two half fields per side: R (c1 left, c2 right) and B (c3 left, c4
%   right); the stub's half domain [0, Lh] maps to x = xi (left) and
%   x = 2 Lh - xi (right). Coordinates: dn_R = (int_left R - int_right R)/Lh,
%   dn_B likewise, z = int_0^1 cos(2 pi eta) B d(eta). The integrator is a
%   stub that never time-steps (opts.integrator replaces it in the tests).
P = struct();
P.id = 'stub_two_pool';
P.field_names = {'R', 'B'};
P.n_fields = 2;
P.sol_rows = [1 5; 3 7];
P.profile_rows = [1 3; 2 4];
P.coord_labels = {'dn_R', 'dn_B', 'cos moment of B'};
P.coord_names = {'dn_R', 'dn_B', 'zB'};
P.mass_names = {'N_R', 'N_B'};
P.umax_names = {'RmaxL', 'RmaxR'};
P.lattice_axes = [1 2];
P.mass_pairs = zeros(0, 2);
P.full_profile = @full_profile;
P.observables = @observables;
P.trace_ic = @trace_ic;
P.integrate = @(theta, Y0, entry) error('SA:TEST:INTEGRATE', 'stub post: inject opts.integrator');
P.entry_defaults = @() struct('perturbation', struct('kind', 'none'), 't_out', 0:1:20, 'nx', 101, ...
    'profile_every', 5, 'reltol', 1e-6, 'abstol', 1e-8);
P.plot_nx = 201;
P.aspect = @(c) [1 1 1];
P.files = {mfilename('fullpath')};
P.files{1} = [P.files{1} '.m'];
P.stops = {'horizon_reached'};
P.drift_warn = 1e-4;
P.drift_warn_text = 'stub drift warning';
end

function F = full_profile(theta, xi, YL, YR, nx)
Lh = theta.Lh; L = 2 * Lh;
xi = xi(:).';
x = linspace(0, L, nx); im = (nx + 1) / 2;
Y = zeros(size(YL, 1), nx);
for k = 1:size(YL, 1)
    Y(k, 1:im) = interp1(xi, YL(k, :), x(1:im), 'pchip');
    Y(k, im+1:end) = interp1(xi, YR(k, :), L - x(im+1:end), 'pchip');
end
F = struct('x', x, 'eta', x / L, 'Y', Y, 'im', im);
end

function O = observables(theta, x, Y)
Lh = theta.Lh; im = (numel(x) + 1) / 2;
rl = trapz(x(1:im), Y(1, 1:im)); rr = trapz(x(im:end), Y(1, im:end));
bl = trapz(x(1:im), Y(2, 1:im)); br = trapz(x(im:end), Y(2, im:end));
eta = x / (2 * Lh);
O = struct('c', [(rl - rr) / Lh, (bl - br) / Lh, trapz(eta, cos(2 * pi * eta) .* Y(2, :))], ...
    'masses', [rl + rr, bl + br], 'umax', [max(Y(1, 1:im)) max(Y(1, im:end))]);
end

function ic = trace_ic(theta, sol, entry)
F = full_profile(theta, sol.x, sol.y([1 5], :), sol.y([3 7], :), entry.nx);
ic = struct('x', F.x, 'Y_root', F.Y, 'Y', F.Y);
end
