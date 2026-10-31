function P = sap_post(main_variable)
%SAP_POST  Post-processing of the original five-field model (plots, traces).
%
%   P = SAP_POST() / SAP_POST(main_variable) is the post struct
%   sa_post_model('original_polarity', main_variable) returns. Coordinates
%   (plot order): x = dG, y = dP from the byte copy sca_v_obs.m
%   (figures/fig1/code), z = dn2 = int_0^1 cos(2 pi eta) m d(eta) of the
%   main variable m (a '+'-sum of u q v w p, sa_main_variable; default u),
%   with sca_v_obs's mesh and trapezoid: for 'u+q' it is sca_v_obs's dn2
%   bit for bit (Fig 3's coordinate, the default before software 0.10).
%   The full-PDE operator is the byte copy fpeq_model.m
%   (pde_rhs, pde_jacobian) with the run's own theta; the integrator is the
%   Fig 3 trace loop (ode15s, analytic sparse Jacobian, Refine 1, output
%   function with a stationarity stop at max|dY/dt| <= 1e-8 and a mass-drift
%   stop at 1e-6 relative). The output function also prints a progress
%   status line when one is due (sa_progress; screen only).
%
%   Half-domain convention: xi in [0, 1/2], xi = 0 the outer pole of each
%   half, xi = 1/2 the seam; the left half maps to eta = xi, the right half to
%   eta = 1 - xi (pchip on a uniform mesh with a node at eta = 1/2, the seam
%   node taken from the left half), as figures/fig1/code/sca_v_ic.m does.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(main_variable), main_variable = 'u'; end
P = struct();
P.id = 'original_polarity';
P.field_names = {'u', 'q', 'v', 'w', 'p'};
v = sa_main_variable('parse', main_variable, P.field_names, 'main_variable');
P.main_variable = v.text;
P.legacy_main_variable = 'u';           % the label's variable before software 0.10 (rows 1 / 11)
P.legacy_dn2_variable = 'u+q';          % dn2's variable before software 0.10 (sca_v_obs)
P.n_fields = 5;
P.sol_rows = [1 3 5 7 9; 11 13 15 17 19];
P.profile_rows = [1:5; 6:10];           % identity-profile rows of the left / right half
P.coord_labels = {'\Delta G', '\Delta P', sprintf('\\Delta n_2 (%s)', v.text)};
P.coord_names = {'dG', 'dP', 'dn2'};
P.mass_names = {'M_G', 'M_P'};
P.umax_names = {'umaxL', 'umaxR'};
P.lattice_axes = [2 1];                 % dG is lattice axis 2, dP axis 1
P.mass_pairs = [1 3];                   % (u, v): both in the G pool only
P.full_profile = @full_profile;
P.observables = @(theta, x, Y) observables(theta, x, Y, v.index);
P.trace_ic = @trace_ic;
P.integrate = @(theta, Y0, entry) integrate(theta, Y0, entry, v.index);
P.stop_rule = @stop_rule;
P.entry_defaults = @entry_defaults;
P.plot_nx = 1001;
P.aspect = @(c) [1 1/3.6 1/1.4];        % Fig 3b display stretch dG : 3.6 dP : 1.4 dn2
P.files = {fullfile(here, 'sap_post.m'), fullfile(here, 'fpeq_model.m'), fullfile(here, 'sca_v_obs.m')};
P.stops = {'horizon_reached', 'stationary', 'mass_drift_exceeded'};
P.drift_warn = Inf;                     % the integrator stops above 1e-6 instead
P.drift_warn_text = '';
end

% =========================================================================

function e = entry_defaults()
e = struct('perturbation', struct('kind', 'none'), 't_out', [0:0.05:50, 51:1:5000], ...
    'nx', 1001, 'profile_every', 10, 'reltol', 1e-8, 'abstol', 1e-10);
end

function F = full_profile(theta, xi, YL, YR, nx) %#ok<INUSL>
if mod(nx - 1, 2) ~= 0
    error('SAP_POST:NX', 'nx - 1 must be even so a node sits at eta = 1/2 (got nx = %d).', nx);
end
xi = xi(:).';
if ~(abs(xi(1)) < 1e-14 && abs(xi(end) - 0.5) < 1e-14)
    error('SAP_POST:DOMAIN', 'the half-domain mesh is not [0, 1/2].');
end
eta = linspace(0, 1, nx);
im  = (nx+1)/2;
Y = zeros(size(YL, 1), nx);
for sp = 1:size(YL, 1)
    Y(sp, 1:im)     = interp1(xi, YL(sp,:), eta(1:im),     'pchip');
    Y(sp, im+1:end) = interp1(xi, YR(sp,:), 1-eta(im+1:end),'pchip');
end
F = struct('x', eta, 'eta', eta, 'Y', Y, 'im', im);
end

function O = observables(theta, x, Y, mv) %#ok<INUSL>
S = sca_v_obs(x, Y);
im = (numel(x) + 1) / 2;
O = struct('c', [S.d_G S.d_P dn2(x, Y, mv)], 'masses', [S.M_G S.M_P], ...
    'umax', [max(Y(1, 1:im)) max(Y(1, im:end))]);
end

function z = dn2(eta, Y, mv)
% sca_v_obs's dn2 expression with the main variable's rows mv in place of
% u + q (rows [1 2] give its value bit for bit)
z = trapz(eta, cos(2*pi*eta) .* sa_main_variable('sum', Y, mv));
end

function ic = trace_ic(theta, sol, entry)
xi = sol.x(:).';
YL = sol.y([1 3 5 7 9], :);
YR = sol.y([11 13 15 17 19], :);
F = full_profile(theta, xi, YL, YR, entry.nx);
ic = struct('x', F.x, 'Y_root', F.Y, 'Y', perturb(F.Y, F.eta, entry.perturbation, [1 3]));
end

function Y = perturb(Y, eta, pert, pairs)
if strcmp(pert.kind, 'none'), return; end
a = pert.amplitude; m = pert.mode;
for k = 1:size(pairs, 1)
    i = pairs(k, 1); j = pairs(k, 2);
    d = a * min(Y(i, :), Y(j, :)) .* cos(m * pi * eta);
    Y(i, :) = Y(i, :) + d;
    Y(j, :) = Y(j, :) - d;
end
end

function [k, d] = stop_rule(m, M0, r)
% Fig 3's stops, checked after every output frame: the relative drift of
% either pool mass above 1e-6, then stationarity max|dY/dt| <= 1e-8
k = ''; d = '';
dG = abs(m(1)-M0(1))/M0(1); dP = abs(m(2)-M0(2))/M0(2);
if dG > 1e-6 || dP > 1e-6
    k = 'mass_drift_exceeded';
    d = sprintf('relative drift %.3e (G) %.3e (P) above 1e-6', dG, dP);
elseif r <= 1e-8
    k = 'stationary';
    d = sprintf('max|dY/dt| = %.3e <= 1e-8', r);
end
end

function fr = integrate(theta, Y0, entry, mv)
th = theta;
nx = size(Y0, 2);
eta = linspace(0, 1, nx);
im = (nx + 1) / 2;
every = entry.profile_every;
S = struct('t', [], 'coord', zeros(0, 3), 'masses', zeros(0, 2), 'umax', zeros(0, 2), ...
    'rhs_inf', [], 'profile_frames', [], 'stop', '', 'stop_detail', '');
prof = {};
zlast = [];
push(0, Y0);
M0 = S.masses(1, :);
w0 = tic;
opts = odeset('RelTol', entry.reltol, 'AbsTol', entry.abstol, ...
    'Jacobian', @(t,z) fpeq_model('pde_jacobian', reshape(z,nx,5).', th.L, th), ...
    'OutputFcn', @outfcn, 'Refine', 1);
ode15s(@(t,z) reshape(fpeq_model('pde_rhs', reshape(z,nx,5).', th.L, th).',[],1), ...
       entry.t_out, reshape(Y0.',[],1), opts);
if isempty(S.stop), S.stop = 'horizon_reached'; end
fr = S;
fr.wall_s = toc(w0);
fr.mass_drift_rel = abs(S.masses - M0) ./ M0;
fr.Y_final = zlast;
if isempty(fr.profile_frames) || fr.profile_frames(end) ~= numel(S.t)
    prof{end+1} = zlast;
    fr.profile_frames(end+1) = numel(S.t);
end
fr.profiles = cat(3, prof{:});
fr.x = eta;

    function push(t, Y)
        O = sca_v_obs(eta, Y);
        F = fpeq_model('pde_rhs', Y, th.L, th);
        S.t(end+1) = t;
        S.coord(end+1, :) = [O.d_G O.d_P dn2(eta, Y, mv)];
        S.masses(end+1, :) = [O.M_G O.M_P];
        S.umax(end+1, :) = [max(Y(1,1:im)) max(Y(1,im:end))];
        S.rhs_inf(end+1) = max(abs(F(:)));
        if mod(numel(S.t) - 1, every) == 0
            prof{end+1} = Y;
            S.profile_frames(end+1) = numel(S.t);
        end
        zlast = Y;
    end
    function st = outfcn(tv, zv, flag)
        st = 0;
        if ~isempty(flag), return; end
        for c = 1:numel(tv)
            if tv(c) <= S.t(end) + 1e-14, continue; end
            push(tv(c), reshape(zv(:,c), nx, 5).');
            if sa_progress('due')            % screen only: a status line during a long trace
                sa_progress('status', struct('within', sprintf('simulated t = %.6g of %.6g', tv(c), ...
                    entry.t_out(end))));
            end
            [k, d] = stop_rule(S.masses(end, :), M0, S.rhs_inf(end));
            if ~isempty(k)
                S.stop = k; S.stop_detail = d;
                st = 1; return
            end
        end
    end
end
