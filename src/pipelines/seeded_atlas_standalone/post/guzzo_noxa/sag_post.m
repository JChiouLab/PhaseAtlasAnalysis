function P = sag_post(main_variable)
%SAG_POST  Post-processing of the M. xanthus model (plots, traces).
%
%   P = SAG_POST() / SAG_POST(main_variable) is the post struct
%   sa_post_model('guzzo_noxa', main_variable) returns.
%   Coordinates (plot order): x = dn_R, y = dn_B (the grid axes: left-
%   minus-right half masses of R = Rm + Rc and B = Bm + Bc divided by the
%   half length Lh, trapezoid split at the middle node; fraction of the total
%   Delta n/M = dn/(2 rho)), z = dn2 = int_0^1
%   cos(2 pi eta) m d(eta) with eta = x / L (the paper's Delta n2 of Eq. 10 = L * dn2),
%   m the main variable (a '+'-sum
%   of Rm Rc Bm Bc, sa_main_variable; default Bm).
%   Full-model trace: pdepe (m = 0) with the structure of Guzzo_noXA.m (the copy
%   beside this file, kept as the reference: pdefun c = 1, f = D u_x,
%   zero-flux bcfun, its tolerances) and the source terms of gnx_reaction
%   (the model definition in src/models/guzzo_noxa), with the run's own
%   theta. The initial values are given at the mesh nodes exactly. No drift
%   or stationarity stop: pdepe runs to the horizon (stop horizon_reached),
%   or returns fewer times than requested (stop pdepe_incomplete, as
%   Guzzo_noXA.m guards). The masses are recorded at every frame.
%
%   Half-domain convention: x in [0, Lh] pole -> middle; Y_L(x) = Y(x),
%   Y_R(x) = Y(L - x); pchip on a uniform mesh with a node at x = Lh,
%   the middle node from the left half.

here = fileparts(mfilename('fullpath'));
core = fullfile(fileparts(fileparts(fileparts(fileparts(fileparts(here))))), 'src', 'models', 'guzzo_noxa');
if exist(core, 'dir') ~= 7
    core = fullfile(fileparts(fileparts(fileparts(fileparts(here)))), 'models', 'guzzo_noxa');
end
if ~contains([pathsep path pathsep], [pathsep core pathsep]), addpath(core); end
if nargin < 1 || isempty(main_variable), main_variable = 'Bm'; end
P = struct();
P.id = 'guzzo_noxa';
P.field_names = {'Rm', 'Rc', 'Bm', 'Bc'};
v = sa_main_variable('parse', main_variable, P.field_names, 'main_variable');
P.main_variable = v.text;
P.legacy_main_variable = 'Bm';          % the label's variable of older versions (rows 5 / 13)
P.legacy_dn2_variable = 'Bm';           % the third coordinate's variable of older versions
P.n_fields = 4;
P.sol_rows = [1 3 5 7; 9 11 13 15];
P.profile_rows = [1:4; 5:8];            % identity-profile rows of the left / right half
P.coord_labels = {'\Delta n_R', '\Delta n_B', sprintf('\\Delta n_2 (%s)', v.text)};
P.coord_names = {'dn_R', 'dn_B', ['cos_moment_' strrep(v.text, '+', '_')]};
P.mass_names = {'N_R', 'N_B'};
P.umax_names = {'BmmaxL', 'BmmaxR'};
P.lattice_axes = [1 2];
P.velocity = struct('records', {{'velocity_R_slope', 'velocity_B_slope'}}, 'sign', 1);
                                        % d(dn_R)/dt, d(dn_B)/dt (= -2 x the left-half sources), the sheet PDFs' flux
P.mass_pairs = [1 2; 3 4];              % (Rm, Rc) and (Bm, Bc)
P.full_profile = @full_profile;
P.observables = @(theta, x, Y) observables(theta, x, Y, v.index);
P.trace_ic = @trace_ic;
P.integrate = @(theta, Y0, entry) integrate(theta, Y0, entry, v.index);
P.entry_defaults = @entry_defaults;
P.plot_nx = 1201;
P.aspect = @aspect;
P.files = {fullfile(here, 'sag_post.m'), fullfile(here, 'Guzzo_noXA.m'), fullfile(core, 'gnx_reaction.m')};
P.stops = {'horizon_reached', 'pdepe_incomplete'};
P.drift_warn = 1e-4;                    % no drift stop: a WARNING in trace_summary.txt
P.drift_warn_text = ['pdepe does not conserve the trapezoid mass exactly; read the trace with ' ...
    'that in mind.'];
end

% =========================================================================

function e = entry_defaults()
e = struct('perturbation', struct('kind', 'none'), 't_out', linspace(0, 200, 601), ...
    'nx', 241, 'profile_every', 10, 'reltol', 1e-6, 'abstol', 1e-8);
end

function a = aspect(c)
% equal dn axes; the third axis scaled to the same display height
r = max([max(c(:, 1)) - min(c(:, 1)), max(c(:, 2)) - min(c(:, 2)), eps]);
z = max(max(c(:, 3)) - min(c(:, 3)), eps);
a = [1 1 z / r];
end

function F = full_profile(theta, xi, YL, YR, nx)
if mod(nx - 1, 2) ~= 0
    error('SAG_POST:NX', 'nx - 1 must be even so a node sits at the seam (got nx = %d).', nx);
end
L = theta.L;
xi = xi(:).';
x = linspace(0, L, nx);
im = (nx + 1) / 2;
Y = zeros(size(YL, 1), nx);
for sp = 1:size(YL, 1)
    Y(sp, 1:im)     = interp1(xi, YL(sp, :), x(1:im), 'pchip');
    Y(sp, im+1:end) = interp1(xi, YR(sp, :), L - x(im+1:end), 'pchip');
end
F = struct('x', x, 'eta', x / L, 'Y', Y, 'im', im);
end

function O = observables(theta, x, Y, mv)
% mv: the main variable's rows (Bm alone, [3], gives the older value bit for bit)
L = theta.L; Lh = L / 2;
im = (numel(x) + 1) / 2;
R = Y(1, :) + Y(2, :); B = Y(3, :) + Y(4, :);
nRL = trapz(x(1:im), R(1:im)); nRR = trapz(x(im:end), R(im:end));
nBL = trapz(x(1:im), B(1:im)); nBR = trapz(x(im:end), B(im:end));
eta = x / L;
O = struct('c', [(nRL - nRR) / Lh, (nBL - nBR) / Lh, trapz(eta, cos(2 * pi * eta) .* sa_main_variable('sum', Y, mv))], ...
    'masses', [nRL + nRR, nBL + nBR], 'umax', [max(Y(3, 1:im)) max(Y(3, im:end))]);
end

function ic = trace_ic(theta, sol, entry)
xi = sol.x(:).';
YL = sol.y([1 3 5 7], :);
YR = sol.y([9 11 13 15], :);
F = full_profile(theta, xi, YL, YR, entry.nx);
ic = struct('x', F.x, 'Y_root', F.Y, 'Y', perturb(F.Y, F.eta, entry.perturbation, [1 2; 3 4]));
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

function fr = integrate(theta, Y0, entry, mv)
L = theta.L;
nx = size(Y0, 2);
x = linspace(0, L, nx);
t = entry.t_out;
D = [theta.DM; theta.DC; theta.DM; theta.DC];
w0 = tic;
opts = odeset('RelTol', entry.reltol, 'AbsTol', entry.abstol);
sol = pdepe(0, @(xq, tq, u, dudx) pdefun(u, dudx, D, theta), @(xq) Y0(:, node_of(x, xq)), ...
    @bcfun, x, t, opts);
nt = size(sol, 1);
stop = 'horizon_reached'; detail = '';
if nt < numel(t)
    stop = 'pdepe_incomplete';
    detail = sprintf('pdepe returned %d of %d output times (last t = %.6g)', nt, numel(t), t(nt));
end
fr = struct('t', t(1:nt), 'coord', zeros(nt, 3), 'masses', zeros(nt, 2), 'umax', zeros(nt, 2), ...
    'rhs_inf', nan(1, nt), 'profile_frames', [], 'stop', stop, 'stop_detail', detail);
prof = {};
for k = 1:nt
    Y = reshape(sol(k, :, :), nx, 4).';
    O = observables(theta, x, Y, mv);
    fr.coord(k, :) = O.c; fr.masses(k, :) = O.masses; fr.umax(k, :) = O.umax;
    if mod(k - 1, entry.profile_every) == 0 || k == nt
        prof{end+1} = Y; %#ok<AGROW>
        fr.profile_frames(end+1) = k;
    end
end
fr.wall_s = toc(w0);
fr.mass_drift_rel = abs(fr.masses - fr.masses(1, :)) ./ fr.masses(1, :);
fr.Y_final = reshape(sol(nt, :, :), nx, 4).';
fr.profiles = cat(3, prof{:});
fr.x = x;
end

function k = node_of(x, xq)
k = find(x == xq, 1);
if isempty(k)
    error('SAG_POST:IC', 'pdepe asked for an initial value off the mesh (x = %.17g).', xq);
end
end

function [c, f, s] = pdefun(u, dudx, D, theta)
c = ones(4, 1);
f = D .* dudx;
[sR, sB] = gnx_reaction(u, theta);
s = [ sR; -sR; ...
      sB; -sB ];
end

function [pl, ql, pr, qr] = bcfun(~, ~, ~, ~, ~)
pl = zeros(4,1); ql = ones(4,1);
pr = zeros(4,1); qr = ones(4,1);
end
