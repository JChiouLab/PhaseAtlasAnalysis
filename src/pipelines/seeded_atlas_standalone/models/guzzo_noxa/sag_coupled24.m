function varargout = sag_coupled24(action, varargin)
%SAG_COUPLED24  Guzzo_noXA coupled forced paired-half BVP, 24 states (unclipped).
%
%   Both halves run pole -> seam on
%   xi in [0, Lh], Lh = L/2: Y_L(xi) = Y(x = xi), Y_R(xi) = Y(x = L - xi)
%   (myxo bvp24_noxa_coupled.m:18-22).
%
%   State (24 rows)
%     1:8    left half   Rm Rm' Rc Rc' Bm Bm' Bc Bc'
%     9:16   right half  same order
%     17:20  shifts      sigma_RL sigma_BL sigma_RR sigma_BR  (sigma' = 0)
%     21:24  accumulators I_RL I_BL I_RR I_BR
%   Rows 1:16 are myxo's y1..y16, rows 21:24 myxo's y17..y20 and rows 17:20
%   myxo's bvp4c parameters q = [sRL sBL sRR sBR]. Reactions are gnx_reaction
%   (unclipped); the shift enters the MEMBRANE species of each pool
%   (bvp24_noxa_coupled.m:55-63).
%
%   Actions
%     F = ode(x, Y, theta)                24 x N right-hand side (vectorized)
%     r = bc(ya, yb, theta, tgt)          24 boundary residuals, myxo's order
%     t = targets(theta, dn_R, dn_B)      half-pool targets nRL nRR nBL nBR
%     g = pair_guess(xg, YL4, YR4)        24-row cold guess (myxo's builder)
%     D = diagnostics(sol, theta, tgt)    seam / mass / shift / velocity values
%     L = layout()
%     n = guard_count() ; guard_reset()   singularity-guard hits
%
%   Singularity guard (as in the uncoupled solver): if any evaluated Bm_L or Bm_R is
%   <= -K/2 the whole right-hand side is NaN and the hit is counted.

persistent n_guard
if isempty(n_guard), n_guard = 0; end
switch action
    case 'ode'
        [F, hit] = ode(varargin{:});
        n_guard = n_guard + hit;
        varargout{1} = F;
    case 'bc',          varargout{1} = bc(varargin{:});
    case 'targets',     varargout{1} = targets(varargin{:});
    case 'pair_guess',  varargout{1} = pair_guess(varargin{:});
    case 'diagnostics', varargout{1} = diagnostics(varargin{:});
    case 'layout',      varargout{1} = layout();
    case 'guard_count', varargout{1} = n_guard;
    case 'guard_reset', n_guard = 0; varargout{1} = 0;
    otherwise, error('SAG:COUPLED24', 'Unknown sag_coupled24 action %s.', action);
end
end

% =========================================================================

function L = layout()
L = struct('n_states', 24, ...
    'left', 1:8, 'right', 9:16, ...
    'conc_left', [1 3 5 7], 'conc_right', [9 11 13 15], ...
    'deriv_left', [2 4 6 8], 'deriv_right', [10 12 14 16], ...
    'sigma', 17:20, 'sigma_names', {{'sigma_RL', 'sigma_BL', 'sigma_RR', 'sigma_BR'}}, ...
    'acc', 21:24, 'acc_names', {{'I_RL', 'I_BL', 'I_RR', 'I_BR'}}, ...
    'species', {{'Rm', 'Rc', 'Bm', 'Bc'}}, ...
    'myxo_rows', [1:16, 21:24], ...
    'convention', ['xi in [0, Lh] pole -> seam; Y_L(xi) = Y(x = xi); ' ...
                   'Y_R(xi) = Y(x = L - xi); dimensional um']);
end

% -------------------------------------------------------------------------
% The corrector's right-hand side. Unclipped: no clamp of any state may
% appear in this subfunction (source scan).
function [F, hit] = ode(~, Y, theta)
hit = 0;
if any(Y(5, :) <= -0.5 * theta.K) || any(Y(13, :) <= -0.5 * theta.K)
    hit = 1;
    F = nan(24, size(Y, 2));
    return
end
[sRL, sBL] = gnx_reaction(Y([1 3 5 7], :), theta);
[sRR, sBR] = gnx_reaction(Y([9 11 13 15], :), theta);
z = zeros(1, size(Y, 2));
F = [Y(2, :);  -(sRL + Y(17, :)) / theta.DM;
     Y(4, :);  -(-sRL)           / theta.DC;
     Y(6, :);  -(sBL + Y(18, :)) / theta.DM;
     Y(8, :);  -(-sBL)           / theta.DC;
     Y(10, :); -(sRR + Y(19, :)) / theta.DM;
     Y(12, :); -(-sRR)           / theta.DC;
     Y(14, :); -(sBR + Y(20, :)) / theta.DM;
     Y(16, :); -(-sBR)           / theta.DC;
     z; z; z; z;
     Y(1, :) + Y(3, :);  Y(5, :) + Y(7, :);
     Y(9, :) + Y(11, :); Y(13, :) + Y(15, :)];
end

% -------------------------------------------------------------------------
function r = bc(ya, yb, theta, tgt)
% 24 rows in myxo's order (bvp24_noxa_coupled.m:81-93): 8 pole no-flux,
% 4 seam continuity, 4 seam flux matches, 4 accumulator starts, 4 targets.
r = [ya(2);  ya(4);  ya(6);  ya(8);
     ya(10); ya(12); ya(14); ya(16);
     yb(1) - yb(9);
     yb(3) - yb(11);
     yb(5) - yb(13);
     yb(7) - yb(15);
     theta.DM * (yb(2) + yb(10));
     theta.DC * (yb(4) + yb(12));
     theta.DM * (yb(6) + yb(14));
     theta.DC * (yb(8) + yb(16));
     ya(21); ya(22); ya(23); ya(24);
     yb(21) - tgt.nRL;  yb(22) - tgt.nBL;
     yb(23) - tgt.nRR;  yb(24) - tgt.nBR];
end

% -------------------------------------------------------------------------
function t = targets(theta, dnR, dnB)
% (rho + dn/2)*Lh on the left, (rho - dn/2)*Lh on the right
% (bvp24_noxa_coupled.m:96-102 at rho = 1; model page general form).
Lh = theta.Lh;
t = struct('nRL', (theta.rho_R + dnR / 2) * Lh, 'nRR', (theta.rho_R - dnR / 2) * Lh, ...
           'nBL', (theta.rho_B + dnB / 2) * Lh, 'nBR', (theta.rho_B - dnB / 2) * Lh);
end

% -------------------------------------------------------------------------
function g = pair_guess(xg, GL, GR)
% Myxo's guess builder (bvp24_noxa_coupled.m:116-128): concentration rows from
% the two halves (each pole -> seam), derivative rows by gradient on xg,
% accumulators cumsum(.)*dx, zero shifts (here states, there parameters).
xg = xg(:).';
n = numel(xg);
dx = xg(2) - xg(1);
y = zeros(24, n);
for k = 1:4
    y(2*k - 1, :)     = GL(k, :);
    y(2*k, :)         = gradient(GL(k, :), xg);
    y(8 + 2*k - 1, :) = GR(k, :);
    y(8 + 2*k, :)     = gradient(GR(k, :), xg);
end
y(21, :) = cumsum(GL(1, :) + GL(2, :)) * dx;
y(22, :) = cumsum(GL(3, :) + GL(4, :)) * dx;
y(23, :) = cumsum(GR(1, :) + GR(2, :)) * dx;
y(24, :) = cumsum(GR(3, :) + GR(4, :)) * dx;
g = struct('x', xg, 'y', y);
end

% -------------------------------------------------------------------------
function D = diagnostics(sol, theta, tgt)
% Values of the module checks, recomputed from the returned
% iterate (x, y, yp, the bvp4c interpolant). Never trusts a corrector flag.
x = sol.x; y = sol.y; yp = sol.yp;
ya = y(:, 1); yb = y(:, end);
D = struct();
D.n_mesh = numel(x);
D.seam_conc_mismatch = [yb(1) - yb(9), yb(3) - yb(11), yb(5) - yb(13), yb(7) - yb(15)];
D.seam_flux = [theta.DM * (yb(2) + yb(10)), theta.DC * (yb(4) + yb(12)), ...
               theta.DM * (yb(6) + yb(14)), theta.DC * (yb(8) + yb(16))];
n = [tgt.nRL tgt.nBL tgt.nRR tgt.nBR];                 % accumulator order 21:24
D.targets = n;
D.acc_end = yb(21:24).';
D.acc_start = ya(21:24).';
D.acc_rel = abs(D.acc_end - n) ./ n;
D.acc_start_max = max(abs(D.acc_start));
% exact integral of bvp4c's cubic Hermite interpolant of each half-pool density
rowsets = {[1 3], [5 7], [9 11], [13 15]};
I1 = zeros(1, 4); I2 = zeros(1, 4);
xm = (x(1:end-1) + x(2:end)) / 2;
[Ym, YPm] = deval(sol, xm);
xb = reshape([x(1:end-1); xm], 1, []); xb(end+1) = x(end);
for p = 1:4
    r = rowsets{p};
    v = sum(y(r, :), 1); s = sum(yp(r, :), 1);
    I1(p) = hermite_integral(x, v, s);
    vb = reshape([v(1:end-1); sum(Ym(r, :), 1)], 1, []); vb(end+1) = v(end); %#ok<AGROW>
    sb = reshape([s(1:end-1); sum(YPm(r, :), 1)], 1, []); sb(end+1) = s(end); %#ok<AGROW>
    I2(p) = hermite_integral(xb, vb, sb);
end
D.mass_dense = I1;
D.mass_dense_bisected = I2;
D.mass_dense_rel = abs(I1 - n) ./ n;
D.mass_refine_rel = abs(I2 - I1) ./ abs(I1);
% shifts and the shift-flux identities
s = ya(17:20).';
D.sigma = s;                                           % RL BL RR BR
D.s_max = max(abs(s));
D.antisym = [abs(s(1) + s(3)), abs(s(2) + s(4))];
D.antisym_max = max(D.antisym);
J = [theta.DM * yb(2) + theta.DC * yb(4), theta.DM * yb(6) + theta.DC * yb(8), ...
     theta.DM * yb(10) + theta.DC * yb(12), theta.DM * yb(14) + theta.DC * yb(16)];
D.J = J;                                               % J_RL J_BL J_RR J_BR
D.vel_identity = abs(J + theta.Lh * s);               % every half pool (the velocity
D.vel_discrepancy_max = max(D.vel_identity);           % check needs the right half too)
D.J_max = max(abs(J));
D.vel_from_slope = (2 / theta.Lh) * J(1:2);            % = (4/L) J, myxo velocity
D.vel_from_shift = -2 * s(1:2);
end

function I = hermite_integral(x, v, s)
% Exact integral of the C1 cubic Hermite interpolant with node values v and
% slopes s: sum over intervals of h/2 (v_k + v_k+1) + h^2/12 (s_k - s_k+1).
h = diff(x);
I = sum(h / 2 .* (v(1:end-1) + v(2:end)) + h.^2 / 12 .* (s(1:end-1) - s(2:end)));
end
