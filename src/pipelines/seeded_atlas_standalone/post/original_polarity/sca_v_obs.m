function O = sca_v_obs(eta, Y)
%SCA_V_OBS  Reduced observables of a full-domain state on a UNIFORM mesh with a
%   node at eta = 1/2.
%
%     d_G, d_P : left-minus-right half masses, shared-node symmetric
%                trapezoid split at eta = 1/2 -- the same convention
%                as the atlas, so the trace is measured in the
%                coordinates the atlas is drawn in.
%     M_G, M_P : totals, composite trapezoid (the quantity the semidiscrete
%                operator conserves exactly).
%     dn2      : Delta-n_2 = int_0^1 cos(2*pi*eta) * (u + q) d(eta),
%                the shape coordinate, computed here on u + q (an older
%                convention; the paper's Fig. 3 uses u alone). The paper's
%                Delta n2 = M_BG * dn2 (total GEF = 1 in this frame),
%                written in the coupled dimensionless frame (x = L*eta).

nx = numel(eta);
im = (nx+1)/2;
assert(abs(eta(im) - 0.5) < 1e-12, 'sca_v_obs:seam','no node at eta = 1/2');

G  = sum(Y(1:4,:), 1);      % u + q + v + w
Pm = Y(2,:) + Y(5,:);       % q + p

gL = trapz(eta(1:im),  G(1:im));   gR = trapz(eta(im:end),  G(im:end));
pL = trapz(eta(1:im), Pm(1:im));   pR = trapz(eta(im:end), Pm(im:end));

O = struct();
O.d_G = gL - gR;   O.d_P = pL - pR;
O.M_G = gL + gR;   O.M_P = pL + pR;
O.dn2 = trapz(eta, cos(2*pi*eta) .* (Y(1,:) + Y(2,:)));
O.u_max = max(Y(1,:));  O.u_min = min(Y(1,:));
O.minimum = min(Y, [], 'all');
[~, ip] = max(Y(1,:));  O.u_peak_eta = eta(ip);
O.convention = ['shared-node symmetric trapezoid split at eta=1/2 ' ...
                '(the atlas split convention); dn2 = int_0^1 cos(2 pi eta) (u + q) d eta'];
end
