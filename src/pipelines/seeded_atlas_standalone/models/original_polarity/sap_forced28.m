function varargout = sap_forced28(action, varargin)
%CHA_FORCED28  Forced 28-state paired-half quasi-steady BVP (unclipped).
%
%   Pole-to-seam xi in [0,H], H = 1/2. Each half carries its OWN shift pair and
%   its OWN mass accumulators (design section 5.3). State layout and the 28
%   boundary rows are frozen in DESIGN_LOCAL_P3.md sections 1-2.
%
%   Actions
%     ode(xi,y,th)                 -> 28x1 dy/dxi
%     bc(ya,yb,th,tgt)             -> 28x1 boundary residual (solver form)
%     bc_labelled(ya,yb,th,tgt)    -> .res (28x1), .label{28}, .group{28}
%     targets(th,d_P,d_G)          -> struct n_GL,n_PL,n_GR,n_PR (design 5.1)
%     from24(sol24,th)             -> 28-row guess from a Packet-2 24-state sol
%     pair_guess(xi,YL5,YR5,th)    -> 28-row guess from two half profiles
%     swap(sol)                    -> reflected 28-state (pure block swap)
%     swap_targets(tgt)            -> reflected targets
%     reconstruct(sol,neta)        -> full-domain [eta;Y]
%     diagnostics(sol,th,tgt)      -> seam/mass/shift/velocity diagnostics
%
%   Experiment-local; nothing here is canonical.

switch action
    case 'ode',          varargout{1} = ode28(varargin{:});
    case 'bc',           varargout{1} = bc28(varargin{:});
    case 'bc_labelled',  varargout{1} = bc28_labelled(varargin{:});
    case 'targets',      varargout{1} = targets(varargin{:});
    case 'from24',       varargout{1} = from24(varargin{:});
    case 'pair_guess',   varargout{1} = pair_guess(varargin{:});
    case 'swap',         varargout{1} = swapsol(varargin{:});
    case 'swap_targets', varargout{1} = swaptgt(varargin{:});
    case 'reconstruct',  [varargout{1:nargout}] = reconstruct(varargin{:});
    case 'diagnostics',  varargout{1} = diagnostics(varargin{:});
    case 'layout',       varargout{1} = layout();
    otherwise, error('sap_forced28:action','Unknown action %s', action);
end
end

% =========================================================================

function L = layout()
L = struct('n_states',28,'H',0.5, ...
    'left',1:10,'right',11:20, ...
    'conc_left',[1 3 5 7 9],'conc_right',[11 13 15 17 19], ...
    'dconc_left',[2 4 6 8 10],'dconc_right',[12 14 16 18 20], ...
    'sGL',21,'sPL',22,'sGR',23,'sPR',24, ...
    'IGL',25,'IPL',26,'IGR',27,'IPR',28, ...
    'species',{{'u','q','v','w','p'}}, ...
    'convention','xi in [0,H] pole-to-seam; Y_L(xi)=Y(eta=xi); Y_R(xi)=Y(eta=1-xi)');
end

function dy = ode28(~, y, th)
mu = sap_model('mu', th);
dy = zeros(28,1);
sh = [21 22; 23 24];          % [sG sP] per half
for hh = 1:2
    o  = (hh-1)*10;
    sG = y(sh(hh,1)); sP = y(sh(hh,2));
    u = y(o+1); q = y(o+3); v = y(o+5); w = y(o+7); p = y(o+9);
    Ru = (th.alpha*u^2 + 1)*v - th.beta*u - th.gamma*u*p;   % UNCLIPPED
    Rq =  th.gamma*u*p - th.delta*q;
    Rv = -(th.alpha*u^2 + 1)*v + th.beta*u + th.epsilon*w;
    Rw =  th.delta*q - th.epsilon*w;
    Rp = -th.gamma*u*p + th.delta*q;
    dy(o+1)  = y(o+2);   dy(o+2)  = -(Ru + sG)/mu(1);
    dy(o+3)  = y(o+4);   dy(o+4)  = -Rq/mu(2);
    dy(o+5)  = y(o+6);   dy(o+6)  = -Rv/mu(3);
    dy(o+7)  = y(o+8);   dy(o+8)  = -Rw/mu(4);
    dy(o+9)  = y(o+10);  dy(o+10) = -(Rp + sP)/mu(5);
end
dy(21:24) = 0;
dy(25) = y(1)+y(3)+y(5)+y(7);      % I_GL'
dy(26) = y(9)+y(3);                % I_PL'
dy(27) = y(11)+y(13)+y(15)+y(17);  % I_GR'
dy(28) = y(19)+y(13);              % I_PR'
end

function r = bc28(ya, yb, th, tgt)
B = bc28_labelled(ya, yb, th, tgt);
r = B.res;
end

function B = bc28_labelled(ya, yb, th, tgt)
mu = sap_model('mu', th);
sp = {'u','q','v','w','p'};
res = zeros(28,1); lab = cell(28,1); grp = cell(28,1);
k = 0;
for s = 1:5
    k=k+1; res(k)=ya(2*s);       lab{k}=sprintf('outer_noflux_L_%s',sp{s}); grp{k}='outer_left';
end
for s = 1:5
    k=k+1; res(k)=ya(10+2*s);    lab{k}=sprintf('outer_noflux_R_%s',sp{s}); grp{k}='outer_right';
end
for s = 1:5
    k=k+1; res(k)=yb(2*s-1)-yb(10+2*s-1); lab{k}=sprintf('seam_conc_%s',sp{s}); grp{k}='seam_conc';
end
for s = 1:5
    k=k+1; res(k)=mu(s)*(yb(2*s)+yb(10+2*s)); lab{k}=sprintf('seam_flux_%s',sp{s}); grp{k}='seam_flux';
end
k=k+1; res(k)=ya(25); lab{k}='I_GL(0)=0'; grp{k}='init_integral';
k=k+1; res(k)=ya(26); lab{k}='I_PL(0)=0'; grp{k}='init_integral';
k=k+1; res(k)=ya(27); lab{k}='I_GR(0)=0'; grp{k}='init_integral';
k=k+1; res(k)=ya(28); lab{k}='I_PR(0)=0'; grp{k}='init_integral';
k=k+1; res(k)=yb(25)-tgt.n_GL; lab{k}='I_GL(H)=n_GL'; grp{k}='half_mass';
k=k+1; res(k)=yb(26)-tgt.n_PL; lab{k}='I_PL(H)=n_PL'; grp{k}='half_mass';
k=k+1; res(k)=yb(27)-tgt.n_GR; lab{k}='I_GR(H)=n_GR'; grp{k}='half_mass';
k=k+1; res(k)=yb(28)-tgt.n_PR; lab{k}='I_PR(H)=n_PR'; grp{k}='half_mass';
assert(k==28,'sap_forced28:count','expected 28 rows, assembled %d',k);
B = struct('res',res,'label',{lab},'group',{grp},'n',k);
end

function t = targets(th, d_P, d_G)
% design section 5.1 with M_G = 1, M_P = rho
MG = 1; MP = th.rho;
t = struct('d_P', d_P, 'd_G', d_G, 'M_G', MG, 'M_P', MP, ...
    'n_GL', (MG + d_G)/2, 'n_GR', (MG - d_G)/2, ...
    'n_PL', (MP + d_P)/2, 'n_PR', (MP - d_P)/2, ...
    'derivation','DESIGN_LOCAL_P3 section 3: n_{c,L/R} = (M_c +/- d_c)/2');
end

function t = swaptgt(t0)
t = t0;
t.n_GL = t0.n_GR; t.n_GR = t0.n_GL;
t.n_PL = t0.n_PR; t.n_PR = t0.n_PL;
t.d_G  = -t0.d_G; t.d_P  = -t0.d_P;
t.derivation = [t0.derivation ' ; reflected (block swap of targets)'];
end

function g = from24(sol24, ~)
% Expand an accepted Packet-2 24-state solution to a 28-state guess: the single
% shared shift pair becomes two equal pairs, and the two global accumulators are
% replaced by four per-half accumulators rebuilt by cumtrapz of each half.
xi = sol24.x(:).'; n = numel(xi);
y = zeros(28, n);
y(1:20, :) = sol24.y(1:20, :);
y(21,:) = sol24.y(21,:); y(22,:) = sol24.y(22,:);
y(23,:) = sol24.y(21,:); y(24,:) = sol24.y(22,:);
y(25,:) = cumtrapz(xi, y(1,:)+y(3,:)+y(5,:)+y(7,:));
y(26,:) = cumtrapz(xi, y(9,:)+y(3,:));
y(27,:) = cumtrapz(xi, y(11,:)+y(13,:)+y(15,:)+y(17,:));
y(28,:) = cumtrapz(xi, y(19,:)+y(13,:));
g = struct('x', xi, 'y', y);
end

function g = pair_guess(xi, YL5, YR5, ~)
xi = xi(:).'; n = numel(xi);
y = zeros(28, n);
y([1 3 5 7 9], :)      = YL5;
y([11 13 15 17 19], :) = YR5;
for s = 1:5
    y(2*s, :)    = gradient(YL5(s,:), xi);
    y(10+2*s, :) = gradient(YR5(s,:), xi);
end
y(21:24,:) = 0;
y(25,:) = cumtrapz(xi, sum(YL5(1:4,:),1));
y(26,:) = cumtrapz(xi, YL5(5,:)+YL5(2,:));
y(27,:) = cumtrapz(xi, sum(YR5(1:4,:),1));
y(28,:) = cumtrapz(xi, YR5(5,:)+YR5(2,:));
g = struct('x', xi, 'y', y);
end

function g = swapsol(sol)
% Reflection eta -> 1-eta: PURE block swap (DESIGN_LOCAL_P3 section 5).
g = struct('x', sol.x(:).', 'y', sol.y);
g.y(1:10,:)  = sol.y(11:20,:);
g.y(11:20,:) = sol.y(1:10,:);
g.y(21:22,:) = sol.y(23:24,:);
g.y(23:24,:) = sol.y(21:22,:);
g.y(25:26,:) = sol.y(27:28,:);
g.y(27:28,:) = sol.y(25:26,:);
end

function [eta, Yf] = reconstruct(sol, neta)
if nargin < 2 || isempty(neta), neta = 601; end
H = 0.5;
eta = linspace(0,1,neta);
Yf = zeros(5, neta);
iL = eta <= H; iR = ~iL;
cl = @(z) min(max(z, sol.x(1)), sol.x(end));
YL = deval(sol, cl(eta(iL)));
YR = deval(sol, cl(1 - eta(iR)));
Yf(:, iL) = YL([1 3 5 7 9], :);
Yf(:, iR) = YR([11 13 15 17 19], :);
end

function D = diagnostics(sol, th, tgt)
mu = sap_model('mu', th);
H  = 0.5;
xi = sol.x(:).';
ya = sol.y(:,1); yb = sol.y(:,end);
B = bc28_labelled(ya, yb, th, tgt);

D = struct();
D.n_mesh = numel(xi);
D.bc_residual = B.res; D.bc_label = {B.label}; D.bc_group = {B.group};
D.bc_max_abs = max(abs(B.res));
for g = {'outer_left','outer_right','seam_conc','seam_flux','init_integral','half_mass'}
    m = strcmp(B.group, g{1});
    D.(['max_' g{1}]) = max(abs(B.res(m)));
end

% ---- shifts ------------------------------------------------------------
D.s_GL = sol.y(21,1); D.s_PL = sol.y(22,1);
D.s_GR = sol.y(23,1); D.s_PR = sol.y(24,1);
D.s_all = [D.s_GL D.s_PL D.s_GR D.s_PR];
D.s_max = max(abs(D.s_all));
D.s_constant_check = max(max(abs(sol.y(21:24,:) - sol.y(21:24,1))));
D.antisym_G = abs(D.s_GL + D.s_GR);
D.antisym_P = abs(D.s_PL + D.s_PR);
D.antisym_max = max(D.antisym_G, D.antisym_P);
D.antisym_rel = D.antisym_max / max(D.s_max, realmin);

% ---- seam --------------------------------------------------------------
D.seam_left_values  = sol.y([1 3 5 7 9], end).';
D.seam_right_values = sol.y([11 13 15 17 19], end).';
D.seam_conc_mismatch = D.seam_left_values - D.seam_right_values;
D.seam_dY_left  = sol.y([2 4 6 8 10], end).';
D.seam_dY_right = sol.y([12 14 16 18 20], end).';
D.j_from_left  = -(mu(:).') .* D.seam_dY_left;    % j = -mu Y_L'(H)
D.j_from_right =  (mu(:).') .* D.seam_dY_right;   % j = +mu Y_R'(H)
D.j_mismatch = D.j_from_left - D.j_from_right;
cG = [1 1 1 1 0]; cP = [0 1 0 0 1];
D.J_G = cG * D.j_from_left(:);
D.J_P = cP * D.j_from_left(:);

% ---- integrated shift/flux identities (design section 6) ---------------
% J_G = H s_GL, J_P = H s_PL ; and per-half c^T mu Y_h'(H) = -H s_{c,h}
D.identity_JG = D.J_G - H*D.s_GL;
D.identity_JP = D.J_P - H*D.s_PL;
D.cG_mu_dYL = cG * (mu(:).*D.seam_dY_left(:));
D.cP_mu_dYL = cP * (mu(:).*D.seam_dY_left(:));
D.cG_mu_dYR = cG * (mu(:).*D.seam_dY_right(:));
D.cP_mu_dYR = cP * (mu(:).*D.seam_dY_right(:));
D.identity_half_GL = D.cG_mu_dYL + H*D.s_GL;
D.identity_half_PL = D.cP_mu_dYL + H*D.s_PL;
D.identity_half_GR = D.cG_mu_dYR + H*D.s_GR;
D.identity_half_PR = D.cP_mu_dYR + H*D.s_PR;
D.identity_max = max(abs([D.identity_JG D.identity_JP ...
    D.identity_half_GL D.identity_half_PL D.identity_half_GR D.identity_half_PR]));

% ---- projected velocity, two independent routes ------------------------
D.vel_flux  = [-2*D.J_P, -2*D.J_G];        % [ddt d_P, ddt d_G]
D.vel_shift = [-2*H*D.s_PL, -2*H*D.s_GL];
D.vel_discrepancy = abs(D.vel_flux - D.vel_shift);
D.vel_max = max(abs(D.vel_flux));
D.vel_discrepancy_max = max(D.vel_discrepancy);
D.vel_discrepancy_rel = D.vel_discrepancy_max / max(D.vel_max, realmin);

% ---- masses: G6a accumulators (authoritative), G6b dense, G6c native ---
D.acc_endpoints = [yb(25) yb(26) yb(27) yb(28)];
D.targets = [tgt.n_GL tgt.n_PL tgt.n_GR tgt.n_PR];
D.G6a_residual = max(abs(D.acc_endpoints - D.targets));
D.d_G = yb(25) - yb(27);            % authoritative: accumulator endpoints
D.d_P = yb(26) - yb(28);
D.M_G = yb(25) + yb(27);
D.M_P = yb(26) + yb(28);

[D.dense_20001, D.dense_40001, D.G6b_residual, D.G6b_refine] = dense_masses(sol, D.targets);
D.G6b_ndense = [20001 40001];

nl = [trapz(xi, sum(sol.y([1 3 5 7],:),1)), trapz(xi, sol.y(9,:)+sol.y(3,:)), ...
      trapz(xi, sum(sol.y([11 13 15 17],:),1)), trapz(xi, sol.y(19,:)+sol.y(13,:))];
D.native_trapz = nl;
D.G6c_residual = max(abs(nl - D.targets));
D.G6c_note = 'native adaptive-node trapezoid: DIAGNOSTIC ONLY (Packet-2 review 4)';

% ---- positivity ---------------------------------------------------------
conc = sol.y([1 3 5 7 9 11 13 15 17 19], :);
D.minimum = min(conc, [], 'all');
D.all_positive = D.minimum > 0;
end

function [m1, m2, resid, refine] = dense_masses(sol, tgt)
m1 = dense_at(sol, 20001);
m2 = dense_at(sol, 40001);
resid  = max(abs(m1 - tgt));
refine = max(abs(m1 - m2));
end

function m = dense_at(sol, nd)
xd = linspace(sol.x(1), sol.x(end), nd);
Y  = deval(sol, xd);
m = [trapz(xd, sum(Y([1 3 5 7],:),1)), trapz(xd, Y(9,:)+Y(3,:)), ...
     trapz(xd, sum(Y([11 13 15 17],:),1)), trapz(xd, Y(19,:)+Y(13,:))];
end
