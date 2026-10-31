function theta = sap_normalize(p)
%SAP_NORMALIZE  Dimensional parameters -> the corrector's normalized theta.
%
%   theta = SAP_NORMALIZE(p) for the original five-field model. Required
%   fields of p (dimensional, as in the uncoupled library's `par`):
%     a b c d e f    reaction rates
%     Dm Dc          membrane / cytosolic diffusion
%     L              FULL domain length (the half domain is L/2)
%     M_BG M_P       total GEF-pool and PAK-pool masses
%   Unknown fields are rejected; every value must be finite and positive.
%
%   Transformations (concentration unit C0 = M_BG/L, time unit
%   1/f, unit domain eta = x/L):
%     alpha = a (M_BG/L)^2 / f     gamma = c (M_BG/L) / f
%     beta = b/f  delta = d/f  epsilon = e/f
%     rho = M_P / M_BG             (masses are (1, rho))
%   The corrector forms mu = D/(f L^2) itself from Dm, Dc, f and L, so those
%   are carried dimensionally. C0 is carried for the library adapter.

need = {'a','b','c','d','e','f','Dm','Dc','L','M_BG','M_P'};
id = 'SAP:NORMALIZE';
if ~(isstruct(p) && isscalar(p))
    error(id, 'parameters must be a scalar struct.');
end
have = fieldnames(p);
missing = setdiff(need, have);
if ~isempty(missing)
    error(id, 'parameters are missing: %s.', strjoin(missing, ', '));
end
extra = setdiff(have, need);
if ~isempty(extra)
    error(id, ['parameters carry unknown field(s) %s; the original model takes ' ...
        'exactly %s.'], strjoin(extra', ', '), strjoin(need, ', '));
end
for k = 1:numel(need)
    v = p.(need{k});
    if ~(isa(v,'double') && isreal(v) && isscalar(v) && isfinite(v) && v > 0)
        error(id, 'parameter %s must be a finite positive double scalar.', need{k});
    end
end
if p.M_P >= p.M_BG
    error(id, 'M_P (%g) must be smaller than M_BG (%g): rho = M_P/M_BG < 1.', p.M_P, p.M_BG);
end

C0 = p.M_BG / p.L;
theta = struct();
theta.alpha   = p.a * C0^2 / p.f;
theta.beta    = p.b / p.f;
theta.gamma   = p.c * C0 / p.f;
theta.delta   = p.d / p.f;
theta.epsilon = p.e / p.f;
theta.rho     = p.M_P / p.M_BG;
theta.f       = p.f;
theta.Dm      = p.Dm;
theta.Dc      = p.Dc;
theta.L       = p.L;
theta.C0      = C0;
theta.M_BG    = p.M_BG;
theta.convention = ['unit domain eta = x/L; concentrations / C0 (C0 = M_BG/L); ' ...
    'time * f; masses (1, rho); mu = D/(f L^2)'];
end
