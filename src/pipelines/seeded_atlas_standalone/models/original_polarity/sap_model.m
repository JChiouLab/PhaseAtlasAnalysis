function varargout = sap_model(action, varargin)
%SAP_MODEL  Reaction terms and diffusion vector of the S. pombe model (Eq. 5).
%
%   Only mu and reaction are provided; parameter normalization lives in
%   sap_normalize.m, which derives theta from the user's parameters.
%
%   The two local functions below (MU_VEC and REACTION) are a verbatim copy
%   of the original model definition.
%
%   Actions
%     mu(th)          -> 5x1 diffusion vector at th.L
%     reaction(Y,th)  -> 5 x nx unclipped reaction

switch action
    case 'mu',       varargout{1} = mu_vec(varargin{:});
    case 'reaction', varargout{1} = reaction(varargin{:});
    otherwise
        error('sap_model:action', ...
            ['Unknown action %s. This file provides ' ...
             'mu and reaction only; parameters come from ' ...
             'the pipeline settings.'], action);
end
end

% ---- verbatim copy of the original model definition -----------------------
function mu = mu_vec(th)
m = th.Dm/(th.f*th.L^2); c = th.Dc/(th.f*th.L^2);
mu = [m; m; c; c; c];
end

function R = reaction(Y, th)
% UNCLIPPED. No max(0,.) anywhere.
u = Y(1,:); q = Y(2,:); v = Y(3,:); w = Y(4,:); p = Y(5,:);
a = th.alpha; b = th.beta; g = th.gamma; d = th.delta; e = th.epsilon;
Ru = (a*u.^2 + 1).*v - b*u - g*u.*p;
Rq =  g*u.*p - d*q;
Rv = -(a*u.^2 + 1).*v + b*u + e*w;
Rw =  d*q - e*w;
Rp = -g*u.*p + d*q;
R = [Ru; Rq; Rv; Rw; Rp];
end
% ---- end verbatim block -------------------------------------------------
