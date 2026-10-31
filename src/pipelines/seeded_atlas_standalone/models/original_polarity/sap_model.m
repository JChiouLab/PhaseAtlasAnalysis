function varargout = sap_model(action, varargin)
%SAP_MODEL  PARTIAL copy of the campaign cha_model.m: mu and reaction only.
%
%   DECLARED PARTIAL COPY (milestone-2 packet section 2, "Copy only mu/reaction
%   from cha_model.m"). The source file's frozen theta is the h04 L = 7.90 point
%   and must never be used by this pipeline, so the 'theta' action and every
%   other action of the source are deliberately absent. Parameter normalization
%   lives in sap_normalize.m, which derives theta from the user's parameters.
%
%   The two local functions below (MU_VEC and REACTION) are BYTE-IDENTICAL to
%   lines 51-66 of the source; sa_test_copies verifies that line range verbatim.
%   Source: data/experiments/coupled_half_domain_atlas_20260811/matlab/cha_model.m
%
%   Actions
%     mu(th)          -> 5x1 diffusion vector at th.L
%     reaction(Y,th)  -> 5 x nx unclipped reaction

switch action
    case 'mu',       varargout{1} = mu_vec(varargin{:});
    case 'reaction', varargout{1} = reaction(varargin{:});
    otherwise
        error('sap_model:action', ...
            ['Unknown or deliberately omitted action %s. This file is a PARTIAL ' ...
             'copy carrying mu and reaction only; the source''s frozen theta is ' ...
             'never used by the standalone pipeline.'], action);
end
end

% ---- VERBATIM from cha_model.m lines 51-66 ------------------------------
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
