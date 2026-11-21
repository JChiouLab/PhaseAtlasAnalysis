function varargout = sa_identity(action, varargin)
%SA_IDENTITY  Root-identity helpers around the copied root-store kernel.
%
%   K = SA_IDENTITY('store_settings', M, theta, S, field_scale)
%       Kernel settings: comparison_mesh from the module, field_scale column,
%       similarity_tolerance and metric_id from the settings.
%   c = SA_IDENTITY('candidate', M, sol, target, attempt_id, observed, V, S, route, K)
%       Kernel candidate: mesh = the n_dense uniform grid over the
%       returned domain, profile = deval of the identity rows on it; the raw
%       mesh and profile stay in the attempt record. coordinate = the integer
%       lattice index as doubles; target_key = the lattice key.
%       route (default 'cold'): the discovery's route.
%       identity.candidate_mesh = 'comparison' evaluates the profile by
%       deval directly on the store's comparison mesh K.comparison_mesh (the
%       kernel's interpolation then is the identity); 'dense' (default) keeps
%       the n_dense grid.
%   P = SA_IDENTITY('pair', M, solA, solB, K, S)
%       Report distances for one pair: scaled_rms_v1 (kernel formula
%       on the comparison mesh; 'comparison': deval directly on that
%       mesh, 'dense': the kernel's linear interpolation of the dense grid),
%       relative-sup over the identity fields and
%       over all states (per component
%       max|a-b| / max(max|a|, max|b|), maximum over components, on a common
%       uniform grid of identity.relsup_nodes points), and their ratio.
%   a = SA_IDENTITY('amplification', profile, field_scale)
%       Per-field max|value| / field_scale.
%   [mesh, profile] = SA_IDENTITY('returned', M, sol, S, K)
%       The candidate's mesh and profile alone (the same expressions
%       'candidate' uses, so the bits are the same): the refill of a light
%       in-memory shard (sa_store_io 'full') recomputes a discovery's
%       returned_mesh / returned_profile from its attempt's sol with it.

switch action
    case 'store_settings', varargout{1} = store_settings(varargin{:});
    case 'candidate',      varargout{1} = candidate(varargin{:});
    case 'pair',           varargout{1} = pair(varargin{:});
    case 'amplification',  varargout{1} = amplification(varargin{:});
    case 'profile_distance', varargout{1} = profile_distance(varargin{:});
    case 'returned',       [varargout{1}, varargout{2}] = returned(varargin{:});
    otherwise, error('SA:IDENTITY', 'Unknown sa_identity action %s.', action);
end
end

function K = store_settings(M, theta, S, fs)
K = struct('comparison_mesh', M.comparison_mesh(theta, S), ...
    'field_scale', fs.values(:), ...
    'similarity_tolerance', S.identity.similarity_tolerance, ...
    'metric_id', S.identity.metric_id);
end

function c = candidate(M, sol, target, attempt_id, observed, V, S, route, K)
if nargin < 8 || isempty(route), route = 'cold'; end
if nargin < 9, K = []; end
[mesh, profile] = returned(M, sol, S, K);
c = struct('target_key', target.key, 'coordinate', double(target.index), ...
    'mesh', mesh, 'profile', profile, ...
    'seed_family', target.family, 'observed_family', observed, 'route', route, ...
    'attempt_id', attempt_id, 'valid', logical(V.valid), 'validity', V);
end

function [mesh, profile] = returned(M, sol, S, K)
% The kernel candidate's mesh and profile ('candidate' and the light-shard
% refill both use this one function)
if strcmp(S.identity.candidate_mesh, 'comparison')
    mesh = K.comparison_mesh;
else
    mesh = linspace(sol.x(1), sol.x(end), S.validity.n_dense);
end
profile = deval(sol, mesh, [M.identity_fields.row]);
end

function P = pair(M, solA, solB, K, S)
rows = [M.identity_fields.row];
if strcmp(S.identity.candidate_mesh, 'comparison')
    pa = deval(solA, K.comparison_mesh, rows);
    pb = deval(solB, K.comparison_mesh, rows);
else
    dense = @(s) deval(s, linspace(s.x(1), s.x(end), S.validity.n_dense), rows);
    pa = interp_to(solA, dense(solA), K.comparison_mesh, S);
    pb = interp_to(solB, dense(solB), K.comparison_mesh, S);
end
d = (pa - pb) ./ K.field_scale;
P = struct();
P.scaled_rms = sqrt(sum(d(:).^2) / numel(d));
P.relsup_identity = relsup(solA, solB, rows, S.identity.relsup_nodes);
P.relsup_all = relsup(solA, solB, 1:size(solA.y, 1), S.identity.relsup_nodes);
P.ratio = P.scaled_rms / P.relsup_identity;
end

function p = interp_to(sol, Y, cm, S)
mesh = linspace(sol.x(1), sol.x(end), S.validity.n_dense);
p = zeros(size(Y, 1), numel(cm));
for k = 1:size(Y, 1)
    p(k, :) = interp1(mesh, Y(k, :), cm, 'linear');   % the kernel's interpolation
end
end

function d = relsup(solA, solB, rows, n)
if abs(solA.x(1) - solB.x(1)) > 1e-12 || abs(solA.x(end) - solB.x(end)) > 1e-12
    d = NaN; return
end
xs = linspace(solA.x(1), solA.x(end), n);
A = deval(solA, xs, rows);
B = deval(solB, xs, rows);
d = 0;
for c = 1:numel(rows)
    s = max(max(abs(A(c, :))), max(abs(B(c, :))));
    if s > 0
        d = max(d, max(abs(A(c, :) - B(c, :))) / s);
    end
end
end

function d = profile_distance(M, solA, solB, K)
% scaled_rms_v1 between two iterates on the comparison mesh (deval; used for
% the report-only edge_large_change flag; iterates may sit at
% different coordinates, so this is never an identity distance).
rows = [M.identity_fields.row];
dd = (deval(solA, K.comparison_mesh, rows) - deval(solB, K.comparison_mesh, rows)) ./ K.field_scale;
d = sqrt(sum(dd(:).^2) / numel(dd));
end

function a = amplification(profile, fs)
a = max(abs(profile), [], 2) ./ fs(:);
end
