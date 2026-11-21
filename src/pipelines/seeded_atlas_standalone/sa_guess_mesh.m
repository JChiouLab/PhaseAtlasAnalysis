function varargout = sa_guess_mesh(action, varargin)
%SA_GUESS_MESH  Guess mesh for every route that passes a stored iterate as a
%   guess.
%
%   g = SA_GUESS_MESH('from', sol, S)
%       The guess a continuation edge, substep, restart level (option
%       'myxo_restart') or the diagonal/deep rescue rungs pass to the corrector
%       (the coupled bvp4c solve) from the stored iterate sol (all states,
%       the uniform sources s_L, s_R ('shifts') and accumulators included).
%   g = SA_GUESS_MESH('secant', P1, P2, S)
%       The secant guess 2*P1 - P2 evaluated on the guess mesh of P1.
%   x = SA_GUESS_MESH('mesh', sol, S)       the guess mesh alone
%   t = SA_GUESS_MESH('tag', S)              '' for 'parent', else
%                                            ' [<rule> <nodes>]' (guess_source)
%
%   Rule = S.continuation.guess_mesh, N = S.continuation.guess_nodes:
%     'parent'    (DEFAULT) the parent's own mesh:
%                 g = struct('x', sol.x, 'y', sol.y); secant
%                 2*deval(P1, P1.x) - deval(P2, P1.x).
%     'uniform'   N equally spaced nodes on [sol.x(1), sol.x(end)] (the
%                 shape of both models' cold guess meshes); y = deval(sol, x).
%     'subsample' the parent's mesh coarsened to at most N of its own nodes,
%                 density-preserving: n <= N keeps the parent mesh; n > N keeps
%                 the nodes round(linspace(1, n, N)) (exactly N distinct
%                 nodes, both ends included), y = sol.y at those nodes (no
%                 interpolation).
%   The guess carries x and y only (the form both models pass to bvp4c);
%   bvp4c recomputes yp. Deterministic: a pure function of the stored iterate
%   and the settings (no clock, no random number). Model-agnostic: only
%   sol.x, sol.y and deval are used. Multipoint iterates (repeated x) are
%   refused under 'uniform' and 'subsample' (fail closed; no current model
%   is multipoint). Cold calls and refinement do not use this function.
%   (i) Replay never rebuilds a guess (a replayed attempt is taken
%   from its file); the rule affects live calls only, so resume equivalence
%   is tested on stub and real kill-resume runs. (ii)
%   guess_nodes is inert under 'parent' but still affects results: a resume with it
%   changed is refused (SA:RESUME:E_SETTINGS), as for a changed rule.
%   (iii) 'uniform' ENLARGES a parent with n < N nodes to N nodes;
%   'subsample' never enlarges (n <= N passes the parent mesh unchanged).

switch action
    case 'from',   varargout{1} = from_sol(varargin{:});
    case 'secant', varargout{1} = secant(varargin{:});
    case 'mesh',   [varargout{1}, varargout{2}] = mesh_of(varargin{:});
    case 'tag',    varargout{1} = tag(varargin{:});
    otherwise, error('SA:GUESSMESH', 'Unknown sa_guess_mesh action %s.', action);
end
end

function [rule, N] = rule_of(S)
if nargin < 1 || isempty(S)
    rule = 'parent'; N = Inf; return
end
rule = S.continuation.guess_mesh;
N = S.continuation.guess_nodes;
end

function t = tag(S)
[rule, N] = rule_of(S);
if strcmp(rule, 'parent'), t = ''; else, t = sprintf(' [%s %d]', rule, N); end
end

function [x, idx] = mesh_of(sol, S)
% idx: node indices of sol.x (subsample / parent) or [] (uniform)
[rule, N] = rule_of(S);
n = numel(sol.x);
switch rule
    case 'parent'
        x = sol.x; idx = 1:n;
    case 'uniform'
        no_multipoint(sol);
        x = linspace(sol.x(1), sol.x(end), N); idx = [];
    case 'subsample'
        no_multipoint(sol);
        if n <= N
            idx = 1:n;
        else
            idx = round(linspace(1, n, N));
        end
        x = sol.x(idx);
    otherwise
        error('SA:GUESSMESH', 'Unknown continuation.guess_mesh %s.', rule);
end
end

function g = from_sol(sol, S)
if nargin < 2, S = []; end
[rule, ~] = rule_of(S);
if strcmp(rule, 'parent')
    g = struct('x', sol.x, 'y', sol.y);          % the parent's mesh, bit for bit
    return
end
[x, idx] = mesh_of(sol, S);
if isempty(idx)
    y = deval(sol, x);
else
    y = sol.y(:, idx);
end
g = struct('x', x, 'y', y);
end

function g = secant(P1, P2, S)
if nargin < 3, S = []; end
x1 = mesh_of(P1, S);
g = struct('x', x1, 'y', 2 * deval(P1, x1) - deval(P2, x1));   % secant form on the rule's mesh
end

function no_multipoint(sol)
if any(diff(sol.x) == 0)
    error('SA:GUESSMESH:E_MULTIPOINT', ['continuation.guess_mesh other than ''parent'' is ' ...
        'defined for two-point iterates only; this iterate repeats an x node (multipoint).']);
end
end
