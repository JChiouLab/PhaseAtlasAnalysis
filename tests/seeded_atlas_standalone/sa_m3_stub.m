function varargout = sa_m3_stub(action, varargin)
%SA_M3_STUB  Scripted stub corrector for the zero-solve tests of the search
%   (sa_t24 to sa_t40). Runs on the stub_two_pool fixture model (4 fields, 2 pools,
%   shifts as bvp4c parameters, half domain [0, 1.5]); never calls bvp4c.
%
%   p = SA_M3_STUB('params')                 R = B = 4, Lh = 1.5
%   s = SA_M3_STUB('settings', run_dir, box) lattice step 0.65 (13 x 13
%       nodes, |k| = 6 is the ring: outside the stub library's native range),
%       cold at the anchor only (X_X), refinement off, continuation + rescue
%       on, checkpoints every 5 calls, no disk floor; the explicit
%       rescue ladder (diagonal, secant, deep; ring true)
%   st = SA_M3_STUB('state', ...)            handle state; name/value pairs:
%       'fail'   fn(d, branch, exact) -> 'S' (valid) | 'R' (returned, maxres
%                1: refused) | 'N' (no return); default always 'S'
%       'seq'    char sequence consumed first, one letter per call (S/R/N)
%       'kill_at' call number at which the corrector throws SA:TEST:KILL
%                (a kill mid-call: journal start written, no attempt file)
%       'force'  fn(d, branch) -> branch actually returned ('' = keep)
%       'log'    {}: every guess the corrector receives is appended (sa_t42_am8)
%   c = SA_M3_STUB('corrector', st)          the injected corrector
%   n = SA_M3_STUB('calls', st)
%   b = SA_M3_STUB('branch_of', guess)       'A' or 'B'
%
%   Returned iterate at target n (half pools) and branch b:
%       c_j(x) = n_j (1 + a_b cos(k x)), k = pi/Lh, a_A = 0.1, a_B = 0.3,
%   an exact solution of the stub ODE c'' = -k^2 (c - s), s = n (bvp4c
%   parameters), with exact slopes (201 nodes). The branch follows the
%   guess (amplitude of row 1; a constant cold guess with y(1,1) = 2, the Y
%   left type, gives B), so roots propagate along continuation. 'exact' is
%   true when the guess's first values of all four fields equal the
%   target's own values (a linear secant extrapolation does; a neighbour's
%   iterate does not).

switch action
    case 'params', varargout{1} = struct('R', 4, 'B', 4, 'Lh', 1.5);
    case 'settings', varargout{1} = settings(varargin{:});
    case 'state', varargout{1} = make_state(varargin{:});
    case 'corrector'
        st = varargin{1};
        varargout{1} = @(th, g, t, o) stub_call(th, g, t, o, st);
    case 'calls', varargout{1} = varargin{1}('n');
    case 'branch_of', varargout{1} = branch_of(varargin{1});
    case 'sol', varargout{1} = make_sol(varargin{:});
    otherwise, error('SA:M3STUB', 'Unknown action %s.', action);
end
end

function s = settings(run_dir, box)
if nargin < 2, box = []; end
s = struct('model', 'stub_two_pool', 'run_dir', run_dir);
s.lattice = struct('step', 0.65);
s.families = {'X_X'};
s.targets = struct('policy', 'explicit', 'list', [0 0]);
s.validity = struct('refine_enabled', false);
s.routes = struct('cold', true, 'continuation', true, 'rescue', true);
s.continuation = struct('box', box);
% The explicit original rescue ladder, so the earlier stub tests keep their
% meaning under the mode-layer defaults (secant only, no
% ring rescue)
s.rescue = struct('rungs', {{'diagonal', 'secant', 'deep'}}, 'ring', true);
s.budget = struct('max_calls', 5000);
s.checkpoint = struct('every_calls', 5, 'min_free_gib', 0);
end

function st = make_state(varargin)
st = containers.Map('KeyType', 'char', 'ValueType', 'any');
st('n') = 0; st('seq') = ''; st('kill_at') = Inf;
st('fail') = @(d, b, ex) 'S';
st('force') = @(d, b) '';
for k = 1:2:numel(varargin)
    st(varargin{k}) = varargin{k+1};
end
end

function b = branch_of(g)
r = g.y(1, :);
amp = (max(r) - min(r)) / max(abs(max(r) + min(r)), eps);
if amp > 0.2
    b = 'B';
elseif amp < 1e-12
    if g.y(1, 1) == 2, b = 'B'; else, b = 'A'; end
else
    b = 'A';
end
end

function sol = make_sol(th, n, b)
k = pi / th.Lh;
x = linspace(0, th.Lh, 201);
if b == 'B', a = 0.3; else, a = 0.1; end
y = zeros(8, numel(x)); yp = y;
for j = 1:4
    y(2*j-1, :) = n(j) * (1 + a * cos(k * x));
    y(2*j, :) = -n(j) * a * k * sin(k * x);
    yp(2*j-1, :) = y(2*j, :);
    yp(2*j, :) = -n(j) * a * k^2 * cos(k * x);
end
sol = struct('solver', 'bvp4c', 'x', x, 'y', y, 'yp', yp, 'parameters', n(:), ...
    'stats', struct('maxres', 1e-9, 'nmeshpoints', numel(x)));
end

function [sol, info] = stub_call(th, g, t, o, st) %#ok<INUSD>
n = st('n') + 1; st('n') = n;
if isKey(st, 'log')                 % guess-mesh tests: every guess received, in call order
    st('log') = [st('log'), {g}];
end
if n == st('kill_at')
    error('SA:TEST:KILL', 'stub kill at call %d', n);
end
b = branch_of(g);
fb = feval(st('force'), t.d, b);
if ~isempty(fb), b = fb; end
if b == 'B', a = 0.3; else, a = 0.1; end
exact = all(abs(g.y([1 3 5 7], 1) - t.n(:) * (1 + a)) <= 1e-9 * t.n(:));
seq = st('seq');
if ~isempty(seq)
    c = seq(1); st('seq') = seq(2:end);
else
    c = feval(st('fail'), t.d, b, exact);
end
info = struct('class', 'STUB', 'warning_id', '');
switch c
    case 'S', sol = make_sol(th, t.n, b);
    case 'R', sol = make_sol(th, t.n, b); sol.stats.maxres = 1;
    case 'N', sol = [];
end
end
