function varargout = sa_deval(varargin)
%SA_DEVAL  Exact, vectorized replica of the toolbox deval for one-region bvp4c
%   solutions (M6 packet v0.2 section 3.2; RC-6.2, RC-6.3).
%
%   [Sx, Spx] = SA_DEVAL(sol, xint)
%   [Sx, Spx] = SA_DEVAL(sol, xint, idx)
%   Same outputs as R2025b deval.m (sha256 c8920dd29039...) bit for bit:
%   the same interval search (matlab.internal.math.discretize, as deval);
%   the stored values at mesh nodes; at t(end) the value y(:,end) and the
%   slope of the last interval at s = 1; inside an interval the ntrp3h
%   formula with the same per-element operation order:
%       h = t(k+1)-t(k); s = (tint-t(k))./h; s2 = s.*s; s3 = s.*s2;
%       slope = (y1-y0)./h; c = 3*slope - 2*yp0 - yp1; d = yp0 + yp1 - 2*slope;
%       yint = y0 + (h.*d.*s3 + h.*c.*s2 + h.*yp0.*s);
%       ypint = yp0 + (3*d.*s2 + 2*c.*s)
%   (deval's outer products (h*d)*s3 are rank one: one product per element,
%   so the elementwise form is identical per element). Nodes and t(end) use
%   the node/end abscissa itself, as deval does (s = +0 at a node, s = 1 at
%   t(end)). Unsorted input is sorted and the order restored; only the
%   requested rows are computed. The interpolation itself is unchanged, so
%   deval's cost (a search over the remaining points per interval) is gone.
%
%   The fast path returns double only and NEVER warns (sap_p28_solve reads
%   lastwarn after its diagnostics call deval). Everything else falls back
%   to the toolbox deval, which returns the same output or raises the same
%   error identifier (RC-6.3):
%     - a first argument that is not a scalar struct (incl. the swapped
%       (xint, sol) order);
%     - a struct without a solver field (the toolbox warns
%       MATLAB:deval:MissingSolverField and sets lastwarn) or a solver other
%       than 'bvp4c'; an idata field (idxNonNegative);
%     - sol.x that is not a finite, real, non-sparse double row of >= 2
%       strictly increasing points (repeated mesh points = interfaces;
%       decreasing x); sol.y / sol.yp not real double of the matching size;
%     - xint empty, not a real non-sparse double vector, or with a NaN,
%       Inf or out-of-range entry;
%     - idx not a nonempty real double vector of integers inside 1..n
%       (including 0: the toolbox raises MATLAB:badsubscript);
%     - more than two outputs or a wrong argument count.
%     - the deval mode of this MATLAB is not 'replica' (milestone 9, E3):
%       the replica is used only when the toolbox deval.m and
%       private/ntrp3h.m are the reference files it was proven against and
%       the self-check case (e) passed (sa_bvp4c_setup generation record,
%       cached at run start); otherwise every call goes to the toolbox deval.
%   sa_deval lives in solver/, which has no private/deval.m, so its own
%   deval is the toolbox's (pinned by T61 and the runtime evidence).
%   Reached through the one-line private forwarders deval.m in the package
%   root's private/ and in models/<id>/private/ (tests/ and the uncoupled
%   pipeline still reach the toolbox deval).

nout = max(1, nargout);
if ~fast_ok(nargin, nout, varargin)
    [varargout{1:nout}] = deval(varargin{:});
    return
end
sol = varargin{1};
t = sol.x; y = sol.y; yp = sol.yp;
n = size(y, 1);
N = numel(t);
if nargin < 3, idx = (1:n).'; else, idx = varargin{3}; idx = idx(:); end
tint = varargin{2};
tint = tint(:).';
m = numel(tint);
[tint, ord] = sort(tint);
k = matlab.internal.math.discretize(full(tint), t, false);
atend = tint == t(N);
k(atend) = N - 1;
isnode = ~atend & (tint == t(k));
te = tint;
te(isnode) = t(k(isnode));        % deval evaluates a node at t(bottom) itself
te(atend) = t(N);                 % and t(end) as t(end) over the last interval
tl = t(k);
h = t(k + 1) - tl;
s = (te - tl) ./ h;
s2 = s .* s;
s3 = s .* s2;
y0 = y(idx, k); y1 = y(idx, k + 1);
m0 = yp(idx, k); m1 = yp(idx, k + 1);
slope = (y1 - y0) ./ h;
c = 3*slope - 2*m0 - m1;
d = m0 + m1 - 2*slope;
Sx = y0 + (h.*d.*s3 + h.*c.*s2 + h.*m0.*s);
Sx(:, isnode) = y0(:, isnode);
if any(atend), Sx(:, atend) = repmat(y(idx, N), 1, nnz(atend)); end
Sx(:, ord) = Sx;
Sx = reshape(Sx, numel(idx), m);
varargout{1} = Sx;
if nout > 1
    Spx = m0 + (3*d.*s2 + 2*c.*s);
    Spx(:, ord) = Spx;
    varargout{2} = reshape(Spx, numel(idx), m);
end
end

function tf = fast_ok(nin, nout, a)
tf = false;
if nin < 2 || nin > 3 || nout > 2, return; end
if ~strcmp(sa_bvp4c_setup('deval_mode'), 'replica'), return; end
sol = a{1};
if ~(isstruct(sol) && isscalar(sol)), return; end
if ~all(isfield(sol, {'solver', 'x', 'y', 'yp'})) || isfield(sol, 'idata'), return; end
if ~(ischar(sol.solver) && strcmp(sol.solver, 'bvp4c')), return; end
x = sol.x; y = sol.y; yp = sol.yp;
if ~(isa(x, 'double') && isreal(x) && ~issparse(x) && isrow(x) && numel(x) >= 2), return; end
if ~(all(isfinite(x)) && all(diff(x) > 0)), return; end
if ~(isa(y, 'double') && isreal(y) && ~issparse(y) && ismatrix(y) && size(y, 1) >= 1 && ...
        size(y, 2) == numel(x)), return; end
if ~(isa(yp, 'double') && isreal(yp) && ~issparse(yp) && isequal(size(yp), size(y))), return; end
xi = a{2};
if ~(isa(xi, 'double') && isreal(xi) && ~issparse(xi) && isvector(xi) && ~isempty(xi)), return; end
if ~(all(isfinite(xi)) && min(xi) >= x(1) && max(xi) <= x(end)), return; end
if nin == 3
    idx = a{3};
    if ~(isa(idx, 'double') && isreal(idx) && ~issparse(idx) && isvector(idx) && ~isempty(idx)), return; end
    if ~(all(isfinite(idx)) && all(idx == fix(idx)) && all(idx >= 1) && all(idx <= size(y, 1))), return; end
end
tf = true;
end
