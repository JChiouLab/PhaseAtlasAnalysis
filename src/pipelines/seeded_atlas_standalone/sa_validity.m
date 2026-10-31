function varargout = sa_validity(action, varargin)
%SA_VALIDITY  The acceptance boundary.
%
%   V = SA_VALIDITY('evaluate', M, R, theta, target, S)
%       Evaluate every check on one returned iterate R (from M.solve: R.sol,
%       R.opts, R.warning_id, R.warning_msg) and return the validity record.
%   [valid, failing] = SA_VALIDITY('verdict', V, M)
%       The verdict rule applied to a freshly evaluated record. Fail closed:
%       a check passes only if its recorded value and tolerance are
%       finite and the declared comparison holds; a missing required check, a
%       duplicated one, or a non-gating copy of a required check is invalid.
%       A STORED verdict is never read on reload, resume or export; use
%       'rederive'.
%   [V, info] = SA_VALIDITY('rederive', M, rec, theta, S, opts_expected, prov)
%       Verdict re-derivation. rec is a stored attempt record (sol, opts,
%       warning_id, warning_msg, validity, target, optional provenance).
%       (i) opts_expected (rebuilt from the effective settings and the
%       attempt's kind) must equal rec.opts (isequal); (ii) the boundary is
%       re-evaluated on the stored iterate and warning; (iii) when rec was
%       written under the same software_version, matlab_version and thread
%       mode as prov, every recomputed check (name, value, tolerance,
%       comparison, gating) must equal the stored one (isequaln) and so must
%       the check list; (iv) the recomputed verdict is returned. An inequality
%       in (i) or (iii) sets info.integrity_ok = false with the reason (the
%       caller stops). Older records are re-derived under the current rules
%       and info reports the stored verdict BESIDE the recomputed one (no
%       integrity failure, no relabel). Stored stats and warning id stay
%       trusted inputs.
%   names = SA_VALIDITY('shared_required')
%   r = SA_VALIDITY('bvp4c_residual', M, sol, theta, opts)   recomputed residual
%
%   SHARED checks, all gating unless noted:
%     convergence_maxres         stats.maxres < RelTol (missing stats -> NaN)
%     convergence_proposed_mesh  stats.nmeshpoints <= NMax
%     convergence_warning        no RelTolNotMet / IllCondJac warning (0 == 0)
%     mesh                       returned mesh points < mesh_reject
%     finite                     number of entries of sol.y, sol.parameters
%                                and sol.yp (when present) that are non-finite
%                                OR have a nonzero imaginary part == 0; every
%                                imaginary entry and every non-finite state or
%                                parameter enters the guard that skips the
%                                positivity, residual and module checks, so a
%                                complex iterate fails closed before min
%     positivity                 signed minimum, over every returned mesh
%                                interval and every module concentration row,
%                                of the EXACT minimum of bvp4c's cubic Hermite
%                                interpolant (node values y, slopes yp,
%                                interior critical points; sa_hermite_min, a
%                                verbatim copy of gnx_hermite_min, source
%                                SHA-256 10b230f7d82fcadbe81f8f369132ead5b80fa1614d94209915cd3804638ae39a);
%                                must be > 0 and finite; non-finite y or yp in
%                                a row gives NaN (refused)
%     positivity_sampled (record only, not gating)  the earlier sampled
%                                minimum (mesh nodes, interval fractions,
%                                n_dense grid); non-finite sample -> NaN
%     residual (NOT gating)       per row i, dense point k:
%                                r_ik = |y'_ik - F_ik| / max(|F_ik|, AbsTol/RelTol),
%                                y' from deval, F from M.ode; value = sup over
%                                all state rows; flagged when > factor*RelTol
%     residual_bvp4c             bvp4c's interval residual norm recomputed
%                                from the stored x, y, yp and M.ode with
%                                bvp4c's scaling (R2025b bvp4c.m:324-418:
%                                Lobatto IIIa midpoint Newton residual and the
%                                two Lobatto points, max(|f|, AbsTol/RelTol),
%                                sqrt(h/2 * sum)); value = max over intervals,
%                                must be < opts.reltol. Gating since the
%                                agreement test passed (residual_bvp4c_gating).
%   MODULE checks come from M.diagnostics; their names must include
%   M.required_checks.

switch action
    case 'evaluate',        varargout{1} = evaluate(varargin{:});
    case 'verdict',         [varargout{1:max(nargout,1)}] = verdict(varargin{:});
    case 'shared_required', varargout{1} = shared_required();
    case 'rederive',        [varargout{1}, varargout{2}] = rederive(varargin{:});
    case 'bvp4c_residual',  varargout{1} = bvp4c_residual(varargin{:});
    case 'residual_bvp4c_gating', varargout{1} = residual_bvp4c_gating();
    otherwise, error('SA:VALIDITY', 'Unknown sa_validity action %s.', action);
end
end

% =========================================================================

function tf = residual_bvp4c_gating()
% Adoption switch of the recomputed bvp4c residual as a gate: adopted after
% |recomputed - stats.maxres| = 0 held on every accepted record tested; now
% a gate and a required shared check.
tf = true;
end

function names = shared_required()
names = {'convergence_maxres','convergence_proposed_mesh','convergence_warning', ...
    'mesh','finite','positivity'};
if residual_bvp4c_gating(), names{end+1} = 'residual_bvp4c'; end
end

function V = evaluate(M, R, theta, target, S)
sol = R.sol;
opts = R.opts;
V = struct('schema', 'sa_validity/1', 'checks', [], 'valid', false, ...
    'failing', {{}}, 'flags', {{}}, 'positivity', struct(), 'residual', struct(), ...
    'convergence_class', 'NO_RETURN', 'module_error', '', 'returned', ~isempty(sol));
checks = rec('convergence_maxres', NaN, opts.reltol, '<', true);
checks(end+1) = rec('convergence_proposed_mesh', NaN, opts.nmax, '<=', true);
checks(end+1) = rec('convergence_warning', NaN, 0, '==', true);
checks(end+1) = rec('mesh', NaN, opts.mesh_reject, '<', true);
checks(end+1) = rec('finite', NaN, 0, '==', true);
checks(end+1) = rec('positivity', NaN, 0, '>', true);
checks(end+1) = rec('positivity_sampled', NaN, NaN, 'record', false);
checks(end+1) = rec('residual', NaN, S.validity.residual_flag_factor * opts.reltol, '<=', false);
checks(end+1) = rec('residual_bvp4c', NaN, opts.reltol, '<', residual_bvp4c_gating());

if ~isempty(sol)
    % ---- convergence (classifier copied with the corrector) -------------
    Cv = M.convergence(sol, opts, R.warning_id, R.warning_msg);
    V.convergence_class = Cv.class;
    checks = setv(checks, 'convergence_maxres', Cv.stats_maxres);
    checks = setv(checks, 'convergence_proposed_mesh', Cv.stats_nmeshpoints);
    if strcmp(R.warning_id, '<unknown>')
        checks = setv(checks, 'convergence_warning', NaN);
    else
        checks = setv(checks, 'convergence_warning', ...
            double(Cv.warn_reltol_not_met || Cv.warn_illcond_jac));
    end
    % ---- mesh and finiteness --------------------------------------------
    checks = setv(checks, 'mesh', numel(sol.x));
    % An entry is bad when non-finite or with a nonzero imaginary part
    badv = @(v) nnz(~isfinite(v(:)) | imag(v(:)) ~= 0);
    nbad_y = badv(sol.y);
    if isfield(sol, 'parameters'), nbad_y = nbad_y + badv(sol.parameters); end
    nbad = nbad_y;
    if isfield(sol, 'yp')
        nbad = nbad + badv(sol.yp);                                % deval's slopes
        nbad_y = nbad_y + nnz(imag(sol.yp(:)) ~= 0);               % guard: complex slopes
    end
    checks = setv(checks, 'finite', nbad);
    if nbad_y == 0
        % positivity is evaluated whenever the states are finite and real, so
        % a non-finite slope is also caught by its own exact-minimum rule
        [checks, V.positivity] = positivity(checks, M, sol, S);
    end
    if nbad == 0
        [checks, V.residual] = residual(checks, M, sol, theta, opts, S);
        if isfield(sol, 'yp')
            checks = setv(checks, 'residual_bvp4c', bvp4c_residual(M, sol, theta, opts));
        end
        % ---- module checks ----------------------------------------------
        try
            mc = M.diagnostics(sol, theta, target, S.model_options);
            for k = 1:numel(mc)
                checks(end+1) = rec(mc(k).name, mc(k).value, mc(k).tolerance, ...
                    mc(k).comparison, mc(k).gating, 'module'); %#ok<AGROW>
            end
        catch ME
            V.module_error = ME.message;
        end
    end
end
for k = 1:numel(checks)
    checks(k).pass = check_passes(checks(k));
end
V.checks = checks;
r = checks(strcmp({checks.name}, 'residual'));
if ~r.pass
    V.flags{end+1} = 'residual_scaled_sup_above_flag';
end
[V.valid, V.failing] = verdict(V, M);
end

function [checks, P] = positivity(checks, M, sol, S)
rows = [M.fields.row];
x = sol.x(:).';
P = struct();
% ---- exact minimum of the cubic Hermite interpolant, per row --------------
P.method = 'exact_hermite (sa_hermite_min)';
ex = nan(1, numel(rows)); xat = nan(1, numel(rows)); kat = nan(1, numel(rows));
if isfield(sol, 'yp')
    for r = 1:numel(rows)
        [ex(r), xat(r), kat(r)] = sa_hermite_min(x, sol.y(rows(r), :), sol.yp(rows(r), :));
    end
end
if any(isnan(ex))
    r = find(isnan(ex), 1); m = NaN;
else
    [m, r] = min(ex);
end
P.min = m;
P.row = M.fields(r).name;
P.xi = xat(r);
P.interval = kat(r);
P.row_minima = ex;
% ---- the sampled minimum, kept as a record (positivity_sampled) -----------
h = diff(x);
fr = S.validity.interval_samples(:);
inner = x(1:end-1) + fr .* h;                 % every interval, every fraction
xd = linspace(x(1), x(end), S.validity.n_dense);
pts = unique([x, inner(:).', xd]);
Y = deval(sol, pts, rows);
% MATLAB min skips NaN: count non-finite samples explicitly and fail closed
bad = ~isfinite(Y);
P.n_points = numel(pts);
P.n_nonfinite = nnz(bad);
if P.n_nonfinite > 0
    [rs, c] = ind2sub(size(Y), find(bad, 1));    % first non-finite sample
    ms = NaN;
else
    [ms, where] = min(Y(:));
    [rs, c] = ind2sub(size(Y), where);
end
P.sampled_min = ms;
P.sampled_row = M.fields(rs).name;
P.sampled_xi = pts(c);
P.nodal_min = min(sol.y(rows, :), [], 'all');
Yd = deval(sol, xd, rows);
if any(~isfinite(Yd(:))), P.uniform_grid_min = NaN; else, P.uniform_grid_min = min(Yd(:)); end
checks = setv(checks, 'positivity', m);
checks = setv(checks, 'positivity_sampled', ms);
end

function [checks, Q] = residual(checks, M, sol, theta, opts, S)
xd = linspace(sol.x(1), sol.x(end), S.validity.n_dense);
[Y, Yp] = deval(sol, xd);
F = M.ode(xd, Y, sol, theta);
floor_ = opts.abstol / opts.reltol;
D = abs(Yp - F);
Rr = D ./ max(abs(F), floor_);
[sup, where] = max(Rr(:));
[r, c] = ind2sub(size(Rr), where);
Q = struct('scaled_sup', sup, 'scaled_rms', sqrt(mean(Rr(:).^2)), ...
    'unscaled_sup', max(D(:)), 'sup_row', r, 'sup_xi', xd(c), ...
    'floor', floor_, 'n_points', numel(xd));
checks = setv(checks, 'residual', sup);
end

function res = bvp4c_residual(M, sol, theta, opts)
% bvp4c's residual (R2025b bvp4c.m colloc_RHS + residual), one
% region, from the stored iterate. F at the nodes and the Lobatto values come
% from M.ode; the Hermite interpolant uses the stored slopes yp, as bvp4c does.
x = sol.x(:).'; y = sol.y; yp = sol.yp;
if numel(x) < 2, res = NaN; return; end
h = diff(x);
threshold = opts.abstol / opts.reltol;
F = M.ode(x, y, sol, theta);
yi = y(:, 1:end-1); yip1 = y(:, 2:end);
Fi = F(:, 1:end-1); Fip1 = F(:, 2:end);
xip05 = (x(1:end-1) + x(2:end)) / 2;
yip05 = (yi + yip1) / 2 - (Fip1 - Fi) .* (h ./ 8);
Fmid = M.ode(xip05, yip05, sol, theta);
Phi = yip1 - yi - (Fip1 + 4 * Fmid + Fi) .* (h / 6);       % Lobatto IIIa residual
lob4 = (1 + sqrt(3/7)) / 2; lob2 = (1 - sqrt(3/7)) / 2;
lobw24 = 49/90; lobw3 = 32/45;
temp = (Phi .* (1.5 ./ h)) ./ max(abs(Fmid), threshold);
r2 = lobw3 * dot(temp, temp, 1);
for lob = [lob2, lob4]
    [yL, ypL] = hermite_at(lob, h, y, yp);
    fL = M.ode(x(1:end-1) + lob * h, yL, sol, theta);
    temp = (ypL - fL) ./ max(abs(fL), threshold);
    r2 = r2 + lobw24 * dot(temp, temp, 1);
end
r = sqrt(abs(h / 2) .* r2);
if any(~isfinite(r)), res = NaN; else, res = max(r); end
end

function [Sx, Spx] = hermite_at(hnode, diffx, y, yp)
% bvp4c's interp_Hermite (same arithmetic order).
N = size(y, 2);
y_N1 = y(:, 1:N-1); y_2N = y(:, 2:N);
yp_N1 = yp(:, 1:N-1); yp_2N = yp(:, 2:N);
rowscal = 1 ./ diffx;
slope = (y_2N - y_N1) .* rowscal;
c = (3*slope - 2*yp_N1 - yp_2N) .* rowscal;
d = (yp_N1 + yp_2N - 2*slope) .* (rowscal.^2);
rowscal2 = hnode .* diffx;
Sx = ((d.*rowscal2 + c).*rowscal2 + yp_N1).*rowscal2 + y_N1;
Spx = (3.*d.*rowscal2 + 2.*c).*rowscal2 + yp_N1;
end

% =========================================================================

function [V, info] = rederive(M, recd, theta, S, opts_expected, prov)
info = struct('integrity_ok', true, 'reason', '', 'same_version', false, ...
    'valid_stored', NaN, 'valid_recomputed', NaN, 'flip', false, 'mismatches', {{}});
if isfield(recd, 'validity') && isstruct(recd.validity) && isfield(recd.validity, 'valid')
    info.valid_stored = double(recd.validity.valid);
end
% (i) options rebuilt from the effective settings and the attempt's kind
if ~isequal(opts_expected, recd.opts)
    info.integrity_ok = false;
    info.reason = 'opts differ from the options rebuilt from the effective settings';
end
% (ii) re-evaluate on the stored iterate and warning
R = struct('sol', recd.sol, 'opts', recd.opts, 'warning_id', recd.warning_id, ...
    'warning_msg', recd.warning_msg);
V = evaluate(M, R, theta, recd.target, S);
info.valid_recomputed = double(V.valid);
info.flip = ~isnan(info.valid_stored) && info.valid_stored ~= info.valid_recomputed;
% (iii) same software, MATLAB version and thread mode: every check identical
if isfield(recd, 'provenance') && isstruct(recd.provenance) && ...
        isequal(recd.provenance, prov)
    info.same_version = true;
    cs = recd.validity.checks; cn = V.checks;
    if numel(cs) ~= numel(cn) || ~isequal({cs.name}, {cn.name})
        info.mismatches{end+1} = 'check list';
    else
        for k = 1:numel(cn)
            for f = {'value', 'tolerance', 'comparison', 'gating'}
                if ~isequaln(cs(k).(f{1}), cn(k).(f{1}))
                    info.mismatches{end+1} = sprintf('%s.%s', cn(k).name, f{1}); %#ok<AGROW>
                end
            end
        end
    end
    if ~isempty(info.mismatches)
        info.integrity_ok = false;
        info.reason = strtrim([info.reason ' same-version check mismatch: ' ...
            strjoin(info.mismatches, ', ')]);
    end
end
end

% =========================================================================

function [valid, failing] = verdict(V, M)
need = [shared_required(), M.required_checks];
failing = {};
if ~isstruct(V) || ~isfield(V, 'checks') || ~isstruct(V.checks) || isempty(V.checks) || ...
        ~all(isfield(V.checks, {'name','value','tolerance','comparison','gating'}))
    valid = false; failing = {'<malformed validity record>'}; return
end
names = {V.checks.name};
for k = 1:numel(need)
    hit = find(strcmp(names, need{k}));
    if numel(hit) ~= 1
        failing{end+1} = sprintf('%s (missing or duplicated)', need{k}); %#ok<AGROW>
    elseif ~isequal(V.checks(hit).gating, true)
        failing{end+1} = sprintf('%s (not gating)', need{k}); %#ok<AGROW>
    end
end
for k = 1:numel(V.checks)
    c = V.checks(k);
    if isequal(c.gating, true) && ~check_passes(c)
        failing{end+1} = c.name; %#ok<AGROW>
    end
end
failing = unique(failing, 'stable');
valid = isempty(failing);
end

function tf = check_passes(c)
v = c.value; t = c.tolerance;
if ~(isnumeric(v) && isscalar(v) && isfinite(v))
    tf = false; return
end
switch c.comparison
    case 'record', tf = true; return      % recorded value, nothing to compare
end
if ~(isnumeric(t) && isscalar(t) && isfinite(t))
    tf = false; return
end
switch c.comparison
    case '<',  tf = v <  t;
    case '<=', tf = v <= t;
    case '>',  tf = v >  t;
    case '==', tf = v == t;
    otherwise, tf = false;
end
end

function r = rec(name, value, tol, cmp, gating, owner)
if nargin < 6, owner = 'shared'; end
r = struct('name', name, 'value', double(value), 'tolerance', double(tol), ...
    'comparison', cmp, 'gating', logical(gating), 'owner', owner, 'pass', false);
end

function checks = setv(checks, name, value)
k = strcmp({checks.name}, name);
checks(k).value = double(value);
end
