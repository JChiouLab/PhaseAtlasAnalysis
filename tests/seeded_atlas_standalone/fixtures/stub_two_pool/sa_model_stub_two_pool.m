function M = sa_model_stub_two_pool()
%SA_MODEL_STUB_TWO_POOL  Test-only second model (RC-11). Never solves.
%
%   Differs from the original model wherever shared code could hard-code it:
%   4 fields (not 10 concentrations), 8 states (not 28), 2 pools named R and B,
%   half domain [0, 1.5] (not [0, 0.5]), library types X and Y (not
%   HSS/HPSS/MPSS), shifts held as bvp4c `parameters` (not states 21:24),
%   identity on 4 rows, field scale of ones. Its "library" lives in memory.
%   State: [c1 c1' c2 c2' c3 c3' c4 c4']; ODE c_j'' = -k^2 (c_j - s_j),
%   k = pi / Lh, s = sol.parameters.

M = struct();
M.id = 'stub_two_pool';
M.version = 'stub_two_pool/1 (test fixture)';
M.unclipped = true;
M.n_states = 8;
M.fields = struct('name', {'c1','c2','c3','c4'}, 'row', {1, 3, 5, 7});
M.identity_fields = M.fields;
M.pools = struct('name', {'R', 'B'});
M.library_types = {'X', 'Y'};
M.domain = @(theta) [0 theta.Lh];
M.normalize = @normalize;
M.pool_targets = @pool_targets;
M.library = @library;
M.pair_guess = @pair_guess;
M.ode = @ode;
M.solve = @solve;
M.convergence = @convergence;
M.diagnostics = @diagnostics;
M.shifts = @(sol) struct('names', {{'s1','s2','s3','s4'}}, 'values', sol.parameters(:).', ...
    'max', max(abs(sol.parameters)));
M.shape_label = @(sol, theta, mo) 'X_X';
M.field_scale = @(theta, lib, S) struct('values', ones(4, 1), 'source', 'stub: ones');
M.comparison_mesh = @(theta, S) linspace(0, theta.Lh, S.identity.n_comparison_nodes);
M.required_checks = {'stub_noflux'};
M.model_options_defaults = struct('note', 'stub');
M.validate_model_options = @(mo) mo;
M.copies = struct('file', {}, 'copy_sha256', {});
end

function th = normalize(p)
th = struct('R', p.R, 'B', p.B, 'Lh', p.Lh);
end

function t = pool_targets(th, d)
v = [(th.R + d(1)) / 2, (th.R - d(1)) / 2, (th.B + d(2)) / 2, (th.B - d(2)) / 2];
t = struct('values', v, 'labels', {{'R_L','R_R','B_L','B_R'}}, ...
    'solver_target', struct('n', v, 'd', d));
end

function out = library(action, varargin)
switch action
    case 'load'
        th = varargin{2};
        res = 11;
        mY = true(res); mY(1:3, :) = false;          % asymmetric: exercises the mirror
        out = struct('types', {{'X','Y'}}, 'res', res, ...
            'mask', struct('X', true(res), 'Y', mY), ...
            'files', struct('type', {}, 'path', {}, 'bytes', {}, 'sha256', {}), ...
            'ax_row', linspace(-0.9 * th.B, 0.9 * th.B, res), ...
            'ax_col', linspace(-0.9 * th.R, 0.9 * th.R, res));
    case 'profiles'
        out = varargin{1};
    case 'cell'
        [lib, ~, d] = varargin{:};
        out = struct('in_range', false, 'cell', [NaN NaN], 'offset', [NaN NaN], 'reason', '');
        if abs(d(1)) > lib.ax_col(end) * (1 + 1e-12) || abs(d(2)) > lib.ax_row(end) * (1 + 1e-12)
            out.reason = 'outside stub range'; return
        end
        [e1, col] = min(abs(lib.ax_col - d(1)));
        [e2, row] = min(abs(lib.ax_row - d(2)));
        out.in_range = true; out.cell = [row col];
        out.offset = [e2 / diff(lib.ax_row(1:2)), e1 / diff(lib.ax_col(1:2))];
end
end

function [g, why] = pair_guess(lib, th, tL, cL, tR, cR, mo) %#ok<INUSD>
x = linspace(0, th.Lh, 41);
y = ones(8, numel(x));
y(1, :) = 1 + strcmp(tL, 'Y');                     % tells the stub corrector the left type
y(2:2:8, :) = 0;
g = struct('x', x, 'y', y); why = '';
end

function F = ode(x, Y, sol, th) %#ok<INUSL>
k = pi / th.Lh;
F = zeros(size(Y));
for j = 1:4
    F(2*j-1, :) = Y(2*j, :);
    F(2*j, :) = -k^2 * (Y(2*j-1, :) - sol.parameters(j));
end
end

function R = solve(theta, guess, target, solver, mo, corrector) %#ok<INUSD>
opts = struct('reltol', solver.reltol, 'abstol', solver.abstol, 'nmax', solver.nmax, ...
    'mesh_reject', solver.mesh_reject);
% M6: the mesh trace as the real adapters do (reset before the corrector,
% take after it; a stub corrector that notes nothing gives NaN)
tok = sa_mesh_trace('reset');
t0 = tic;
try
    [sol, info] = corrector(theta, guess, target, opts);
catch ME
    sa_mesh_trace('reset');
    rethrow(ME);
end
w = toc(t0);
tr = sa_mesh_trace('take', tok, opts.nmax);
reason = '';   % M6 tests only: an explicit reason_m6 passes through (no other stub sets it)
if isfield(info, 'reason_m6'), reason = info.reason_m6; end
R = struct('sol', sol, 'opts', opts, 'warning_id', info.warning_id, 'warning_msg', '', ...
    'corrector_class', info.class, 'corrector_reason', reason, 'wall_s', w, ...
    'corrector_nodal_min', NaN, 'mesh_trace', tr, 'max_proposed_mesh', tr.max_proposed);
end

function C = convergence(sol, opts, wid, wmsg) %#ok<INUSD>
C = struct('class', 'STUB', 'stats_maxres', NaN, 'stats_nmeshpoints', NaN, ...
    'warn_reltol_not_met', strcmp(wid, 'MATLAB:bvp4c:RelTolNotMet'), ...
    'warn_illcond_jac', strcmp(wid, 'MATLAB:bvp4c:IllCondJac'));
if isfield(sol, 'stats')
    C.stats_maxres = sol.stats.maxres; C.stats_nmeshpoints = sol.stats.nmeshpoints;
end
end

function c = diagnostics(sol, theta, target, mo) %#ok<INUSD>
c = struct('name', 'stub_noflux', 'value', max(abs(sol.y(2:2:8, 1))), ...
    'tolerance', 1e-8, 'comparison', '<', 'gating', true, 'owner', 'module');
end
