function M = sa_model_stub_asym()
%SA_MODEL_STUB_ASYM  The stub model with asymmetric pool targets (M7a T68):
%   the R pool's right half target shrinks twice as fast, so the admissible
%   half-widths of axis 1 differ in the two directions.
M = sa_model_stub_two_pool();
M.id = 'stub_asym';
M.version = 'stub_asym/1 (test fixture)';
M.pool_targets = @pool_targets;
end

function t = pool_targets(th, d)
v = [(th.R + d(1)) / 2, (th.R - 2 * d(1)) / 2, (th.B + d(2)) / 2, (th.B - d(2)) / 2];
t = struct('values', v, 'labels', {{'R_L','R_R','B_L','B_R'}}, ...
    'solver_target', struct('n', v, 'd', d));
end
