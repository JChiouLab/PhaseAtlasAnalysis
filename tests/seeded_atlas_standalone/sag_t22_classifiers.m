function rows = sag_t22_classifiers(E) %#ok<INUSD>
%SAG_T22_CLASSIFIERS  The pre-seeding classifier and the myxo
%   comparator on planted pairs, including the row mapping and a planted
%   1-ulp change.
T = 'T22 classifiers';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
A = 1 + rand(RandStream('mt19937ar', 'Seed', 22), 4, 121);
B = A; B(3, 50) = B(3, 50) * (1 + eps);
d = sag_preseed_classify('distance', A, B);
base = struct('returned', true, 'converged', true, 'finite', true, 'valid', true, 'type_only', false, 'd', d);
c = sag_preseed_classify('class', base);
rows = sa_check(rows, T, 'section 5: a 1-ulp change is detected (d > 0) and classed confirmed', ...
    d > 0 && d < 1e-15 && strcmp(c, 'confirmed'), sprintf('d %.3g', d));
cases = {setf(base, 'd', 1e-6), 'confirmed'; setf(base, 'd', 2e-6), 'band'; ...
    setf(base, 'd', 1e-2), 'moved'; setf(setf(base, 'valid', false), 'd', 5e-2), 'moved'; ...
    setf(setf(base, 'type_only', true), 'valid', false), 'moved'; ...
    setf(setf(base, 'valid', false), 'd', 1e-7), 'not_confirmed'; ...
    setf(setf(base, 'converged', false), 'd', NaN), 'not_confirmed'; ...
    setf(setf(base, 'returned', false), 'valid', false), 'not_confirmed'};
got = cellfun(@(r) sag_preseed_classify('class', r), cases(:, 1), 'UniformOutput', false);
rows = sa_check(rows, T, 'section 5 classes: edges 1e-6 (confirmed) and 1e-2 (moved), type-only moved, failures not_confirmed', ...
    isequal(got, cases(:, 2)), strjoin(got', ','));
% row mapping of the myxo comparator
x = linspace(0, 3, 11);
Y = repmat((1:24).', 1, numel(x));
sol = struct('solver', 'bvp4c', 'x', x, 'y', Y, 'yp', zeros(24, numel(x)));
P = sag_myxo_compare('ours_on_myxo', sol, linspace(0, 3, 241));
rows = sa_check(rows, T, 'row mapping: reference prof 1..8 <- ours [1 3 5 7 9 11 13 15]; shift <- rows 17:20', ...
    isequal(P.prof(:, 1).', [1 3 5 7 9 11 13 15]) && isequal(P.shift, 17:20), mat2str(P.prof(:, 1).'));
M8 = 1 + rand(RandStream('mt19937ar', 'Seed', 23), 8, 241);
O8 = M8; O8(7, 100) = O8(7, 100) * (1 + eps);
d0 = sag_myxo_compare('distance', O8, M8);
cl = {sag_myxo_compare('class', d0), sag_myxo_compare('class', 1e-6), sag_myxo_compare('class', 1.1e-6), ...
    sag_myxo_compare('class', 1e-2), sag_myxo_compare('class', NaN)};
rows = sa_check(rows, T, 'section 6c classes: 1-ulp agree (d > 0); 1e-6 agree; band; 1e-2 different; NaN not_computed', ...
    d0 > 0 && isequal(cl, {'agree', 'agree', 'band', 'different', 'not_computed'}), sprintf('d %.3g', d0));
bd = {sag_myxo_compare('band_decision', 1e-7, 5e-3), sag_myxo_compare('band_decision', 1e-5, 1e-7), ...
    sag_myxo_compare('band_decision', 1e-5, 1e-5), sag_myxo_compare('band_decision', 1e-7, 2e-2)};
rows = sa_check(rows, T, 'band decision: reference accuracy / our accuracy / different / different', ...
    isequal(bd, {'same_solution_myxo_accuracy', 'same_solution_our_accuracy', 'different', 'different'}), ...
    strjoin(bd, ','));
ij = {sag_myxo_compare('cell', [0 0]), sag_myxo_compare('cell', [-12 5]), sag_myxo_compare('cell', [13 0])};
rows = sa_check(rows, T, 'lattice index [kR kB] -> reference-atlas cell (13 + kB, 13 + kR); |k| = 13 outside', ...
    isequal(ij{1}, [13 13]) && isequal(ij{2}, [18 1]) && all(isnan(ij{3})), '');
rf = [sag_myxo_compare('reference', M8, true), sag_myxo_compare('reference', -M8, true), ...
    sag_myxo_compare('reference', M8, false)];
rows = sa_check(rows, T, 'reference cell = ok_grid and nodal min > 0', isequal(rf, [true false false]), '');
end

function s = setf(s, f, v)
s.(f) = v;
end
