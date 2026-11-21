function out = sag_myxo_compare(action, varargin)
%SAG_MYXO_COMPARE  Comparison with the reference M. xanthus coupled sheets ("myxo";
%   post-run, zero solves).
%
%   [i, j] = SAG_MYXO_COMPARE('cell', index)        lattice index [kR kB] at
%       step 0.15 -> myxo sheet cell (row i = 13 + kB indexes dn_B, column
%       j = 13 + kR indexes dn_R); NaN outside the 25 x 25 sheet
%   P = SAG_MYXO_COMPARE('ours_on_myxo', sol, xg)
%       our return evaluated at myxo's xg: 8 x numel(xg), rows mapped
%       myxo prof 1..8 -> ours [1 3 5 7 9 11 13 15]; .shift = rows 17:20
%       (myxo shift order sRL sBL sRR sBR)
%   d = SAG_MYXO_COMPARE('distance', Pours, Pmyxo)
%       max over the 8 fields of max_x|ours - myxo| / max_x|myxo|
%   c = SAG_MYXO_COMPARE('class', d)  'agree' (d <= 1e-6) | 'band' |
%       'different' (d >= 1e-2) | 'not_computed' (NaN)
%   c = SAG_MYXO_COMPARE('band_decision', d_ours_tight, d_myxo_tight)
%       'same_solution_myxo_accuracy' if d(ours, ours_tight) <= 1e-6
%       and d(myxo, ours_tight) < 1e-2; 'same_solution_our_accuracy' if
%       d(ours, ours_tight) > 1e-6 and d(myxo, ours_tight) <= 1e-6;
%       'different' otherwise
%   r = SAG_MYXO_COMPARE('reference', prof_ij, ok_ij)
%       myxo reference usable: ok_grid and nodal min > 0 (negative cells
%       solve the clamped problem)

switch action
    case 'cell'
        k = varargin{1};
        i = 13 + k(2); j = 13 + k(1);
        if i < 1 || i > 25 || j < 1 || j > 25, i = NaN; j = NaN; end
        out = [i j];
    case 'ours_on_myxo'
        sol = varargin{1}; xg = varargin{2};
        Y = deval(sol, xg);
        out = struct('prof', Y([1 3 5 7 9 11 13 15], :), 'shift', sol.y(17:20, 1).');
    case 'distance'
        A = varargin{1}; B = varargin{2};       % ours, myxo
        d = 0;
        for k = 1:size(B, 1)
            s = max(abs(B(k, :)));
            if ~(s > 0) || any(~isfinite(A(k, :))), d = NaN; break; end
            d = max(d, max(abs(A(k, :) - B(k, :))) / s);
        end
        out = d;
    case 'class'
        d = varargin{1};
        if ~isfinite(d), out = 'not_computed';
        elseif d <= 1e-6, out = 'agree';
        elseif d >= 1e-2, out = 'different';
        else, out = 'band';
        end
    case 'band_decision'
        a = varargin{1}; b = varargin{2};
        if a <= 1e-6 && b < 1e-2
            out = 'same_solution_myxo_accuracy';
        elseif a > 1e-6 && b <= 1e-6
            out = 'same_solution_our_accuracy';
        else
            out = 'different';
        end
    case 'reference'
        prof = varargin{1}; ok = varargin{2};
        out = logical(ok) && all(isfinite(prof(:))) && min(prof(:)) > 0;
    otherwise
        error('SAG:MYXO', 'Unknown action %s.', action);
end
end
