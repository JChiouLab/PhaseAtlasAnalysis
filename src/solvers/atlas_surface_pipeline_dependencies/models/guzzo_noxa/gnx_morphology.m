function [accepted, info, view, opts] = gnx_morphology(x, y4, mt)
%GNX_MORPHOLOGY  M. xanthus (Guzzo_noXA) view into classify_hpss_morphology (used unchanged).
%
%   [accepted, info, view, opts] = GNX_MORPHOLOGY(x, y4)
%   [accepted, info, view, opts] = GNX_MORPHOLOGY(x, y4, mt)
%   x is 1 x n, y4 is 4 x n with rows Rm Rc Bm Bc (field order). The view
%   struct handed to classify_hpss_morphology has 7 rows: row 1 = Bm (the
%   peak field), rows 3, 5, 7 = Rm, Rc, Bc, rows 2, 4, 6 = 0, so its
%   concentration rows are exactly the four fields.
%   mt holds the per-model thresholds (default
%   gnx_model('descriptor').morphology); per-profile values:
%     amplitude_tol = amplitude_rel * max(Bm)   (relative contrast of the reference M. xanthus code)
%     absolute_tol  = absolute_rel  * max(Bm)
%   accepted = (info.class == "HPSS"). Orientation is not decided here
%   (classify flips internally); the caller checks peak-at-0 separately.

if nargin < 3 || isempty(mt)
    d = gnx_model('descriptor');
    mt = d.morphology;
end
x = x(:).';
n = numel(x);
view = struct('x', x, 'y', zeros(7, n));
view.y(1, :) = y4(3, :);    % Bm
view.y(3, :) = y4(1, :);    % Rm
view.y(5, :) = y4(2, :);    % Rc
view.y(7, :) = y4(4, :);    % Bc
bmax = max(y4(3, :));
opts = struct('amplitude_tol', mt.amplitude_rel * bmax, ...
    'edge_fraction', mt.edge_fraction, ...
    'trough_rel_tol', mt.trough_rel_tol, ...
    'prominence_rel', mt.prominence_rel, ...
    'absolute_tol', mt.absolute_rel * bmax, ...
    'autocorr_min', mt.autocorr_min, ...
    'negative_tol', mt.negative_tol);
info = classify_hpss_morphology(view, opts);
accepted = info.class == "HPSS";
end
