function out = sag_preseed_classify(action, varargin)
%SAG_PRESEED_CLASSIFY  Packet 2C section 5 pre-seeding classifier (zero-solve).
%
%   d = SAG_PRESEED_CLASSIFY('distance', Ystored, Ytight)
%       max over the rows (fields) of max|Ytight - Ystored| / max|Ystored|
%       (2M section 5.3 form); rows are fields, columns mesh points.
%   c = SAG_PRESEED_CLASSIFY('class', r)
%       r: returned, converged, finite (logical), valid (the full v1.2 record
%       verdict of the tight return, for quarters including the 'mpss'
%       record), type_only (the tight return fails ONLY its type checks:
%       morphology/orientation for hpss|quarter, or the mpss shape checks),
%       d (the classifying distance; NaN when not computed).
%       Classes (section 5 table), in this order:
%         moved          converged and finite, and d >= 1e-2 (another
%                        solution), or the tight return fails only its type
%         confirmed      valid and d <= 1e-6
%         band           valid and 1e-6 < d < 1e-2
%         not_confirmed  anything else (tight solve failed or refused)
%   t = SAG_PRESEED_CLASSIFY('thresholds')   [1e-6 1e-2]
%   Seedable in 2C: confirmed and band (band flagged); moved stops C2 for a
%   lead ruling; not_confirmed is excluded (fail closed, section 5 default).

switch action
    case 'thresholds', out = [1e-6 1e-2];
    case 'distance',   out = distance(varargin{:});
    case 'class',      out = classify(varargin{:});
    otherwise, error('SAG:PRESEED', 'Unknown action %s.', action);
end
end

function d = distance(A, B)
d = 0;
for k = 1:size(A, 1)
    s = max(abs(A(k, :)));
    if ~(s > 0) || any(~isfinite(B(k, :))), d = NaN; return; end
    d = max(d, max(abs(B(k, :) - A(k, :))) / s);
end
end

function c = classify(r)
t = [1e-6 1e-2];
if r.returned && r.converged && r.finite && ((isfinite(r.d) && r.d >= t(2)) || r.type_only)
    c = 'moved';
elseif r.valid && isfinite(r.d) && r.d <= t(1)
    c = 'confirmed';
elseif r.valid && isfinite(r.d) && r.d < t(2)
    c = 'band';
else
    c = 'not_confirmed';
end
end
