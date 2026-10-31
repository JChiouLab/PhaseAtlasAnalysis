% validate_halfwave -- compatibility wrapper for the shared HPSS morphology gate.
% Shared by the HPSS sweeper and MPSS constructor; do not fork classifier logic.
function [ok, info] = validate_halfwave(sol, opts)
    if nargin < 2, opts = struct(); end
    info = classify_hpss_morphology(sol, opts);
    ok = info.accepted;
end
