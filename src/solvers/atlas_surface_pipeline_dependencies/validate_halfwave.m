% validate_halfwave -- wrapper for the half-peak (HPSS) morphology check.
% Shared by the half-peak sweeper and the full-peak constructor; keep a single classifier.
function [ok, info] = validate_halfwave(sol, opts)
    if nargin < 2, opts = struct(); end
    info = classify_hpss_morphology(sol, opts);
    ok = info.accepted;
end
