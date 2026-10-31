function [res, st] = sa_m3_mc_run(run_dir, extra, varargin)
%SA_M3_MC_RUN  The T28 multi-seed stub run: families X_X (branch A) and Y_X
%   (branch B) cold at k:0,0, k:2,0 and k:-2,1; box [3 1]; rescue off; at
%   k:1,1 every return is forced onto branch A (a lineage of branch B lands
%   there on the anchor lineage's root: a merge).
if nargin < 2, extra = struct(); end
o = struct('families', {{'X_X', 'Y_X'}}, ...
    'targets', struct('list', [0 0; 2 0; -2 1]), 'routes', struct('rescue', false));
f = fieldnames(extra);
for k = 1:numel(f), o.(f{k}) = extra.(f{k}); end
force = @(dd, b) forced(dd);
[res, st] = sa_m3_run(run_dir, [3 1], 'force', force, 'settings', o, varargin{:});
end

function b = forced(dd)
b = '';
if abs(dd(1) - 0.65) < 1e-12 && abs(dd(2) - 0.65) < 1e-12, b = 'A'; end
end
