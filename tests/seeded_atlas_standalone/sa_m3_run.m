function [res, st] = sa_m3_run(run_dir, box, varargin)
%SA_M3_RUN  One stub run (sa_m3_stub settings) with state name/value pairs;
%   'settings' (a struct merged into the stub settings) is taken out first.
s = sa_m3_stub('settings', run_dir, box);
k = find(strcmp(varargin(1:2:end), 'settings'));
if ~isempty(k)
    o = varargin{2*k};
    varargin(2*k-1:2*k) = [];
    s = merge(s, o);
end
st = sa_m3_stub('state', varargin{:});
p = sa_m3_stub('params');
if isfield(s, 'params_override'), p = s.params_override; s = rmfield(s, 'params_override'); end
res = seeded_atlas_run(p, struct(), s, sa_m3_stub('corrector', st));
end

function a = merge(a, b)
f = fieldnames(b);
for k = 1:numel(f)
    if isstruct(b.(f{k})) && isfield(a, f{k}) && isstruct(a.(f{k}))
        a.(f{k}) = merge(a.(f{k}), b.(f{k}));
    else
        a.(f{k}) = b.(f{k});
    end
end
end
