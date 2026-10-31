function varargout = deval(varargin)
%DEVAL  M6 private forwarder (packet v0.2 section 3.2, Tier 1 + Tier 2).
%   Every deval call of the functions in this model folder (the verbatim
%   copies included, with no byte change) reaches sa_deval, the exact
%   vectorized replica of the toolbox deval; sa_deval falls back to the
%   toolbox deval outside its fast path. Pinned by T61 and recorded in the
%   runtime evidence (seeded_atlas_run).
[varargout{1:max(1, nargout)}] = sa_deval(varargin{:});
end
