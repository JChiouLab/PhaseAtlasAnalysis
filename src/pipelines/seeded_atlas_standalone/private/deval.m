function varargout = deval(varargin)
%DEVAL  M6 private forwarder (packet v0.2 section 3.2, Tier 1 + Tier 2).
%   A file in a private/ folder wins over the path for every function in the
%   parent folder, so every deval call in the package root (sa_validity,
%   sa_identity, sa_store_io, sa_guess_mesh, sa_pilot_report,
%   seeded_atlas_compare, ...) reaches sa_deval, the exact vectorized
%   replica of the toolbox deval, with no byte change to those files.
%   sa_deval falls back to the toolbox deval outside its fast path.
%   tests/ and subfolders still reach the toolbox deval. Pinned by T61 and
%   recorded in the runtime evidence (seeded_atlas_run).
[varargout{1:max(1, nargout)}] = sa_deval(varargin{:});
end
