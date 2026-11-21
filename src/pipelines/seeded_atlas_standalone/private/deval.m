function varargout = deval(varargin)
%DEVAL  Private forwarder: deval calls of the package root go to sa_deval.
%   A file in a private/ folder wins over the path for every function in the
%   parent folder, so every deval call in the package root (sa_validity,
%   sa_identity, sa_store_io, sa_guess_mesh, sa_pilot_report,
%   seeded_atlas_compare, ...) reaches sa_deval, a vectorized replica of the
%   toolbox deval that returns the same values, without editing those files.
%   sa_deval falls back to the toolbox deval outside its fast path.
%   tests/ and subfolders still reach the toolbox deval.
[varargout{1:max(1, nargout)}] = sa_deval(varargin{:});
end
