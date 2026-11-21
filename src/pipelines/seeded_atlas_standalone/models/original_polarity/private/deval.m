function varargout = deval(varargin)
%DEVAL  Private forwarder: deval calls of this model folder go to sa_deval.
%   Every deval call of the functions in this model folder (including the
%   verbatim copies) reaches sa_deval, a vectorized replica of the toolbox
%   deval that returns the same values; sa_deval falls back to the toolbox
%   deval outside its fast path.
[varargout{1:max(1, nargout)}] = sa_deval(varargin{:});
end
