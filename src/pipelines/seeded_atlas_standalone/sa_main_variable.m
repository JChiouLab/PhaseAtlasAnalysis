function varargout = sa_main_variable(action, varargin)
%SA_MAIN_VARIABLE  The "main variable": a sum of a model's fields.
%
%   The main variable m is the membrane variable that the shape label (a mode
%   name per half, sa_mode_name) and the third plot coordinate
%   dn2 = int_0^1 cos(2 pi eta) m d(eta) look at; the defaults are u for the
%   fission yeast model and Bm for the M. xanthus model. On the rescaled
%   cell eta = x/L. The paper's Delta n2 = L*dn2 for the M. xanthus model and
%   M_BG*dn2 for the fission yeast model (whose concentrations are scaled so
%   that total GEF is 1). It is written as a
%   '+'-separated sum of the model's
%   field names, e.g. 'u', 'q' or 'u + q'; whitespace is allowed, every
%   field at most once, no other text. Model-owned setting
%   model_options.main_variable (sa_model_<id> validates it with 'parse').
%
%   v = SA_MAIN_VARIABLE('parse', text, names, what)
%       text: the expression; names: the model's field names (cell row);
%       what: the setting's name for messages (default
%       'model_options.main_variable'). v.text is the canonical form
%       (whitespace removed, the given order kept), v.names the field names
%       in that order, v.index their positions in names. Anything else is
%       refused (SA:SETTINGS:E_FIELD) with the allowed field names.
%   m = SA_MAIN_VARIABLE('sum', Y, rows)
%       Y(rows(1), :) + Y(rows(2), :) + ... in the given order (one row:
%       the row itself, no arithmetic, so a single field is bit-identical
%       to reading the row).
%   s = SA_MAIN_VARIABLE('file_tag', text)
%       The text for file names: 'mv-' followed by the canonical text
%       (e.g. mv-u+q).

switch action
    case 'parse',    varargout{1} = parse(varargin{:});
    case 'sum',      varargout{1} = row_sum(varargin{:});
    case 'file_tag', varargout{1} = ['mv-' varargin{1}];
    otherwise, error('SA:MAIN_VARIABLE', 'Unknown sa_main_variable action %s.', action);
end
end

% =========================================================================

function v = parse(text, names, what)
if nargin < 3, what = 'model_options.main_variable'; end
id = 'SA:SETTINGS:E_FIELD';
allowed = strjoin(names, ', ');
if isstring(text) && isscalar(text), text = char(text); end
if ~(ischar(text) && (isrow(text) || isempty(text)))
    error(id, ['%s must be text: a ''+''-separated sum of this model''s fields (%s), ' ...
        'e.g. ''%s'' or ''%s+%s''.'], what, allowed, names{1}, names{1}, names{2});
end
t = regexprep(text, '\s', '');
parts = strsplit(t, '+');
if isempty(t) || any(cellfun(@isempty, parts))
    error(id, ['%s ''%s'' is not a ''+''-separated sum of this model''s fields (%s), ' ...
        'e.g. ''%s'' or ''%s+%s''.'], what, text, allowed, names{1}, names{1}, names{2});
end
[ok, idx] = ismember(parts, names);
if ~all(ok)
    error(id, '%s ''%s'': %s is not a field of this model (allowed: %s).', what, text, ...
        strjoin(parts(~ok), ', '), allowed);
end
if numel(unique(idx)) ~= numel(idx)
    error(id, '%s ''%s'' names a field more than once (each field at most once; allowed: %s).', ...
        what, text, allowed);
end
v = struct('text', strjoin(parts, '+'), 'names', {parts}, 'index', idx);
end

function m = row_sum(Y, rows)
m = Y(rows(1), :);
for k = 2:numel(rows)
    m = m + Y(rows(k), :);
end
end
