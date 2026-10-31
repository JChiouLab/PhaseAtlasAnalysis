function varargout = sap_seed_shape(action, varargin)
%SAP_SEED_SHAPE  The original model's seed shapes (milestone 8, packet section 4).
%
%   A seed shape is a half's mode name in the LEFT half's drawn orientation,
%   which is the stored pole-to-seam orientation (sa_mode_name): '0', '1+',
%   '2-' are today's library types; '1-' and '2+' are new shapes made on load
%   from the HPSS and MPSS files; 'n+' / 'n-' (n >= 3) come from the mode-n
%   file branch_MODE<n>_latest.mat (build_mode_branch). The right half seeds
%   the same stored shape; its drawn name is the mirror of the shape.
%
%   i = SAP_SEED_SHAPE('resolve', shape)
%       .shape .type (type id) .base (file type) .file .transform .n .sign
%       transform: 'none' | 'reverse' (1-: the HPSS guess samples column-
%       reversed) | 'rotate' (2+: the stored MPSS stitch rotated by a quarter)
%       | 'mode_P' / 'mode_T' (n+ / n-: the mode-n piece stitched P-start /
%       T-start). Unknown names: error SA:SEED:E_SHAPE.
%   s = SAP_SEED_SHAPE('shape_of_type', type)
%   t = SAP_SEED_SHAPE('order', types)      types sorted in the type order
%       (HSS, HPSS, MPSS, M1m, M2p, M3p, M3m, M4p, M4m, ...)
%   m = SAP_SEED_SHAPE('mirror', shape)     the drawn name of a right half
%       seeded with shape: the sign flips for odd n and stays for even n
%   [shapes, skipped] = SAP_SEED_SHAPE('all', directory)
%       every shape the library folder offers, in type order; skipped: the
%       MODE files whose run_status is 'empty' (path, status)
%   f = SAP_SEED_SHAPE('mode_families', types)   the label -> type map built
%       from the active types (left: shape -> type; right: mirror(shape) ->
%       type); at the default types it equals the module's 0.11 map
%   [x, y] = SAP_SEED_SHAPE('stitch', xq, yq, n, start)   n copies of a
%       piece (peak at its left end) on [0, n*xq(end)]: start 'T' (n-, a
%       trough at the pole) or 'P' (n+, a peak at the pole); a copy whose
%       right end is a peak is the mirrored piece, x = (k+1) l - fliplr(xq);
%       otherwise x = k l + xq; the duplicate join node (each copy's last)
%       is dropped. At n = 2, 'T' is build_mpss_branch's reflect_stitch.
%   r = SAP_SEED_SHAPE('cold_rate')  measured cold seconds per call (L = 7, n101)
%   d = SAP_SEED_SHAPE('defaults')   the default seed_modes {'0','1+','2-'}

switch action
    case 'resolve',       varargout{1} = resolve(varargin{:});
    case 'shape_of_type', varargout{1} = shape_of_type(varargin{:});
    case 'order',         varargout{1} = order(varargin{:});
    case 'mirror',        varargout{1} = mirror(varargin{:});
    case 'all',           [varargout{1}, varargout{2}] = all_shapes(varargin{:});
    case 'mode_families', varargout{1} = mode_families(varargin{:});
    case 'stitch',        [varargout{1}, varargout{2}] = stitch(varargin{:});
    case 'cold_rate',     varargout{1} = 3.7;
    case 'defaults',      varargout{1} = {'0', '1+', '2-'};
    otherwise, error('SA:SEED:E_ACTION', 'Unknown sap_seed_shape action %s.', action);
end
end

% =========================================================================

function i = resolve(shape)
[n, sg] = parse(shape);
i = struct('shape', shape, 'type', '', 'base', '', 'file', '', 'transform', 'none', ...
    'n', n, 'sign', sg);
switch shape
    case '0',  i.type = 'HSS';  i.base = 'HSS';
    case '1+', i.type = 'HPSS'; i.base = 'HPSS';
    case '2-', i.type = 'MPSS'; i.base = 'MPSS';
    case '1-', i.type = 'M1m';  i.base = 'HPSS'; i.transform = 'reverse';
    case '2+', i.type = 'M2p';  i.base = 'MPSS'; i.transform = 'rotate';
    otherwise
        i.base = sprintf('MODE%d', n);
        if strcmp(sg, '+')
            i.type = sprintf('M%dp', n); i.transform = 'mode_P';
        else
            i.type = sprintf('M%dm', n); i.transform = 'mode_T';
        end
end
i.file = sprintf('branch_%s_latest.mat', i.base);
end

function [n, sg] = parse(shape)
if ~(ischar(shape) && isrow(shape))
    error('SA:SEED:E_SHAPE', 'A seed shape must be a mode name such as ''0'', ''1+'', ''2-'' or ''3+''.');
end
if strcmp(shape, '0'), n = 0; sg = ''; return; end
tok = regexp(shape, '^([1-9][0-9]*)([+-])$', 'tokens', 'once');
if isempty(tok)
    error('SA:SEED:E_SHAPE', ['Unknown seed shape ''%s'': a shape is ''0'' or ''n+'' / ''n-'' ' ...
        '(n >= 1, no leading zero).'], shape);
end
n = str2double(tok{1}); sg = tok{2};
end

function s = shape_of_type(type)
switch type
    case 'HSS',  s = '0';
    case 'HPSS', s = '1+';
    case 'MPSS', s = '2-';
    case 'M1m',  s = '1-';
    case 'M2p',  s = '2+';
    otherwise
        tok = regexp(type, '^M([1-9][0-9]*)([pm])$', 'tokens', 'once');
        if isempty(tok) || str2double(tok{1}) < 3
            error('SA:SEED:E_SHAPE', 'Unknown seed type %s.', type);
        end
        if strcmp(tok{2}, 'p'), s = [tok{1} '+']; else, s = [tok{1} '-']; end
end
end

function t = order(types)
key = zeros(numel(types), 1);
for k = 1:numel(types)
    switch types{k}
        case 'HSS',  key(k) = 1;
        case 'HPSS', key(k) = 2;
        case 'MPSS', key(k) = 3;
        case 'M1m',  key(k) = 4;
        case 'M2p',  key(k) = 5;
        otherwise
            s = shape_of_type(types{k});
            [n, sg] = parse(s);
            key(k) = 10 * n + strcmp(sg, '-');
    end
end
[~, o] = sort(key);
t = types(o);
end

function m = mirror(shape)
[n, sg] = parse(shape);
m = shape;
if n > 0 && mod(n, 2) == 1
    if strcmp(sg, '+'), m = sprintf('%d-', n); else, m = sprintf('%d+', n); end
end
end

function [shapes, skipped] = all_shapes(directory)
shapes = {'0', '1+', '2-', '1-', '2+'};
skipped = struct('file', {}, 'status', {});
f = dir(fullfile(directory, 'branch_MODE*_latest.mat'));
ns = zeros(1, 0);
for k = 1:numel(f)
    tok = regexp(f(k).name, '^branch_MODE([1-9][0-9]*)_latest\.mat$', 'tokens', 'once');
    if ~isempty(tok), ns(end+1) = str2double(tok{1}); end %#ok<AGROW>
end
for n = sort(ns)
    p = fullfile(directory, sprintf('branch_MODE%d_latest.mat', n));
    Z = load(p, 'run_status');
    if strcmp(Z.run_status, 'empty')
        skipped(end+1) = struct('file', p, 'status', Z.run_status); %#ok<AGROW>
    else
        shapes = [shapes, {sprintf('%d+', n), sprintf('%d-', n)}]; %#ok<AGROW>
    end
end
end

function f = mode_families(types)
L = cell(numel(types), 2); R = cell(numel(types), 2);
for k = 1:numel(types)
    s = shape_of_type(types{k});
    L(k, :) = {s, types{k}};
    R(k, :) = {mirror(s), types{k}};
end
f = struct('left', {L}, 'right', {R});
end

function [x, y] = stitch(xq, yq, n, start)
xq = xq(:).';
l = xq(end);
x = zeros(1, 0); y = zeros(size(yq, 1), 0);
for k = 0:n-1
    if strcmp(start, 'T'), mirrored = mod(k, 2) == 0; else, mirrored = mod(k, 2) == 1; end
    if mirrored
        xk = (k + 1) * l - fliplr(xq); yk = fliplr(yq);
    else
        xk = k * l + xq; yk = yq;
    end
    if k < n - 1
        xk = xk(1:end-1); yk = yk(:, 1:end-1);
    end
    x = [x, xk]; %#ok<AGROW>
    y = [y, yk]; %#ok<AGROW>
end
end
