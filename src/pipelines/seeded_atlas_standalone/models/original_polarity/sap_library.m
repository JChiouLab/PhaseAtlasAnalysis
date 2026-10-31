function varargout = sap_library(action, varargin)
%SAP_LIBRARY  Read-only reader of the original model's uncoupled seed library.
%
%   Generalizes the campaign adapter sca_ad_b (no typed L, rho or M_BG: every
%   scale comes from theta, and the library's own `par` must equal the user's
%   parameters).
%
%   lib = SAP_LIBRARY('load', manifest, theta, parameters)
%       manifest.directory holds branch_{HSS,HPSS,MPSS}_latest.mat. Loads the
%       grids, existence masks and `par` of every type (not the profiles) and
%       records each file's identity (path, bytes, SHA-256). A missing file,
%       a parameter mismatch, unequal grids or an asymmetric grid is an error:
%       a library problem never silently degrades to library-free search.
%   lib = SAP_LIBRARY('load', manifest, theta, parameters, types)
%       (milestone 8) the same, for the active seed types (sap_seed_shape):
%       the three files above are always loaded, in that order; a mode-n
%       file branch_MODE<n>_latest.mat is loaded (and enters lib.files) only
%       when an active type needs it, with the same par and grid checks.
%       A type made on load (1-, 2+, n+, n-) uses its base file's mask.
%       Without types (or with the default types) lib is as in 0.11.
%   lib = SAP_LIBRARY('profiles', lib, types)
%       Loads the stored half profiles (SOL) of the listed types.
%   c = SAP_LIBRARY('cell', lib, theta, d)
%       Native cell of lattice coordinate d = [d_P d_G]: in_range, cell
%       [row col] (nearest index on each axis), offsets in native spacings.
%   h = SAP_LIBRARY('half', lib, theta, type, cell, model_options)
%       One half profile on xi = linspace(0, 1/2, guess_nodes), normalized.
%   [g, why] = SAP_LIBRARY('pair_guess', lib, theta, tL, cL, tR, cR, mo)
%       The 28-row cold guess from two halves (copied corrector pair_guess).
%
%   Native grid: rows index GEF_DIFF (dimensional d_G * M_BG), columns index
%   PAK_DIFF (dimensional d_P * M_BG); both coordinates divide by M_BG
%   (the campaign ruling). Per-type storage differs and is detected, not
%   assumed: HPSS = 14-row bvp4c struct (deval), MPSS = 5-row {x,y} (pchip),
%   HSS = uniform 1x5 state of the representative branch.

switch action
    case 'load',       varargout{1} = do_load(varargin{:});
    case 'profiles',   varargout{1} = do_profiles(varargin{:});
    case 'cell',       varargout{1} = do_cell(varargin{:});
    case 'half',       varargout{1} = do_half(varargin{:});
    case 'pair_guess', [varargout{1:max(nargout,1)}] = do_pair(varargin{:});
    otherwise, error('SAP:LIBRARY', 'Unknown sap_library action %s.', action);
end
end

% =========================================================================

function lib = do_load(manifest, theta, p, active)
id = 'SAP:LIBRARY';
if nargin < 4 || isempty(active), active = {'HSS', 'HPSS', 'MPSS'}; end
if ~(isstruct(manifest) && isscalar(manifest) && isfield(manifest, 'directory'))
    error(id, 'library_manifest must be a struct with a ''directory'' field.');
end
extra = setdiff(fieldnames(manifest), {'directory'});
if ~isempty(extra)
    error(id, 'library_manifest has unknown field(s): %s.', strjoin(extra', ', '));
end
types = {'HSS', 'HPSS', 'MPSS'};
info_of = cellfun(@(t) sap_seed_shape('resolve', sap_seed_shape('shape_of_type', t)), active, ...
    'UniformOutput', false);
info_of = [info_of{:}];
for b = unique({info_of.base}, 'stable')                 % mode-n files after the three
    if ~any(strcmp(b{1}, types)), types{end+1} = b{1}; end %#ok<AGROW>
end
lib = struct('types', {types}, 'directory', manifest.directory, 'files', [], ...
    'res', NaN, 'mask', struct(), 'B', struct(), 'C0', theta.C0, 'L', theta.L, ...
    'M_BG', theta.M_BG, 'axis_row', [], 'axis_col', []);
files = struct('type', {}, 'path', {}, 'bytes', {}, 'sha256', {});
grids = {};
for t = 1:numel(types)
    f = fullfile(manifest.directory, sprintf('branch_%s_latest.mat', types{t}));
    if exist(f, 'file') ~= 2
        if t > 3
            error('SA:SEED:E_FILE', ['Library file missing: %s (seed shapes %s need it; ' ...
                'atlas_build_libraries with opts.max_mode >= %s builds it).'], f, ...
                strjoin({info_of(strcmp({info_of.base}, types{t})).shape}, ', '), types{t}(5:end));
        end
        error(id, 'Library file missing: %s.', f);
    end
    info = dir(f);
    files(end+1) = struct('type', types{t}, 'path', f, 'bytes', info.bytes, ...
        'sha256', sa_sha256(f)); %#ok<AGROW>
    if t > 3
        B = load(f, 'existence_mask', 'PAK_DIFF', 'GEF_DIFF', 'par', 'res', 'run_status');
        if strcmp(B.run_status, 'empty')
            error('SA:SEED:E_EMPTY', ['%s has run_status ''empty'' (the builder found no piece): ' ...
                'seed shapes %s cannot be seeded.'], f, ...
                strjoin({info_of(strcmp({info_of.base}, types{t})).shape}, ', '));
        end
    else
        B = load(f, 'existence_mask', 'PAK_DIFF', 'GEF_DIFF', 'par', 'res');
    end
    check_par(B.par, p, f);
    if ~isequal(size(B.existence_mask), [B.res B.res])
        error(id, '%s: existence_mask is not res x res.', f);
    end
    lib.mask.(types{t}) = logical(B.existence_mask);
    lib.B.(types{t}) = struct('path', f, 'SOL', []);
    grids{end+1} = {B.PAK_DIFF, B.GEF_DIFF, B.res}; %#ok<AGROW>
end
for t = 2:numel(grids)
    if ~isequal(grids{1}, grids{t})
        error(id, 'branch_%s and branch_%s do not share one native grid.', types{1}, types{t});
    end
end
for k = 1:numel(info_of)                     % the active types (made on load: base mask)
    ty = info_of(k).type;
    if strcmp(ty, info_of(k).base), continue; end   % a stored type: as in 0.11
    lib.mask.(ty) = lib.mask.(info_of(k).base);
    lib.B.(ty) = struct('path', lib.B.(info_of(k).base).path, 'SOL', [], ...
        'base', info_of(k).base, 'transform', info_of(k).transform, 'n', info_of(k).n);
end
if ~isequal(active, {'HSS', 'HPSS', 'MPSS'})
    lib.types = active;                       % the active types (seedability loops over them)
    lib.files_types = types;
end
[PAK, GEF, res] = grids{1}{:};
if any(any(diff(PAK, 1, 1))) || any(any(diff(GEF, 1, 2)))
    error(id, 'Native grid is not a tensor grid (PAK along columns, GEF along rows).');
end
lib.axis_col = PAK(1, :);           % dimensional d_P * M_BG, indexed by column
lib.axis_row = GEF(:, 1).';         % dimensional d_G * M_BG, indexed by row
for ax = {lib.axis_col, lib.axis_row}
    a = ax{1};
    if any(diff(a) <= 0) || max(abs(a + fliplr(a))) > 1e-12 * max(abs(a))
        error(id, ['Native axis is not increasing and symmetric about 0; the ' ...
            'discrete mirror (res+1-i, res+1-j) would not be the reflection.']);
    end
end
lib.res = res;
lib.files = files;
end

function check_par(par, p, f)
for k = {'a','b','c','d','e','f','Dm','Dc','L','M_BG','M_P'}
    if ~isfield(par, k{1}) || ~isequal(par.(k{1}), p.(k{1}))
        error('SAP:LIBRARY', ['%s was built at %s = %s, but the run parameters ' ...
            'say %s.'], f, k{1}, num2str(getf(par, k{1}), 17), num2str(p.(k{1}), 17));
    end
end
end

function v = getf(s, f)
if isfield(s, f), v = s.(f); else, v = NaN; end
end

function lib = do_profiles(lib, types)
for t = 1:numel(types)
    ty = types{t};
    if isempty(lib.B.(ty).SOL)
        if isfield(lib.B.(ty), 'transform') && any(strcmp(lib.B.(ty).transform, {'mode_P', 'mode_T'}))
            S = load(lib.B.(ty).path, 'SOL_piece');     % the mode-n piece {x, y}
            lib.B.(ty).SOL = S.SOL_piece;
        else
            S = load(lib.B.(ty).path, 'SOL');
            lib.B.(ty).SOL = S.SOL;
        end
    end
end
end

% =========================================================================

function c = do_cell(lib, theta, d)
dim = d * theta.M_BG;                      % [d_P d_G] -> dimensional
lo = [lib.axis_col(1) lib.axis_row(1)];
hi = [lib.axis_col(end) lib.axis_row(end)];
c = struct('in_range', false, 'cell', [NaN NaN], 'offset', [NaN NaN], 'reason', '');
okP = dim(1) >= lo(1) * (1 + 1e-12) && dim(1) <= hi(1) * (1 + 1e-12);
okG = dim(2) >= lo(2) * (1 + 1e-12) && dim(2) <= hi(2) * (1 + 1e-12);
if ~(okP && okG)
    r = {};
    if ~okP, r{end+1} = 'd_P outside native range'; end
    if ~okG, r{end+1} = 'd_G outside native range'; end
    c.reason = strjoin(r, '+');
    return
end
[eP, col] = min(abs(lib.axis_col - dim(1)));
[eG, row] = min(abs(lib.axis_row - dim(2)));
c.in_range = true;
c.cell = [row col];
c.offset = [eG / abs(lib.axis_row(2) - lib.axis_row(1)), ...
            eP / abs(lib.axis_col(2) - lib.axis_col(1))];
end

% =========================================================================

function H = do_half(lib, theta, type, cell, mo)
H = struct('type', type, 'cell', cell, 'available', false, 'reason', '', ...
    'xi', [], 'Y5', [], 'storage', '');
i = cell(1); j = cell(2);
if ~any(strcmp(type, lib.types))
    error('SAP:LIBRARY', 'Unknown half type %s.', type);
end
if i < 1 || j < 1 || i > lib.res || j > lib.res
    H.reason = 'cell index out of range'; return
end
if ~lib.mask.(type)(i, j)
    H.reason = 'no uncoupled half at this cell'; return
end
SOL = lib.B.(type).SOL;
if isempty(SOL)
    error('SAP:LIBRARY', 'Profiles of type %s are not loaded (sap_library ''profiles'').', type);
end
xi = linspace(0, 0.5, mo.guess_nodes);
x_dim = xi * theta.L;                      % xi = x / L (full domain)
s = SOL{i, j};
tr = 'none';
if isfield(lib.B.(type), 'transform'), tr = lib.B.(type).transform; end
switch tr                                  % milestone 8: shapes made on load
    case 'rotate'                          % 2+: the stored stitch rotated by a quarter
        m = (numel(s.x) + 1) / 2;
        Lq = s.x(m);
        s = struct('x', [s.x(m:end) - Lq, s.x(2:m) + Lq], 'y', [s.y(:, m:end), s.y(:, 2:m)]);
    case {'mode_P', 'mode_T'}              % n+ / n-: the piece stitched P / T-start
        [xs, ys] = sap_seed_shape('stitch', s.x, s.y, lib.B.(type).n, tr(end));
        s = struct('x', xs, 'y', ys);
end
if isempty(s)
    error('SAP:LIBRARY', '%s SOL is empty at existent cell (%d,%d).', type, i, j);
end
if isnumeric(s)                            % HSS: uniform representative state
    if numel(s) ~= 5
        error('SAP:LIBRARY', 'HSS state at (%d,%d) has %d entries, expected 5.', i, j, numel(s));
    end
    Y5_dim = repmat(s(:), 1, numel(xi));
    H.storage = 'uniform_state';
else
    xq = min(max(x_dim, s.x(1)), s.x(end));   % clamp round-off only
    if isfield(s, 'solver') && ~isempty(s.solver)
        Y = deval(s, xq);
        H.storage = 'bvp4c_deval';
    else
        Y = interp1(s.x(:), s.y.', xq(:), 'pchip').';
        H.storage = 'xy_pchip';
    end
    switch size(Y, 1)
        case 14, Y5_dim = Y([1 3 5 7 9], :);
        case 5,  Y5_dim = Y;
        otherwise
            error('SAP:LIBRARY', 'Stored %s profile has %d rows; expected 14 or 5.', type, size(Y, 1));
    end
end
if strcmp(tr, 'reverse')                   % 1-: the HPSS half's samples column-reversed
    Y5_dim = Y5_dim(:, end:-1:1);
end
H.xi = xi;
H.Y5 = Y5_dim / theta.C0;
H.available = true;
end

function [g, why] = do_pair(lib, theta, tL, cL, tR, cR, mo)
HL = do_half(lib, theta, tL, cL, mo);
HR = do_half(lib, theta, tR, cR, mo);
g = []; why = '';
if ~HL.available || ~HR.available
    r = {};
    if ~HL.available, r{end+1} = sprintf('left %s(%d,%d): %s', tL, cL, HL.reason); end
    if ~HR.available, r{end+1} = sprintf('right %s(%d,%d): %s', tR, cR, HR.reason); end
    why = strjoin(r, '; ');
    return
end
g = sap_forced28('pair_guess', HL.xi, HL.Y5, HR.Y5, []);
end
