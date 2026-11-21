function varargout = sag_library(action, varargin)
%SAG_LIBRARY  Read-only reader of the M. xanthus uncoupled libraries (isolated half-cell states).
%
%   lib = SAG_LIBRARY('load', manifest, theta, parameters)
%       manifest = struct('directory', d) plus an optional 'confirmation'
%       (a confirmation file with seedable masks). Unknown fields are
%       refused. Loads the labels, axes, masks and `par` of
%       branch_{HSS,HPSS,MPSS}_latest.mat (not the profiles) and records each
%       file's path, bytes and SHA-256. Fails closed unless:
%       model_id guzzo_noxa; schema uncoupled_branch_library/guzzo_noxa/v1;
%       fields {Rm,Rc,Bm,Bc}; a 12-name state_layout; one model_version
%       across the three files; the 16 named parameters equal field-wise and
%       exactly; res and dn_half_range proven from the axes; dn_R and
%       dn_B equal across types, increasing and mirror-symmetric; Lh = L/2;
%       a res x res existence_mask; the per-type labels (orientation,
%       storage_kind, domain).
%       With a confirmation file, the HPSS and MPSS masks are intersected
%       with its seedable masks (confirmed or band), after its
%       model_version and the library SHA-256 values are checked.
%   lib = SAG_LIBRARY('load', manifest, theta, parameters, types)
%       the same, for the active seed shapes (sag_seed_shape):
%       the three files above are always loaded, in that order; a mode-n
%       file branch_MODE<n>_latest.mat (gnx_build_mode) is loaded (and enters
%       lib.files) only when an active type needs it, with the same label,
%       par, axis, Lh and model_version checks (labels: orientation
%       trough_at_0, storage_kind piece_bvp4c+stitched_xy, domain
%       [0 Lh; 0 L/(2n)]). A missing needed mode file is error SA:SEED:E_FILE,
%       one with run_status 'empty' SA:SEED:E_EMPTY. A type made on load
%       (1-, 2+, n+, n-) uses its base file's mask (after the confirmation
%       intersection). Without types (or with the default types) only the
%       three base files are loaded.
%   lib = SAG_LIBRARY('profiles', lib, types)
%   c = SAG_LIBRARY('cell', lib, theta, d)       d = [dn_R dn_B]
%       Native cell [row col] = [index of dn_B, index of dn_R]; in_range
%       within +-dn_half_range*(1 + 1e-12); offsets in native spacings.
%   h = SAG_LIBRARY('half', lib, theta, type, cell, model_options)
%       Concentrations Rm Rc Bm Bc (4 x guess_nodes) on
%       xg = linspace(0, Lh, guess_nodes), pole -> middle, dispatched on
%       storage_kind: uniform state (HSS; hss_branch 'single_root' requires
%       exactly one root), deval of the bvp4c struct (HPSS, peak at the pole),
%       deval of the quarter struct reflected peak-to-peak (MPSS:
%       half(x) = quarter(|x - L/4|), the writer's rule).
%       Made on load (dispatched on the type's transform):
%       1- (M1m) the HPSS half's samples column-reversed (the guess grid is
%       uniform, so this is the exact mirror image); 2+ (M2p) deval of the
%       MPSS quarter at L/4 - |x - L/4| (peaks at both ends); n+ / n- (M<n>p
%       / M<n>m) deval of the mode-n piece at the fold (sag_seed_shape
%       'fold', piece length L/(2n)): peak-start / trough-start.
%   [g, why] = SAG_LIBRARY('pair_guess', lib, theta, tL, cL, tR, cR, mo)
%       The 24-row cold guess from two halves (sag_coupled24 'pair_guess').
%       Each half is already in its own pole -> middle parameterization, so the
%       right half is not reversed.

switch action
    case 'load',       varargout{1} = do_load(varargin{:});
    case 'profiles',   varargout{1} = do_profiles(varargin{:});
    case 'cell',       varargout{1} = do_cell(varargin{:});
    case 'half',       varargout{1} = do_half(varargin{:});
    case 'pair_guess', [varargout{1:max(nargout, 1)}] = do_pair(varargin{:});
    case 'expected_labels', varargout{1} = expected_labels(varargin{:});
    otherwise, error('SAG:LIBRARY', 'Unknown sag_library action %s.', action);
end
end

% =========================================================================

function E = expected_labels(Lh, L, n)
E = struct();
if nargin >= 3                             % the mode-n file of gnx_build_mode
    E.(sprintf('MODE%d', n)) = struct('orientation', 'trough_at_0', ...
        'storage_kind', 'piece_bvp4c+stitched_xy', 'domain', [0 Lh; 0 L / (2 * n)]);
end
E.HSS  = struct('orientation', 'uniform', 'storage_kind', 'uniform_state', 'domain', [0 Lh]);
E.HPSS = struct('orientation', 'peak_at_0', 'storage_kind', 'bvp4c_struct', 'domain', [0 Lh]);
E.MPSS = struct('orientation', 'interior_peak', 'storage_kind', 'quarter_bvp4c+stitched_xy', ...
    'domain', [0 Lh; 0 L / 4]);
end

function lib = do_load(manifest, theta, p, active) %#ok<INUSL>
id = 'SAG:LIBRARY';
if nargin < 4 || isempty(active), active = {'HSS', 'HPSS', 'MPSS'}; end
if ~(isstruct(manifest) && isscalar(manifest) && isfield(manifest, 'directory'))
    error(id, 'library_manifest must be a struct with a ''directory'' field.');
end
extra = setdiff(fieldnames(manifest), {'directory', 'confirmation'});
if ~isempty(extra)
    error(id, 'library_manifest has unknown field(s): %s.', strjoin(extra', ', '));
end
types = {'HSS', 'HPSS', 'MPSS'};
info_of = cellfun(@(t) sag_seed_shape('resolve', sag_seed_shape('shape_of_type', t)), active, ...
    'UniformOutput', false);
info_of = [info_of{:}];
for b = unique({info_of.base}, 'stable')                 % mode-n files after the three
    if ~any(strcmp(b{1}, types)), types{end+1} = b{1}; end %#ok<AGROW>
end
D = gnx_model('descriptor');
EL = expected_labels(theta.Lh, theta.L);
lib = struct('types', {types}, 'directory', manifest.directory, 'files', [], ...
    'res', NaN, 'mask', struct(), 'B', struct(), 'Lh', theta.Lh, 'L', theta.L, ...
    'axis_row', [], 'axis_col', [], 'model_version', '', 'confirmation', []);
files = struct('type', {}, 'path', {}, 'bytes', {}, 'sha256', {});
axes_ = {};
versions = {};
for t = 1:numel(types)
    ty = types{t};
    f = fullfile(manifest.directory, sprintf('branch_%s_latest.mat', ty));
    if exist(f, 'file') ~= 2
        if t > 3
            error('SA:SEED:E_FILE', ['Library file missing: %s (seed shapes %s need it; ' ...
                'gnx_build_mode with n = %s (max_mode >= %s) builds it).'], f, ...
                strjoin({info_of(strcmp({info_of.base}, ty)).shape}, ', '), ty(5:end), ty(5:end));
        end
        error(id, 'Library file missing: %s.', f);
    end
    info = dir(f);
    files(end+1) = struct('type', ty, 'path', f, 'bytes', info.bytes, ...
        'sha256', sa_sha256(f)); %#ok<AGROW>
    need = {'model_id', 'model_version', 'schema', 'fields', 'state_layout', ...
        'branch_type', 'dn_R', 'dn_B', 'Lh', 'res', 'par', 'existence_mask', ...
        'orientation', 'storage_kind', 'domain'};
    w = whos('-file', f);
    miss = setdiff(need, {w.name});
    if ~isempty(miss)
        error(id, '%s lacks %s; not a guzzo_noxa library.', f, strjoin(miss, ', '));
    end
    if t > 3
        need = [need, {'run_status', 'mode_n', 'piece_length'}]; %#ok<AGROW>
        miss = setdiff(need, {w.name});
        if ~isempty(miss)
            error(id, '%s lacks %s; not a guzzo_noxa mode-n library.', f, strjoin(miss, ', '));
        end
    end
    B = load(f, need{:});
    if t > 3
        nmode = str2double(ty(5:end));
        if strcmp(B.run_status, 'empty')
            error('SA:SEED:E_EMPTY', ['%s has run_status ''empty'' (the builder found no piece): ' ...
                'seed shapes %s cannot be seeded.'], f, ...
                strjoin({info_of(strcmp({info_of.base}, ty)).shape}, ', '));
        end
        if ~isequal(B.mode_n, nmode) || ~isequal(B.piece_length, theta.L / (2 * nmode))
            error(id, '%s: mode_n %g / piece_length %.17g, expected %d / L/(2n).', f, ...
                B.mode_n, B.piece_length, nmode);
        end
        EL = expected_labels(theta.Lh, theta.L, nmode);
    end
    if ~(ischar(B.model_id) && strcmp(B.model_id, D.id))
        error(id, '%s: model_id is not %s.', f, D.id);
    end
    if ~strcmp(B.schema, D.schema)
        error(id, '%s: schema %s, expected %s.', f, B.schema, D.schema);
    end
    if ~isequal(B.fields, D.fields)
        error(id, '%s: fields are not {Rm, Rc, Bm, Bc}.', f);
    end
    if ~(iscell(B.state_layout) && numel(B.state_layout) == 12)
        error(id, '%s: state_layout does not have 12 names.', f);
    end
    if ~strcmp(B.branch_type, ty)
        error(id, '%s: branch_type %s.', f, B.branch_type);
    end
    check_par(B.par, theta, f);
    res = numel(B.dn_R);
    if numel(B.dn_B) ~= res || res ~= theta.res
        error(id, '%s: axis length %d/%d, but the run parameters say res = %d.', ...
            f, numel(B.dn_R), numel(B.dn_B), theta.res);
    end
    if B.dn_R(end) ~= theta.dn_half_range || B.dn_B(end) ~= theta.dn_half_range
        error(id, '%s: axis ends at %.17g, but the run parameters say dn_half_range = %.17g.', ...
            f, B.dn_R(end), theta.dn_half_range);
    end
    if ~isequal(B.Lh, theta.L / 2)
        error(id, '%s: Lh %.17g is not L/2.', f, B.Lh);
    end
    if ~(islogical(B.existence_mask) && isequal(size(B.existence_mask), [res res]))
        error(id, '%s: existence_mask is not a logical res x res.', f);
    end
    e = EL.(ty);
    if ~strcmp(B.orientation, e.orientation) || ~strcmp(B.storage_kind, e.storage_kind) || ...
            ~isequal(B.domain, e.domain)
        error(id, '%s: labels (%s, %s, %s) differ from the expected (%s, %s, %s).', f, ...
            B.orientation, B.storage_kind, mat2str(B.domain), e.orientation, ...
            e.storage_kind, mat2str(e.domain));
    end
    versions{end+1} = B.model_version; %#ok<AGROW>
    lib.mask.(ty) = B.existence_mask;
    lib.B.(ty) = struct('path', f, 'SOL', [], 'SOL_quarter', [], 'n_roots', [], ...
        'storage_kind', B.storage_kind);
    if t > 3, lib.B.(ty).SOL_piece = []; end
    axes_{end+1} = {B.dn_R, B.dn_B}; %#ok<AGROW>
end
if numel(unique(versions)) ~= 1
    error(id, 'The library files carry different model_version values.');
end
for t = 2:numel(axes_)
    if ~isequal(axes_{1}, axes_{t})
        error(id, 'branch_%s and branch_%s do not share one native grid.', types{1}, types{t});
    end
end
[dnR, dnB] = axes_{1}{:};
for a = {dnR, dnB}
    v = a{1};
    if ~(isrow(v) && all(diff(v) > 0)) || max(abs(v + fliplr(v))) > 1e-12 * max(abs(v))
        error(id, ['Native axis is not increasing and symmetric about 0; the discrete ' ...
            'mirror (res+1-i, res+1-j) would not be the reflection.']);
    end
end
lib.axis_col = dnR;                        % columns index dn_R
lib.axis_row = dnB;                        % rows index dn_B
lib.res = numel(dnR);
lib.model_version = versions{1};
if isfield(manifest, 'confirmation') && ~isempty(manifest.confirmation)
    [lib, frec] = apply_confirmation(lib, manifest.confirmation, files);
    files(end+1) = frec;
end
lib.files = files;
for k = 1:numel(info_of)                     % the active types (made on load: base mask)
    ty = info_of(k).type;
    if strcmp(ty, info_of(k).base), continue; end   % a stored type: nothing to do
    bB = lib.B.(info_of(k).base);
    lib.mask.(ty) = lib.mask.(info_of(k).base);
    lib.B.(ty) = struct('path', bB.path, 'SOL', [], 'SOL_quarter', [], 'n_roots', [], ...
        'storage_kind', bB.storage_kind, 'SOL_piece', [], 'base', info_of(k).base, ...
        'transform', info_of(k).transform, 'n', info_of(k).n);
end
if ~isequal(active, {'HSS', 'HPSS', 'MPSS'})
    lib.types = active;                       % the active types (seedability loops over them)
    lib.files_types = types;
end
end

function check_par(par, theta, f)
names = {'kR', 'kRB', 'kB', 'kBB', 'dR', 'dB', 'dBA', 'K', 'KA', 'DM', 'DC', ...
    'r_speed', 'dba_factor', 'L', 'rho_R', 'rho_B'};
for k = 1:numel(names)
    n = names{k};
    if ~isstruct(par) || ~isfield(par, n) || ~isequal(par.(n), theta.(n))
        v = NaN; if isstruct(par) && isfield(par, n), v = par.(n); end
        error('SAG:LIBRARY', ['%s was built at %s = %s, but the run parameters ' ...
            'say %s.'], f, n, num2str(v, 17), num2str(theta.(n), 17));
    end
end
end

function [lib, frec] = apply_confirmation(lib, cfile, files)
id = 'SAG:LIBRARY';
if exist(cfile, 'file') ~= 2
    error(id, 'Confirmation file missing: %s.', cfile);
end
C = load(cfile);
if ~isfield(C, 'schema') || ~strcmp(C.schema, 'sag_preseed_confirmation/1')
    error(id, '%s is not a sag_preseed_confirmation/1 file.', cfile);
end
if ~strcmp(C.model_version, lib.model_version)
    error(id, '%s: model_version differs from the library.', cfile);
end
for ty = {'HPSS', 'MPSS'}
    k = strcmp({files.type}, ty{1});
    if ~strcmp(C.library_sha256.(ty{1}), files(k).sha256)
        error(id, '%s: confirms a different %s library file.', cfile, ty{1});
    end
    m = C.seedable.(ty{1});
    if ~(islogical(m) && isequal(size(m), [lib.res lib.res])) || any(m(:) & ~lib.mask.(ty{1})(:))
        error(id, '%s: %s seedable mask is malformed or exceeds the existence mask.', cfile, ty{1});
    end
    lib.mask.(ty{1}) = lib.mask.(ty{1}) & m;
end
info = dir(cfile);
frec = struct('type', 'confirmation', 'path', cfile, 'bytes', info.bytes, ...
    'sha256', sa_sha256(cfile));
lib.confirmation = struct('path', cfile, 'counts', C.counts);
end

function lib = do_profiles(lib, types)
for t = 1:numel(types)
    ty = types{t};
    B = lib.B.(ty);
    switch B.storage_kind
        case 'uniform_state'
            if isempty(B.SOL)
                S = load(B.path, 'SOL', 'n_roots');
                lib.B.(ty).SOL = S.SOL; lib.B.(ty).n_roots = S.n_roots;
            end
        case 'bvp4c_struct'
            if isempty(B.SOL)
                S = load(B.path, 'SOL'); lib.B.(ty).SOL = S.SOL;
            end
        case 'quarter_bvp4c+stitched_xy'
            if isempty(B.SOL_quarter)
                S = load(B.path, 'SOL', 'SOL_quarter');
                lib.B.(ty).SOL = S.SOL; lib.B.(ty).SOL_quarter = S.SOL_quarter;
            end
        case 'piece_bvp4c+stitched_xy'                 % the mode-n piece
            if isempty(B.SOL_piece)
                S = load(B.path, 'SOL_piece'); lib.B.(ty).SOL_piece = S.SOL_piece;
            end
        otherwise
            error('SAG:LIBRARY', 'Unknown storage_kind %s.', B.storage_kind);
    end
end
end

% =========================================================================

function c = do_cell(lib, theta, d) %#ok<INUSL>
c = struct('in_range', false, 'cell', [NaN NaN], 'offset', [NaN NaN], 'reason', '');
lo = [lib.axis_col(1) lib.axis_row(1)];
hi = [lib.axis_col(end) lib.axis_row(end)];
okR = d(1) >= lo(1) * (1 + 1e-12) && d(1) <= hi(1) * (1 + 1e-12);
okB = d(2) >= lo(2) * (1 + 1e-12) && d(2) <= hi(2) * (1 + 1e-12);
if ~(okR && okB)
    r = {};
    if ~okR, r{end+1} = 'dn_R outside native range'; end
    if ~okB, r{end+1} = 'dn_B outside native range'; end
    c.reason = strjoin(r, '+');
    return
end
[eR, col] = min(abs(lib.axis_col - d(1)));
[eB, row] = min(abs(lib.axis_row - d(2)));
c.in_range = true;
c.cell = [row col];
c.offset = [eB / abs(lib.axis_row(2) - lib.axis_row(1)), ...
            eR / abs(lib.axis_col(2) - lib.axis_col(1))];
end

% =========================================================================

function H = do_half(lib, theta, type, cell, mo)
H = struct('type', type, 'cell', cell, 'available', false, 'reason', '', ...
    'xg', [], 'Y4', [], 'storage', '');
i = cell(1); j = cell(2);
if ~any(strcmp(type, lib.types))
    error('SAG:LIBRARY', 'Unknown half type %s.', type);
end
if i < 1 || j < 1 || i > lib.res || j > lib.res
    H.reason = 'cell index out of range'; return
end
if ~lib.mask.(type)(i, j)
    H.reason = 'no seedable uncoupled half at this cell'; return
end
B = lib.B.(type);
xg = linspace(0, theta.Lh, mo.guess_nodes);
tr = 'none';
if isfield(B, 'transform'), tr = B.transform; end
switch tr                                  % shapes made on load
    case 'rotate'                          % 2+: the MPSS quarter, peaks at both ends
        if isempty(B.SOL_quarter), error('SAG:LIBRARY', 'MPSS profiles are not loaded.'); end
        s = B.SOL_quarter{i, j};
        ellq = theta.L / 4;
        xm = ellq - abs(xg - ellq);
        xm = min(max(xm, s.x(1)), s.x(end));
        Y = deval(s, xm);
        H.xg = xg; H.Y4 = Y([1 3 5 7], :); H.storage = 'quarter_deval_rotated';
        H.available = true;
        return
    case {'mode_P', 'mode_T'}              % n+ / n-: the mode-n piece, P / T-start
        if isempty(B.SOL_piece), error('SAG:LIBRARY', '%s profiles are not loaded.', B.base); end
        s = B.SOL_piece{i, j};
        if isempty(s)
            error('SAG:LIBRARY', '%s SOL_piece is empty at existent cell (%d,%d).', B.base, i, j);
        end
        xm = sag_seed_shape('fold', xg, theta.L / (2 * B.n), tr(end));
        xm = min(max(xm, s.x(1)), s.x(end));
        Y = deval(s, xm);
        H.xg = xg; H.Y4 = Y([1 3 5 7], :); H.storage = sprintf('piece_deval_%s', tr(end));
        H.available = true;
        return
end
switch B.storage_kind
    case 'uniform_state'
        if isempty(B.SOL), error('SAG:LIBRARY', 'HSS profiles are not loaded.'); end
        if ~strcmp(mo.hss_branch, 'single_root')
            error('SAG:LIBRARY', 'hss_branch %s is not implemented.', mo.hss_branch);
        end
        if B.n_roots(i, j) ~= 1
            H.reason = sprintf('HSS cell has %d roots; hss_branch single_root fails closed', ...
                B.n_roots(i, j));
            return
        end
        u = B.SOL{i, j};
        if ~(isnumeric(u) && numel(u) == 4)
            error('SAG:LIBRARY', 'HSS state at (%d,%d) is not a 4-vector.', i, j);
        end
        Y4 = repmat(u(:), 1, numel(xg));
        H.storage = 'uniform_state';
    case 'bvp4c_struct'
        if isempty(B.SOL), error('SAG:LIBRARY', 'HPSS profiles are not loaded.'); end
        s = B.SOL{i, j};
        xq = min(max(xg, s.x(1)), s.x(end));     % round-off only
        Y = deval(s, xq);
        Y4 = Y([1 3 5 7], :);
        H.storage = 'bvp4c_deval';
        if strcmp(tr, 'reverse')           % 1-: the HPSS half's samples column-reversed
            Y4 = Y4(:, end:-1:1);
            H.storage = 'bvp4c_deval_reversed';
        end
    case 'quarter_bvp4c+stitched_xy'
        if isempty(B.SOL_quarter), error('SAG:LIBRARY', 'MPSS profiles are not loaded.'); end
        s = B.SOL_quarter{i, j};
        ellq = theta.L / 4;
        xm = abs(xg - ellq);
        xm = min(max(xm, s.x(1)), s.x(end));
        Y = deval(s, xm);
        Y4 = Y([1 3 5 7], :);
        H.storage = 'quarter_deval_reflected';
end
H.xg = xg;
H.Y4 = Y4;
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
g = sag_coupled24('pair_guess', HL.xg, HL.Y4, HR.Y4);
end
