function out = gnx_model(action, varargin)
%GNX_MODEL  M. xanthus model descriptor and parameter normalization.
%
%   d = GNX_MODEL('descriptor')
%       Model contract: id, fields, BVP state layout, conserved quantities,
%       diffusion, base rates, numerics defaults and validity tolerances.
%   P = GNX_MODEL('parameters', s)
%       Normalizes a user struct. s must carry exactly the fields
%       L, r_speed, dba_factor, rho_R, rho_B, res, dn_half_range
%       (unknown or missing fields are an error, id GNX:MODEL). Returns the
%       normalized struct: the eleven rates kR..DC with the scalings applied
%       (dBA = dba_factor*dBA, kR = r_speed*kR, kRB = r_speed*kRB,
%       dR = r_speed*dR), followed by r_speed, dba_factor,
%       L, rho_R, rho_B, res, dn_half_range. The half-domain is L/2 (not stored).
%   v = GNX_MODEL('version', extra_dirs)
%       SHA-256 over the sorted gnx_*.m files of this folder and of the
%       optional extra folders (cell array): the model_version label.
%   h = GNX_MODEL('sha256_file', path) ;  h = GNX_MODEL('sha256_bytes', uint8)
%
% This folder is shared by the uncoupled and the coupled pipelines.

switch action
    case 'descriptor', out = descriptor();
    case 'parameters', out = normalize(varargin{:});
    case 'version',    out = model_version(varargin{:});
    case 'sha256_file',  out = sha256_file(varargin{1});
    case 'sha256_bytes', out = sha256_bytes(varargin{1});
    otherwise, error('GNX:MODEL', 'Unknown gnx_model action %s.', action);
end
end

% -------------------------------------------------------------------------
function d = descriptor()
d = struct();
d.id = 'guzzo_noxa';
d.schema = 'uncoupled_branch_library/guzzo_noxa/v1';
d.fields = {'Rm', 'Rc', 'Bm', 'Bc'};
d.state_layout = {'Rm', 'dRm', 'Rc', 'dRc', 'Bm', 'dBm', 'Bc', 'dBc', ...
                  'sigmaR', 'sigmaB', 'mR', 'mB'};
d.conc_rows = [1 3 5 7];          % BVP rows of Rm Rc Bm Bc
d.deriv_rows = [2 4 6 8];
d.sigma_rows = [9 10];
d.acc_rows = [11 12];
d.C = [1 1 0 0; 0 0 1 1];
d.base = base_rates();            % base rates (before the r_speed and dba_factor scalings)
d.D_names = {'DM', 'DC', 'DM', 'DC'};
d.param_fields_required = {'L', 'r_speed', 'dba_factor', 'rho_R', 'rho_B', ...
                           'res', 'dn_half_range'};
d.axis_names = {'dn_R', 'dn_B'};  % columns index dn_R, rows index dn_B (density differences, see README)
% numerics (bvp4c settings and mesh sizes)
d.numerics = struct('n_init', 121, 'RelTol', 1e-8, 'AbsTol', 1e-6, ...
    'NMax', 20000, 'n_dense', 20001, 'guard_factor', 0.5, ...
    'morph_nodes', 121);
% tolerances of the validity checks
d.tol = struct('mass_rel', 1e-7, 'hss_route_rel', 1e-12, 'hss_resid_rel', 1e-12, ...
    'hss_imag', 1e-10, 'contrast', 0.2, 'mpss_slope', 1e-2, 'mpss_sym', 1e-6, ...
    'mpss_peak_lo', 0.4, 'mpss_peak_hi', 0.6, 'acc_end_rel', 1e-12);
% thresholds of the profile-shape (morphology) classification
d.morphology = struct('amplitude_rel', 0.2, 'edge_fraction', 0.10, ...
    'trough_rel_tol', 0.02, 'prominence_rel', 0.05, 'absolute_rel', 1e-3, ...
    'autocorr_min', 0.5, 'negative_tol', 0);
d.sweep = struct('max_calls', 20000, 'amp_fracs', [0.9 0.7 0.5], ...
    'widths', [0.15 0.3 0.5 1.0 1.5], 'rescue_passes', 3);
end

function b = base_rates()
% Base rates (min^-1; concentrations in units of the mean total densities). With
% r_speed = 0.35 and dba_factor = 1 the scalings of normalize() give Table 1:
% kr, krb, dr = 0.35*kR, 0.35*kRB, 0.35*dR (0.0100, 0.365, 0.0332 before rounding).
b = struct();
b.kR  = 0.028624;
b.kRB = 1.0423;
b.kB  = 8.1255;
b.kBB = 3.5216;
b.dR  = 0.094962;
b.dB  = 5.4381;
b.dBA = 210.01;
b.K   = 0.11440;
b.KA  = 0.14714;
b.DM  = 0.2;
b.DC  = 120.0;
end

% -------------------------------------------------------------------------
function P = normalize(s)
id = 'GNX:MODEL';
if ~(isstruct(s) && isscalar(s))
    error(id, 'parameters must be a scalar struct.');
end
req = {'L', 'r_speed', 'dba_factor', 'rho_R', 'rho_B', 'res', 'dn_half_range'};
f = fieldnames(s);
unknown = setdiff(f, req);
missing = setdiff(req, f);
if ~isempty(unknown)
    error(id, 'Unknown parameter field(s): %s.', strjoin(unknown', ', '));
end
if ~isempty(missing)
    error(id, 'Missing parameter field(s): %s.', strjoin(missing, ', '));
end
for k = 1:numel(req)
    v = s.(req{k});
    if ~(isnumeric(v) && isreal(v) && isscalar(v) && isfinite(v))
        error(id, 'Parameter %s must be a finite real scalar.', req{k});
    end
end
if s.L <= 0 || s.r_speed <= 0 || s.dba_factor < 0 || s.rho_R <= 0 || s.rho_B <= 0
    error(id, 'Need L, r_speed, rho_R, rho_B > 0 and dba_factor >= 0.');
end
if s.res < 3 || s.res ~= round(s.res)
    error(id, 'res must be an integer >= 3.');
end
if s.dn_half_range <= 0 || s.dn_half_range >= 2*min(s.rho_R, s.rho_B)
    error(id, 'dn_half_range must lie in (0, 2*min(rho_R, rho_B)) so every half mass is positive.');
end
b = base_rates();
P = b;
% scalings: dBA by dba_factor (the signal factor f), kR, kRB, dR by r_speed
P.dBA = s.dba_factor * b.dBA;
P.kR  = s.r_speed * b.kR;
P.kRB = s.r_speed * b.kRB;
P.dR  = s.r_speed * b.dR;
P.r_speed = double(s.r_speed);
P.dba_factor = double(s.dba_factor);
P.L = double(s.L);
P.rho_R = double(s.rho_R);
P.rho_B = double(s.rho_B);
P.res = double(s.res);
P.dn_half_range = double(s.dn_half_range);
P = orderfields(P, [fieldnames(b); {'r_speed'; 'dba_factor'; 'L'; 'rho_R'; ...
    'rho_B'; 'res'; 'dn_half_range'}]);
end

% -------------------------------------------------------------------------
function v = model_version(extra_dirs)
if nargin < 1, extra_dirs = {}; end
dirs = [{fileparts(mfilename('fullpath'))}, extra_dirs(:)'];
lines = {};
for k = 1:numel(dirs)
    L = dir(fullfile(dirs{k}, 'gnx_*.m'));
    for q = 1:numel(L)
        lines{end+1} = sprintf('%s %s', L(q).name, ...
            sha256_file(fullfile(L(q).folder, L(q).name))); %#ok<AGROW>
    end
end
lines = sort(lines);
v = sha256_bytes(uint8(strjoin(lines, newline)));
end

function h = sha256_file(p)
fid = fopen(p, 'r');
if fid < 0, error('GNX:MODEL', 'Cannot open %s.', p); end
c = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
h = sha256_bytes(bytes);
end

function h = sha256_bytes(bytes)
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes(:));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2)', 1, []));
end
