function P = sa_post_model(id, main_variable)
%SA_POST_MODEL  Model-specific post-processing for plots and full-PDE traces.
%
%   P = SA_POST_MODEL(id) puts post/<id>/ on the path when it exists and
%   returns the model's post struct from its factory (original_polarity:
%   sap_post; guzzo_noxa: sag_post; any other id: sa_post_<id>, which must
%   be on the path, as test stubs are). post/ is outside the run's code
%   hashes: editing it never blocks a resume; its comment-free token digest
%   (sa_code_digest) enters the trace and figure keys instead.
%   P = SA_POST_MODEL(id, main_variable): the third coordinate dn2 of that
%   main variable (sa_main_variable; '' = the module's default). Only a
%   module whose factory takes the argument has a main variable (then P
%   also carries main_variable, legacy_main_variable and
%   legacy_dn2_variable, the label's and dn2's variables before software
%   0.10); asking another module for one is refused (SA:POST:E_MAIN_VARIABLE).
%
%   Members of P
%     id, n_fields, field_names           the half-domain fields
%     sol_rows                            2 x n_fields: the solver rows of the
%                                         left / right half concentrations
%     profile_rows                        2 x n_fields: the same rows in the
%                                         stored identity profile
%     coord_labels                        {c1 c2 c3} axis labels (plot order)
%     lattice_axes                        [a1 a2]: lattice axis of c1 and c2
%     mass_pairs                          k x 2 field pairs (a, b) whose sum is
%                                         one pool (the cos perturbation)
%     full_profile(theta, xi, YL, YR, nx) -> struct x, eta, Y, im: the full
%                                         domain from the two halves (left at
%                                         x = xi, right at x = L - xi; pchip;
%                                         the seam node from the left half)
%     observables(theta, x, Y)            -> struct c (1 x 3), masses (1 x 2),
%                                         umax (1 x 2, the first field per
%                                         half)
%     trace_ic(theta, sol, entry)         -> struct x, Y (after the
%                                         perturbation), Y_root (before)
%     integrate(theta, Y0, entry)         -> frames struct (see seeded_atlas_trace)
%     entry_defaults()                    the trace entry defaults
%     plot_nx                             mesh of the plot's coordinates
%     aspect(c)                           DataAspectRatio for coordinates c
%     files                               the post files of this model (digests)

here = fileparts(mfilename('fullpath'));
if ~(ischar(id) && isrow(id) && ~isempty(regexp(id, '^[a-z][a-z0-9_]*$', 'once')))
    error('SA:POST:E_ID', 'A model id must be lower-case text.');
end
folder = fullfile(here, id);
if exist(folder, 'dir') == 7 && ~contains([pathsep path pathsep], [pathsep folder pathsep])
    addpath(folder);
end
switch id
    case 'original_polarity', fname = 'sap_post';
    case 'guzzo_noxa',        fname = 'sag_post';
    otherwise,                fname = ['sa_post_' id];
end
if isempty(which(fname))
    error('SA:POST:E_NOT_FOUND', 'No post module for model %s (%s is not on the path).', id, fname);
end
if nargin < 2, main_variable = ''; end
if nargin(fname) >= 1
    P = feval(fname, main_variable);
elseif isempty(main_variable)
    P = feval(fname);
else
    error('SA:POST:E_MAIN_VARIABLE', 'The post module of model %s has no main variable (asked for %s).', ...
        id, main_variable);
end
P.files = [{fullfile(here, 'sa_post_model.m'), fullfile(here, 'sa_code_digest.m'), ...
    fullfile(fileparts(here), 'sa_main_variable.m'), fullfile(fileparts(here), 'sa_mode_name.m')}, P.files];
P.family_order = @family_order;
P.family_colour = @family_colour;
end

% =========================================================================

function o = family_order(f)
% The legend order: the family names written before software 0.11
% (HPSS_HSS, ...) in their fixed order; then mode names (sa_mode_name) by
% total mode (the two halves' mode numbers added, 0 for a flat half), then
% by name; then any other label alphabetically.
fixed = {'HPSS_HSS', 'HSS_HPSS', 'HSS_HSS', 'HPSS_HPSS', 'HSS_MPSS', 'MPSS_HSS', 'HPSS_MPSS', ...
    'MPSS_HPSS', 'MPSS_MPSS'};
rest = setdiff(f, fixed);                        % sorted by name
tot = nan(1, numel(rest));
for k = 1:numel(rest)
    n = sa_mode_name('parse', rest{k});
    if ~isempty(n), tot(k) = sum(n); end
end
m = ~isnan(tot);
modes = rest(m);
[~, j] = sortrows([tot(m).', (1:nnz(m)).']);    % ties keep the name order
o = [fixed(ismember(fixed, f)), modes(j), rest(~m)];
end

function c = family_colour(name)
% One colour per family name, the same in every model. A mode name that
% corresponds to a family name written before software 0.11 keeps that
% family's colour (1+/0 = HPSS_HSS, 0/1- = HSS_HPSS, 0/0 = HSS_HSS, the
% three Fig 3 superfamily colours; 2-/0 = MPSS_HSS, 0/2- = HSS_MPSS, 1+/2- =
% HPSS_MPSS, 2-/1- = MPSS_HPSS, 1+/1- = HPSS_HPSS, 2-/2- = MPSS_MPSS); a
% single peak at the cell middle (1-/1+) is gold, 'unlabelled' grey; any
% other name gets a generated colour (generated_colour).
switch name
    case {'HPSS_HSS', '1+/0'},   h = '2166ac';
    case {'HSS_HPSS', '0/1-'},   h = 'e66101';
    case {'HSS_HSS', '0/0'},     h = '00a79d';
    case {'HPSS_HPSS', '1+/1-'}, h = '7b3294';
    case {'HSS_MPSS', '0/2-'},   h = 'b2182b';
    case {'MPSS_HSS', '2-/0'},   h = 'f4a582';
    case {'HPSS_MPSS', '1+/2-'}, h = '1b7837';
    case {'MPSS_HPSS', '2-/1-'}, h = 'a6dba0';
    case {'MPSS_MPSS', '2-/2-'}, h = '525252';
    case '1-/1+',                h = 'c9a227';
    case 'unlabelled',           h = 'bdbdbd';
    otherwise, c = generated_colour(name); return
end
c = hex(h);
end

function c = generated_colour(name)
% A colour fixed by the name alone (the same in every run and session): the
% hue from the SHA-256 of the name, saturation 0.65, value 0.75; the hue is
% stepped by the golden angle until the colour lies at least 0.25 (RGB
% distance) from every fixed colour above, grey included.
fixed = {'2166ac', 'e66101', '00a79d', '7b3294', 'b2182b', 'f4a582', '1b7837', 'a6dba0', '525252', ...
    'c9a227', 'bdbdbd'};
F = cell2mat(cellfun(@hex, fixed(:), 'UniformOutput', false));
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(name), 'int8'));
d = double(typecast(md.digest(), 'uint8'));
h = (d(1) * 256 + d(2)) / 65536;
for k = 1:64
    c = hsv2rgb([h 0.65 0.75]);
    if min(sqrt(sum((F - c) .^ 2, 2))) >= 0.25, return; end
    h = mod(h + 0.381966, 1);
end
end

function c = hex(h)
c = hex2dec({h(1:2), h(3:4), h(5:6)}).' / 255;
end
