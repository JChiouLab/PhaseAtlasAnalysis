function name = sa_run_name(parameters, settings)
%SA_RUN_NAME  Propose a run_dir folder NAME for one parameter set, following
%   the per-parameter folder naming of the half-cell library runs.
%
%   name = SA_RUN_NAME(parameters, settings)
%       <model>_<mode>_<tag><value>_..._<grid>_<sha8>. Every scalar numeric
%       field of `parameters` is shown, in field order, under its own name without
%       underscores, with its value printed as the library folder names
%       print it: %g, '.' -> 'p', '-' -> 'm', '+' dropped. For the fission-
%       yeast parameter set this is exactly the library folder's parameter
%       text (..._L7_MBG14_MP3p444); this file stays model-free, so no
%       per-model tag table is kept here. Fields that would take the
%       parameter part past 150 characters are left out of the name. <grid>
%       is the coupled grid of settings.lattice, printed the same way: n<count>
%       when lattice.points_per_axis is set (n51; one count per conserved
%       protein joined by 'x', n25x99), otherwise step<step> (step0p02;
%       per-axis steps step0p01x0p05); without a settings.lattice field there
%       is no grid part. The caller passes the grid the run will use (the
%       resolved settings). sha8 is the first 8 hex digits of the SHA-256 of
%       the full parameter text (model, mode and every field at %.17g; not
%       the grid), so two parameter sets that print alike still get different
%       names. The name enters no solve, row or key. A proposal only: the
%       caller decides the
%       parent folder; a changed scientific setting is a new run, so reusing
%       a name for other settings is refused by seeded_atlas_run (non-empty
%       run_dir).
%       When settings.import.source_dir is set, the name gains
%       _imp<stamp8>, the first 8 hex digits of the source run's stamp
%       (sa_import 'stamp'), so a same-grid follow-up never collides with
%       its source. With import.source_stopped true (a killed source) the
%       eight digits are those of the SHA-256 of the source's
%       settings_effective.mat instead, which stays the same when the
%       source is resumed later, so the default name of the run does too.

if ~isstruct(parameters) || ~isscalar(parameters)
    error('SA:RUN_NAME', 'parameters must be a scalar struct.');
end
model = settings.model;
mode = 'smoke';
if isfield(settings, 'mode'), mode = settings.mode; end
max_shown = 150;
f = fieldnames(parameters);
shown = {};
full = sprintf('model=%s;mode=%s;', model, mode);
for k = 1:numel(f)
    v = parameters.(f{k});
    if isnumeric(v) || islogical(v)
        full = [full sprintf('%s=%s;', f{k}, mat2str(double(v), 17))]; %#ok<AGROW>
        if isscalar(v)
            t = [strrep(f{k}, '_', '') regexprep(sprintf('%.6g', double(v)), {'\.', '-', '\+'}, {'p', 'm', ''})];
            if numel(strjoin([shown, {t}], '_')) <= max_shown
                shown{end+1} = t; %#ok<AGROW>
            end
        end
    elseif ischar(v)
        full = [full sprintf('%s=%s;', f{k}, v)]; %#ok<AGROW>
    end
end
grid = grid_part(settings);
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(full), 'int8'));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
imp = {};
if isfield(settings, 'import') && isstruct(settings.import) && isfield(settings.import, 'source_dir') && ...
        ~isempty(settings.import.source_dir)
    st = sa_import('name_stamp', settings.import);   % the folder stamp (a stopped source: its settings file)
    imp = {['imp' st(1:8)]};
end
name = regexprep(strjoin([{model, mode}, shown, grid, {h(1:8)}, imp], '_'), '[^A-Za-z0-9_]', '');
end

function g = grid_part(settings)
% {} (no lattice field), {'n<count>'} or {'step<step>'}; values as the
% parameters are printed, per-axis values joined by 'x'
g = {};
if ~isfield(settings, 'lattice') || ~isstruct(settings.lattice), return; end
lat = settings.lattice;
if isfield(lat, 'points_per_axis') && ~isempty(lat.points_per_axis)
    g = {['n' values_text(lat.points_per_axis)]};
elseif isfield(lat, 'step') && ~isempty(lat.step)
    g = {['step' values_text(lat.step)]};
end
end

function t = values_text(v)
t = strjoin(arrayfun(@(x) regexprep(sprintf('%.6g', double(x)), {'\.', '-', '\+'}, {'p', 'm', ''}), ...
    reshape(v, 1, []), 'UniformOutput', false), 'x');
end
