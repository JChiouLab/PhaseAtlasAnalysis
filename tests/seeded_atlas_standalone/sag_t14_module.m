function rows = sag_t14_module(E)
%SAG_T14_MODULE  sa_model('guzzo_noxa') validates; file naming and shadowing;
%   the shared core resolves from src/models/guzzo_noxa and is unchanged in
%   content (no gnx_* file added; model_version still the libraries'); the
%   two copies equal their sources except the function-name line; the
%   uncoupled pipeline's noxa folder is not on the path.
T = 'T14 module';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
M = sa_model('guzzo_noxa');
rows = sa_check(rows, T, 'sa_model(''guzzo_noxa'') resolves and validates (closed member list)', ...
    strcmp(M.id, 'guzzo_noxa') && M.n_states == 24 && M.unclipped, M.version);
L = dir(fullfile(E.module, '*.m'));
names = {L.name};
bad = names(~startsWith(names, 'sag_') & ~strcmp(names, 'sa_model_guzzo_noxa.m'));
rows = sa_check(rows, T, 'every module file is sag_* (besides the factory)', isempty(bad), ...
    sprintf('%d files: %s', numel(names), strjoin(names, ', ')));
shadow = {};
for k = 1:numel(names)
    [~, nm] = fileparts(names{k});
    w = which(nm, '-all');
    if numel(w) ~= 1 || ~startsWith(w{1}, E.module), shadow{end+1} = nm; end %#ok<AGROW>
end
rows = sa_check(rows, T, 'no module file shadows or is shadowed (which -all: one entry, in the module)', ...
    isempty(shadow), strjoin(shadow, ', '));
wm = which('gnx_model'); wr = which('gnx_reaction');
rows = sa_check(rows, T, 'shared core gnx_model / gnx_reaction resolve from src/models/guzzo_noxa', ...
    startsWith(wm, E.core) && startsWith(wr, E.core) && numel(which('gnx_model', '-all')) == 1 && ...
    numel(which('gnx_reaction', '-all')) == 1, sprintf('%s; %s', wm, wr));
g = dir(fullfile(E.core, 'gnx_*.m'));
rows = sa_check(rows, T, 'no gnx_* file added to the shared core', ...
    isequal(sort({g.name}), {'gnx_jacobian.m', 'gnx_model.m', 'gnx_reaction.m'}), strjoin({g.name}, ', '));
% M8 amended (by path): milestone 8 adds gnx_build_mode.m (the mode-n builder, packet v0.2 section 3.2)
% to the uncoupled noxa folder, which changes the folder fingerprint; the libraries' fingerprint is
% recomputed over the folder's gnx_* files without it (a scratch copy; content-only hash).
fp = fullfile(E.scratch, sprintf('t14_fp_%s', datestr(now, 'yyyymmdd_HHMMSSFFF'))); mkdir(fp); %#ok<DATST,TNOW1>
gf = dir(fullfile(E.solver_noxa, 'gnx_*.m'));
for k = 1:numel(gf)
    if ~strcmp(gf(k).name, 'gnx_build_mode.m'), copyfile(fullfile(gf(k).folder, gf(k).name), fp); end
end
v = gnx_model('version', {fp});
rows = sa_check(rows, T, 'recomputed model_version (without the M8 gnx_build_mode.m) equals the libraries'' 5d5cc282d4492314...', ...
    startsWith(v, E.model_version), v);
pe = strsplit(path, pathsep);
rows = sa_check(rows, T, 'the uncoupled pipeline''s noxa folder is not on the path', ...
    ~any(strcmp(pe, E.solver_noxa)) && isempty(which('gnx_solve')) && isempty(which('gnx_hermite_min')), ...
    strjoin(pe(contains(pe, 'noxa')), ' | '));
src = {fullfile(E.solver_noxa, 'gnx_hermite_min.m'), fullfile(E.original_module, 'sap_bvp_convergence.m')};
cpy = {fullfile(E.module, 'sag_hermite_min.m'), fullfile(E.module, 'sag_bvp_convergence.m')};
pinned = {'10b230f7d82fcadbe81f8f369132ead5b80fa1614d94209915cd3804638ae39a', ...
          '791036185a281e2aedacb7bc0a466006eb25f3ee6317ec5b3e46c4f9557a8437'};
for k = 1:2
    a = strsplit(fileread(src{k}), newline); b = strsplit(fileread(cpy{k}), newline);
    same_rest = numel(a) == numel(b) && isequal(a(2:end), b(2:end));
    [~, sn] = fileparts(src{k}); [~, cn] = fileparts(cpy{k});
    first_ok = strcmp(strrep(a{1}, sn, cn), b{1});
    rows = sa_check(rows, T, sprintf('copy %s = %s except the function-name line; source SHA-256 pinned', cn, sn), ...
        same_rest && first_ok && strcmp(sa_sha256(src{k}), pinned{k}), ...
        sprintf('source sha %s', sa_sha256(src{k})));
end
cf = {M.copies.file};
rows = sa_check(rows, T, 'module copies record the two copied files with their SHA-256', ...
    isequal(cf, {'sag_hermite_min.m', 'sag_bvp_convergence.m'}) && ...
    strcmp(M.copies(1).copy_sha256, sa_sha256(cpy{1})), strjoin(cf, ', '));
end
