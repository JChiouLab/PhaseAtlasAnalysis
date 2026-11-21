function rows = sa_t11_stub_model(E)
%SA_T11_STUB_MODEL  The whole shared path (settings, lattice,
%   seedability, targets, validity, root store, outputs, refinement) runs on a
%   stub second model with 4 fields, 2 pools, an 8-row state, domain [0, 1.5]
%   and shifts as bvp4c parameters. Any hard-coded 28/10/[0,0.5]/y(21:24) in
%   shared code fails here. Also checks the runtime path holds no
%   data/experiments folder, both in the test process and in the runtime
%   evidence the run records in settings_effective.mat. Zero solves: the corrector is a stub.
T = 'T11 stub model';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
p = path;
entries = strsplit(p, pathsep);
rows = sa_check(rows, T, 'runtime MATLAB path has no data/experiments folder', ...
    ~any(contains(entries, fullfile('data', 'experiments'))), ...
    sprintf('%d path entries checked', numel(entries)));
M = sa_model('stub_two_pool');
rows = sa_check(rows, T, 'stub module passes the closed member check', strcmp(M.id, 'stub_two_pool'));

d = fullfile(E.scratch, sprintf('t11_stub_run_%s', datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<DATST,TNOW1>
settings = struct('model', 'stub_two_pool', 'run_dir', d, ...
    'lattice', struct('step', 0.2), 'budget', struct('max_calls', 60));
params = struct('R', 1, 'B', 2, 'Lh', 1.5);
res = seeded_atlas_run(params, struct(), settings, @stub_corrector);

rows = sa_check(rows, T, 'run completes', strcmp(res.stop_reason, 'completed'), ...
    sprintf('%d targets, %d calls, %d roots, census %s', res.n_targets, res.calls, res.n_roots, ...
    strjoin(cellfun(@(n, c) sprintf('%s=%d', n, c), res.census.names, num2cell(res.census.counts), ...
    'UniformOutput', false), ' ')));
A = res.attempts;
cold = A(strcmp({A.kind}, 'cold'));
rows = sa_check(rows, T, 'one cold call per declared target', numel(cold) == res.n_targets);
bad = cold(strcmp({cold.family}, 'Y_X') & strcmp({cold.target_key}, 'k:0,0'));
rows = sa_check(rows, T, 'poisoned return (profile dips below zero) refused for positivity', ...
    numel(bad) == 1 && strcmp(bad.outcome, 'invalid') && contains(bad.failing, 'positivity'), ...
    sprintf('%s', strjoin({bad.failing}, '|')));
ok = all(arrayfun(@(r) size(r.representative_profile, 1) == 4 && ...
    r.comparison_mesh(end) == 1.5, res.store.roots));
rows = sa_check(rows, T, 'roots hold 4 identity rows on the module''s [0, 1.5] comparison mesh', ...
    ok && res.n_roots > 0);
F = load(fullfile(sa_layout('state', d), cold(1).file));
rows = sa_check(rows, T, 'attempt file keeps the 8-row returned state and parameter shifts', ...
    size(F.sol.y, 1) == 8 && numel(F.sol.parameters) == 4 && isfield(F.sol, 'yp'));
ref = res.refinement;
rows = sa_check(rows, T, 'refinement ran per family and reports distance 0 for an identical return', ...
    ~isempty(ref) && all(strcmp({ref.status}, 'refined')) && all([ref.scaled_rms] == 0), ...
    sprintf('%d refinement rows', numel(ref)));
% roots are stored as per-node shards roots/<node>.mat, and
% the run folder holds journal.txt and queue.mat
need = {'settings_effective.mat','attempts.mat','attempts.csv','summary.txt','identity_pairs.csv', ...
    'journal.txt','queue.mat'};
have = cellfun(@(f) exist(fullfile(d, f), 'file') == 2 || exist(fullfile(sa_layout('state', d), f), 'file') == 2, need);   % state files are in internal/
rows = sa_check(rows, T, 'all outputs written', all(have) && exist(fullfile(sa_layout('state', d), 'roots'), 'dir') == 7, ...
    strjoin(need(~have), ','));
Sh = sa_store_io('load_run_shards', sa_layout('state', d), M);
nr = sum(cellfun(@(k) numel(Sh(k).roots), Sh.keys));
rows = sa_check(rows, T, 'root shards reload and re-validate (re-derivation)', nr == res.n_roots);
Eff = load(fullfile(sa_layout('state', d), 'settings_effective.mat'));
rows = sa_check(rows, T, 'effective settings saved with software and model versions', ...
    isfield(Eff.settings, 'software_version') && strcmp(Eff.model.id, 'stub_two_pool'));
% The summary counts cross-route identity from root MEMBERSHIP
% (roots whose discoveries span >= 2 routes; different-route discovery pairs
% inside one root); refinement returns are never discoveries. Cold-only: 0.
P = res.pairs;
tol = Eff.settings.identity.similarity_tolerance;
nref0 = sum(strcmp({P.kind}, 'refinement') & [P.scaled_rms] <= tol);
nmulti = 0; npairs = 0;
for k = Sh.keys
    st = Sh(k{1});
    for r = 1:numel(st.roots)
        rt = {st.roots(r).discoveries.route};
        nmulti = nmulti + (numel(unique(rt)) >= 2);
        for a = 1:numel(rt), for b = a+1:numel(rt), npairs = npairs + ~strcmp(rt{a}, rt{b}); end, end
    end
end
txt = fileread(fullfile(d, 'summary.txt'));
% the summary wording carries no development labels
tok = regexp(txt, 'cross-route identity by root membership: (\d+) root\(s\)[^;]*; (\d+) different-route', 'tokens', 'once');
rows = sa_check(rows, T, ['summary membership counts = independent counts ' ...
    '(refinement pairs excluded); tolerance remains provisional; no "is met"'], ...
    ~isempty(tok) && str2double(tok{1}) == nmulti && str2double(tok{2}) == npairs && ...
    contains(txt, 'remains PROVISIONAL') && ~contains(txt, 'is met'), ...
    sprintf('summary %s, independent %d/%d, refinement pairs at <= tol %d', ...
    strjoin(tok, '/'), nmulti, npairs, nref0));
Rt = Eff.runtime;
pe = Rt.path_entries;
hit = contains(pe, fullfile('data', 'experiments')) | contains(pe, 'coupled_half_domain_atlas');
rows = sa_check(rows, T, ['run records runtime evidence (path, pwd, which) and no recorded ' ...
    'path entry contains data/experiments or coupled_half_domain_atlas'], ...
    iscellstr(pe) && numel(pe) > 1 && ~any(hit) && ischar(Rt.pwd) && ~isempty(Rt.pwd) && ...
    endsWith(Rt.which_seeded_atlas_run, fullfile('seeded_atlas_standalone', 'seeded_atlas_run.m')) && ...
    Rt.corrector_injected && strcmp(Rt.corrector_name, 'stub_corrector') && ...
    endsWith(Rt.corrector_file, 'sa_t11_stub_model.m'), ...
    sprintf('%d entries, %d offending; pwd %s; which %s; corrector %s (%s)', numel(pe), nnz(hit), ...
    Rt.pwd, Rt.which_seeded_atlas_run, Rt.corrector_name, Rt.corrector_file));
got = '';
try, seeded_atlas_run(params, struct(), settings, @stub_corrector); catch ME, got = ME.identifier; end
rows = sa_check(rows, T, 'a non-empty run_dir is refused (no silent overwrite)', strcmp(got, 'SA:RUN'), got);
rmpath(E.stub);
end

function [sol, info] = stub_corrector(th, guess, target, opts) %#ok<INUSD>
k = pi / th.Lh;
x = linspace(0, th.Lh, 2001);         % fine enough for the residual gate
A = target.n(:);
Bamp = 0.2 * ones(4, 1);
if all(target.d == 0) && guess.y(1, 1) == 2            % left type Y at the anchor
    Bamp(2) = 3 * A(2);                                  % profile dips below zero
end
y = zeros(8, numel(x)); yp = y;
for j = 1:4
    y(2*j-1, :) = A(j) + Bamp(j) * cos(k * x);
    y(2*j, :) = -Bamp(j) * k * sin(k * x);
    yp(2*j-1, :) = y(2*j, :);
    yp(2*j, :) = -k^2 * Bamp(j) * cos(k * x);
end
sol = struct('solver', 'bvp4c', 'x', x, 'y', y, 'yp', yp, 'parameters', A, ...
    'stats', struct('maxres', 1e-9, 'nmeshpoints', numel(x)));
info = struct('class', 'STUB', 'warning_id', '');
end
