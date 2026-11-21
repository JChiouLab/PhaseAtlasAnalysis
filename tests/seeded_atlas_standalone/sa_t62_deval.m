function rows = sa_t62_deval(E, varargin)
%SA_T62_DEVAL  sa_deval is bit-identical to the toolbox deval. Zero solves.
%   rows = SA_T62_DEVAL(E)
%     - planted cases on synthetic one-region bvp4c structs: nodes, t(end),
%       duplicates, unsorted and decreasing xint, a single point, a column
%       xint, row subsets (repeated, reordered), slope requested or not, -0,
%       a tiny interval; each also asserted to stay on the FAST path (the
%       profiler shows no toolbox deval call);
%     - fallback cases: every listed unsupported input gives the
%       same output (uint64, class, size) or the same error identifier as
%       deval, and goes through the toolbox (profiler); the missing-solver
%       case reproduces the toolbox warning (lastwarn) only through the
%       fallback; the interface case its NonuniqueSolution warning;
%     - the fast path returns double and leaves lastwarn untouched.
%   X = SA_T62_DEVAL('sweep', spec)       the stored-record sweep
%     spec.runs: struct array (name, dir, limit) of run folders; for every
%     attempt file with a returned iterate, at every point set the pipeline
%     uses: positivity points (mesh nodes, interval fractions, the n_dense
%     grid; field rows), the n_dense grid with slopes (all rows), the
%     n_dense grid on the identity rows (sa_identity dense), the comparison
%     mesh (identity rows), interval midpoints with slopes (sag_coupled24),
%     the shape-label grid, the clamped reconstruction grids (sap_forced28),
%     the subsample meshes at N = 401 / 482 and the secant form
%     (another record's subsample mesh), and a random permutation (unsorted);
%     spec.libraries: struct array (name, model, params, manifest): every
%     bvp4c struct of the loaded library at the library xq sets (uniform,
%     clamped; the reflected quarter grid). uint64 equality incl. NaN and
%     -0; any differing bit is a failure. spec.out: csv path (per record).
if ischar(E) && strcmp(E, 'sweep')
    rows = sweep(varargin{:});
    return
end
T = 'T62 sa_deval';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
rng(62, 'twister');
% ---- planted cases (fast path) ------------------------------------------------------------
s = synth(41, 6, [0 1.5]);
x = s.x;
mid = (x(1:end-1) + x(2:end)) / 2;
cases = {
    'nodes',                 {s, x},                          2
    't(end) alone',          {s, x(end)},                     2
    'nodes and t(end), rows 2:3', {s, x, [2 3]},              2
    'duplicates',            {s, [x(3) x(3) mid(5) mid(5) x(end) x(end)]}, 2
    'unsorted',              {s, mid(randperm(numel(mid)))},  2
    'decreasing',            {s, fliplr(sort([x mid]))},      2
    'single interior point', {s, mid(7)},                     2
    'column xint',           {s, mid(:)},                     2
    'row subset repeated',   {s, mid, [3 1 3]},               2
    'value only',            {s, [x mid]},                    1
    'dense grid',            {s, linspace(x(1), x(end), 2001)}, 2
    'minus zero',            {s, [-0 mid(1) 0]},              2
    };
st = synth(12, 3, [0 1]); st.x(6) = st.x(5) + 1e-13;            % a tiny interval
cases(end+1, :) = {'tiny h', {st, [st.x(5) st.x(5) + 4e-14 st.x(6) (st.x(6) + st.x(7)) / 2]}, 2};
sz = synth(9, 4, [0 2]); sz.yp(:, 3) = -0; sz.y(:, 4) = -0;     % -0 slopes and values at nodes
cases(end+1, :) = {'-0 slopes and values', {sz, [sz.x (sz.x(1:end-1) + sz.x(2:end)) / 2]}, 2};
okp = true; det = {};
for k = 1:size(cases, 1)
    [eq, fast, why] = compare(cases{k, 2}, cases{k, 3});
    okp = okp && eq && fast;
    if ~(eq && fast), det{end+1} = sprintf('%s: eq %d fast %d %s', cases{k, 1}, eq, fast, why); end %#ok<AGROW>
end
rows = sa_check(rows, T, sprintf(['planted cases (%d: nodes, t(end), duplicates, unsorted, decreasing, single point, ' ...
    'column, row subsets, value only, dense, -0, tiny h, -0 slopes): bit-identical to deval (uint64, class, size) ' ...
    'on the fast path'], size(cases, 1)), okp, strjoin(det, ' | '));
% ---- fallback cases ------------------------------------------------------------------
nos = rmfield(s, 'solver');
dd = synth(20, 3, [0 1]); dd.solver = 'dde23';
b5 = s; b5.solver = 'bvp5c';
ifc = s; ifc.x = [s.x(1:20) s.x(20) s.x(21:end-1)];               % a repeated mesh point (interface)
ifc.y = [s.y(:, 1:20) s.y(:, 20) + 0.5 s.y(:, 21:end-1)]; ifc.yp = [s.yp(:, 1:20) s.yp(:, 20) s.yp(:, 21:end-1)];
dec = s; dec.x = fliplr(s.x); dec.y = fliplr(s.y); dec.yp = fliplr(s.yp);
col = s; col.x = s.x(:);
sgl = s; sgl.x = single(s.x);
nn = s; nn.idata = struct('idxNonNegative', 2);
fb = {
    'solver other than bvp4c (dde23)', {dd, mid(1:5) / 1.5}
    'solver bvp5c without idata',      {b5, mid(1:3)}
    'struct without a solver field',   {nos, mid}
    'swapped (xint, sol) order',       {mid, s}
    'empty xint',                      {s, zeros(1, 0)}
    'idx 0',                           {s, mid, 0}
    'idx above n',                     {s, mid, 7}
    'idx negative',                    {s, mid, -1}
    'repeated mesh point',             {ifc, [ifc.x(20) mid(3)]}
    'decreasing x',                    {dec, mid}
    'column sol.x',                    {col, mid}
    'single sol.x',                    {sgl, mid}
    'idata.idxNonNegative',            {nn, mid}
    'single xint',                     {s, single(mid)}
    'sparse xint',                     {s, sparse(mid)}
    'NaN xint',                        {s, [mid(1) NaN]}
    'out-of-range xint',               {s, [mid(1) x(end) + 1]}
    'matrix xint',                     {s, reshape(mid(1:6), 2, 3)}
    'logical idx',                     {s, mid, logical([1 0 1 0 0 0])}
    };
okf = true; det = {};
for k = 1:size(fb, 1)
    [eq, fast, why] = compare(fb{k, 2}, 2);
    if ~(eq && ~fast), okf = false; det{end+1} = sprintf('%s: eq %d fast %d %s', fb{k, 1}, eq, fast, why); end %#ok<AGROW>
end
% the warning cases: lastwarn after the call is the toolbox's
lastwarn('', ''); try, deval(nos, mid); catch, end; [~, w0] = lastwarn;
lastwarn('', ''); try, sa_deval(nos, mid); catch, end; [~, w1] = lastwarn;
lastwarn('', ''); try, deval(ifc, ifc.x(20)); catch, end; [~, w2] = lastwarn;
lastwarn('', ''); try, sa_deval(ifc, ifc.x(20)); catch, end; [~, w3] = lastwarn;
okw = strcmp(w0, 'MATLAB:deval:MissingSolverField') && strcmp(w1, w0) && strcmp(w3, w2) && ...
    strcmp(w2, 'MATLAB:deval:NonuniqueSolution');
rows = sa_check(rows, T, sprintf(['fallback cases (%d, other solver, bvp5c, no solver field, swapped order, empty ' ...
    'xint, idx 0 / above n / negative, interface, decreasing x, column or single sol.x, idxNonNegative, single, sparse, ' ...
    'NaN, out-of-range, matrix xint, logical idx): same output or same error identifier as deval, through the toolbox; ' ...
    'the missing-solver warning (%s) and the interface warning (%s) reproduced only through the fallback'], ...
    size(fb, 1), w0, w2), okf && okw, strjoin(det, ' | '));
% ---- fast path: double, no warning ------------------------------------------------------------------
lastwarn('sentinel message', 'SA:T62:sentinel');
[a, b] = sa_deval(s, [x mid fliplr(mid)], [1 4]);
[wm, wi] = lastwarn;
rows = sa_check(rows, T, 'the fast path returns double and leaves lastwarn untouched (a planted sentinel survives)', ...
    isa(a, 'double') && isa(b, 'double') && strcmp(wi, 'SA:T62:sentinel') && strcmp(wm, 'sentinel message'), ...
    sprintf('%s / %s; lastwarn %s', class(a), class(b), wi));
end

% =========================================================================

function s = synth(N, n, dom)
x = sort(dom(1) + (dom(2) - dom(1)) * rand(1, N - 2));
x = [dom(1) x dom(2)];
s = struct('solver', 'bvp4c', 'x', x, 'y', randn(n, N), 'yp', randn(n, N), ...
    'stats', struct('maxres', 1e-9, 'nmeshpoints', N));
end

function [eq, fast, why] = compare(args, nout)
% toolbox deval vs sa_deval: same outputs (uint64, class, size) or the same
% error identifier; fast = sa_deval did not reach the toolbox deval
why = '';
[A, ea] = call(@deval, args, nout);
profile('clear'); profile('on');
[B, eb] = call(@sa_deval, args, nout);
profile('off');
pi_ = profile('info');
names = {pi_.FunctionTable.FileName};
fast = ~any(endsWith(names, fullfile('funfun', 'deval.m')));
if ~isempty(ea) || ~isempty(eb)
    eq = strcmp(ea, eb);
    if ~eq, why = sprintf('errors %s | %s', ea, eb); end
    return
end
eq = numel(A) == numel(B);
for k = 1:min(numel(A), numel(B))
    eq = eq && same_bits(A{k}, B{k});
end
if ~eq, why = 'outputs differ'; end
end

function [out, e] = call(f, args, nout)
out = cell(1, nout); e = '';
try
    [out{1:nout}] = f(args{:});
catch ME
    e = ME.identifier; if isempty(e), e = ['<no id> ' ME.message]; end
end
end

function tf = same_bits(a, b)
tf = strcmp(class(a), class(b)) && isequal(size(a), size(b)) && issparse(a) == issparse(b);
if ~tf, return; end
if isa(a, 'double')
    tf = isequal(typecast(full(a(:)), 'uint64'), typecast(full(b(:)), 'uint64'));
elseif isa(a, 'single')
    tf = isequal(typecast(a(:), 'uint32'), typecast(b(:), 'uint32'));
else
    tf = isequaln(a, b);
end
end

% =========================================================================
% the stored-record sweep

function X = sweep(spec)
X = struct('n_records', 0, 'n_evals', 0, 'n_mismatch', 0, 'mismatches', {{}}, 't_toolbox', 0, 't_sa', 0, ...
    'per_run', struct('name', {}, 'records', {}, 'evals', {}, 'mismatch', {}, 't_toolbox', {}, 't_sa', {}), ...
    'lib', struct('name', {}, 'structs', {}, 'evals', {}, 'mismatch', {}));
fid = -1;
if isfield(spec, 'out') && ~isempty(spec.out)
    fid = fopen(spec.out, 'w');
    fprintf(fid, 'run,file,mesh,evals,mismatch,t_toolbox_s,t_sa_s\n');
end
rng(620, 'twister');
for r = 1:numel(spec.runs)
    d = spec.runs(r).dir;
    Ef = load(fullfile(sa_layout('state', d), 'settings_effective.mat'), 'settings', 'model', 'store_settings', 'theta');
    S = sa_settings('legacy_fill', Ef.settings);
    M = sa_model(Ef.model.id);
    frows = [M.fields.row]; irows = [M.identity_fields.row];
    cm = Ef.store_settings.comparison_mesh;
    mo = S.model_options;
    f = dir(fullfile(sa_layout('state', d), 'attempts', 'A*.mat'));
    names = sort({f.name});
    names = names(~cellfun(@isempty, regexp(names, '^A\d{6}\.mat$', 'once')));
    names = names(1:min(numel(names), spec.runs(r).limit));
    pr = struct('name', spec.runs(r).name, 'records', 0, 'evals', 0, 'mismatch', 0, 't_toolbox', 0, 't_sa', 0);
    prev = [];
    for k = 1:numel(names)
        Z = load(fullfile(sa_layout('state', d), 'attempts', names{k}), 'sol');
        sol = Z.sol;
        if isempty(sol), continue; end
        sets = point_sets(sol, S, mo, frows, irows, cm, prev);
        prev = sol;
        ne = 0; nm = 0; t1 = 0; t2 = 0;
        for q = 1:size(sets, 1)
            args = [{sol}, sets{q, 2}];
            nout = sets{q, 3};
            tic; [A, ea] = call(@deval, args, nout); t1 = t1 + toc;
            tic; [B, eb] = call(@sa_deval, args, nout); t2 = t2 + toc;
            ok = strcmp(ea, eb);
            if ok && isempty(ea)
                for j = 1:nout, ok = ok && same_bits(A{j}, B{j}); end
            end
            ne = ne + 1;
            if ~ok
                nm = nm + 1;
                X.mismatches{end+1} = sprintf('%s %s set %s', spec.runs(r).name, names{k}, sets{q, 1});
            end
        end
        pr.records = pr.records + 1; pr.evals = pr.evals + ne; pr.mismatch = pr.mismatch + nm;
        pr.t_toolbox = pr.t_toolbox + t1; pr.t_sa = pr.t_sa + t2;
        if fid > 0
            fprintf(fid, '%s,%s,%d,%d,%d,%.4f,%.4f\n', spec.runs(r).name, names{k}, numel(sol.x), ne, nm, t1, t2);
        end
    end
    X.per_run(end+1) = pr;
    X.n_records = X.n_records + pr.records; X.n_evals = X.n_evals + pr.evals; X.n_mismatch = X.n_mismatch + pr.mismatch;
    X.t_toolbox = X.t_toolbox + pr.t_toolbox; X.t_sa = X.t_sa + pr.t_sa;
    fprintf('T62 sweep %s: %d records, %d evaluations, %d differing; toolbox %.1f s, sa_deval %.1f s\n', ...
        pr.name, pr.records, pr.evals, pr.mismatch, pr.t_toolbox, pr.t_sa);
end
if fid > 0, fclose(fid); end
if isfield(spec, 'libraries')
    for k = 1:numel(spec.libraries)
        X.lib(end+1) = lib_sweep(spec.libraries(k));
        X.n_mismatch = X.n_mismatch + X.lib(end).mismatch;
        X.n_evals = X.n_evals + X.lib(end).evals;
        fprintf('T62 sweep library %s: %d bvp4c structs, %d evaluations, %d differing\n', X.lib(end).name, ...
            X.lib(end).structs, X.lib(end).evals, X.lib(end).mismatch);
    end
end
end

function sets = point_sets(sol, S, mo, frows, irows, cm, prev)
x = sol.x(:).';
h = diff(x);
fr = S.validity.interval_samples(:);
inner = x(1:end-1) + fr .* h;
xd = linspace(x(1), x(end), S.validity.n_dense);
pts = unique([x, inner(:).', xd]);
xm = (x(1:end-1) + x(2:end)) / 2;
sets = {
    'positivity points (field rows)', {pts, frows}, 1
    'n_dense grid with slopes', {xd}, 2
    'n_dense grid, field rows', {xd, frows}, 1
    'n_dense grid, identity rows (sa_identity dense)', {xd, irows}, 1
    'interval midpoints with slopes', {xm}, 2
    'unsorted permutation with slopes', {xd(randperm(numel(xd), min(2001, numel(xd))))}, 2
    };
if numel(cm) >= 2 && cm(1) >= x(1) && cm(end) <= x(end)
    sets(end+1, :) = {'comparison mesh, identity rows', {cm, irows}, 1};
end
if isfield(mo, 'shape_nodes')
    sets(end+1, :) = {'shape-label grid', {linspace(x(1), x(end), mo.shape_nodes)}, 1};
end
eta = linspace(0, 1, 601); cl = @(z) min(max(z, x(1)), x(end));
sets(end+1, :) = {'clamped reconstruction grids', {[cl(eta(eta <= 0.5)), cl(1 - eta(eta > 0.5))]}, 1};
for N = [401 482]
    n = numel(x);
    if n <= N, idx = 1:n; else, idx = round(linspace(1, n, N)); end
    sets(end+1, :) = {sprintf('AM-8 subsample N = %d', N), {x(idx)}, 1}; %#ok<AGROW>
    if ~isempty(prev) && abs(prev.x(1) - x(1)) == 0 && abs(prev.x(end) - x(end)) == 0
        np = numel(prev.x);
        if np <= N, ip = 1:np; else, ip = round(linspace(1, np, N)); end
        sets(end+1, :) = {sprintf('secant: previous record''s subsample N = %d', N), {prev.x(ip)}, 1}; %#ok<AGROW>
    end
end
end

function L = lib_sweep(spec)
L = struct('name', spec.name, 'structs', 0, 'evals', 0, 'mismatch', 0);
M = sa_model(spec.model);
theta = M.normalize(spec.params);
lib = M.library('load', spec.manifest, theta, spec.params);
lib = M.library('profiles', lib, M.library_types);
mo = M.model_options_defaults;
gn = 241; if isfield(mo, 'guess_nodes'), gn = mo.guess_nodes; end
sols = collect(lib, 0, {});
for k = 1:numel(sols)
    s = sols{k};
    a = s.x(1); b = s.x(end);
    cl = @(z) min(max(z, a), b);
    sets = {};
    if isfield(theta, 'Lh')
        xg = linspace(0, theta.Lh, gn);
        sets{end+1} = cl(xg); %#ok<AGROW>
        sets{end+1} = cl(abs(xg - theta.L / 4)); %#ok<AGROW>
    end
    if isfield(theta, 'L')
        sets{end+1} = cl(linspace(0, 0.5, gn) * theta.L); %#ok<AGROW>
    end
    sets{end+1} = cl(linspace(a, b, gn)); %#ok<AGROW>
    for q = 1:numel(sets)
        [A, ea] = call(@deval, {s, sets{q}}, 1);
        [B, eb] = call(@sa_deval, {s, sets{q}}, 1);
        ok = strcmp(ea, eb) && (~isempty(ea) || same_bits(A{1}, B{1}));
        L.evals = L.evals + 1; L.mismatch = L.mismatch + ~ok;
    end
    L.structs = L.structs + 1;
end
end

function sols = collect(v, depth, sols)
% every bvp4c-like struct (fields solver, x, y, yp) inside a loaded library
if depth > 6, return; end
if isstruct(v)
    if isscalar(v) && all(isfield(v, {'x', 'y', 'yp', 'solver'}))
        sols{end+1} = v; return
    end
    for e = 1:numel(v)
        f = fieldnames(v);
        for q = 1:numel(f), sols = collect(v(e).(f{q}), depth + 1, sols); end
    end
elseif iscell(v)
    for e = 1:numel(v), sols = collect(v{e}, depth + 1, sols); end
end
end
