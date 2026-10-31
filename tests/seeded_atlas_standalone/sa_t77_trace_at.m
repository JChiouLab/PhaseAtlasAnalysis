function rows = sa_t77_trace_at(E)
%SA_T77_TRACE_AT  Milestone 7c item 2: a trace started from a figure position
%   (entry.at), the nearest root on the figure, its staleness rule, and the
%   figure's data tips. Zero solves: stub runs, a stub integrator.
T = 'T77 trace by position';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d0 = sa_m3_testdir(E, 't77');
d = fullfile(d0, 'stub');
sa_m3_run(d, [1 1]);
cnt = containers.Map({'n'}, {0});
o = struct('integrator', @(theta, Y0, entry) stub_integrate(theta, Y0, entry, cnt));
po = seeded_atlas_plot(d);
csv = read_plot_csv(fullfile(d, 'plots', 'roots_3d.csv'));
Q = seeded_atlas_plot(d, struct('coordinates_only', true));
rows = sa_check(rows, T, ['the coordinates the trace uses are the plot''s own (coordinates_only equals the plot ' ...
    'output and roots_3d.csv bit for bit)'], isequal(Q.coord, po.coord) && isequal(Q.root_id, po.root_id) && ...
    isequal(Q.coord, csv.C) && isequal(Q.root_id, csv.id), sprintf('%d roots', numel(Q.root_id)));
rng = max(Q.coord, [], 1) - min(Q.coord, [], 1);
Pst = sa_post_model('stub_two_pool'); asp = Pst.aspect(Q.coord); ext = rng ./ asp;
rng2 = max(rng, 0.01 * max(ext) * asp);
rows = sa_check(rows, T, ['the axis ranges are the plotted ranges (max - min over the roots), an axis below 1 % ' ...
    'of the largest extent on the figure raised to that 1 %'], isequal(Q.range, rng2), mat2str(Q.range, 4));

% ---- exact and nearby positions ------------------------------------------------------------------------
j = ceil(numel(Q.root_id) / 2);
eX = struct('name', 'X', 'at', csv.C(j, :));
out = evalc('aX = seeded_atlas_trace(d, eX, o);');
rows = sa_check(rows, T, 'exact plotted coordinates pick that root (distance 0); the root, coordinates and distance are printed', ...
    aX.ran && strcmp(aX.at_resolution.root_id, Q.root_id{j}) && aX.at_resolution.scaled_distance == 0 && ...
    contains(out, Q.root_id{j}) && contains(out, 'scaled') && contains(out, 'raw per axis'), strtrim(out));
ok_all = true; worst = '';
for k = 1:numel(Q.root_id)
    p = Q.coord(k, :) + [0.013 -0.007 0.011] .* Q.range;
    [~, kb] = min(sqrt(sum(((Q.coord - p) ./ Q.range) .^ 2, 2)));          % independent nearest
    evalc('ak = seeded_atlas_trace(d, struct(''name'', sprintf(''N%d'', k), ''at'', p), setfield(o, ''plot'', false));');
    if ~strcmp(ak.at_resolution.root_id, Q.root_id{kb}), ok_all = false; worst = sprintf('k %d', k); end
    raw = Q.coord(kb, :) - p;
    if max(abs(ak.at_resolution.raw_difference - raw)) > 0, ok_all = false; worst = sprintf('raw k %d', k); end
end
rows = sa_check(rows, T, ['a point near each root picks the nearest root on the figure (an independent scaled ' ...
    'distance agrees for every root); raw per-axis differences reported'], ok_all, worst);
Zx = load(fullfile(d, 'traces', 'X', 'trace.mat'));
sx = fileread(fullfile(d, 'traces', 'X', 'trace_summary.txt'));
rows = sa_check(rows, T, 'trace.mat keeps at_resolution and the resolved root id in its key; trace_summary.txt names the position', ...
    strcmp(Zx.key.root_id, Q.root_id{j}) && isequal(Zx.key.entry.at, csv.C(j, :)) && isempty(Zx.key.entry.root) && ...
    strcmp(Zx.at_resolution.root_id, Q.root_id{j}) && contains(sx, 'started from at'));
n0 = cnt('n');
aX2 = seeded_atlas_trace(d, eX, o);
rows = sa_check(rows, T, 'an unchanged at entry whose nearest root is unchanged is reused (zero integration)', ...
    ~aX2.ran && cnt('n') == n0, aX2.reason);

% ---- entry validation; root entries still work ------------------------------------------------------------------
ids = cell(1, 4);
bad = {struct('name', 'B1', 'root', Q.root_id{1}, 'at', [0 0 0]), struct('name', 'B2'), ...
    struct('name', 'B3', 'at', [1 2]), struct('name', 'B4', 'at', [1 NaN 2])};
for k = 1:4
    try, seeded_atlas_trace(d, bad{k}, o); catch ME, ids{k} = ME.identifier; end
end
aR = seeded_atlas_trace(d, struct('name', 'R', 'root', Q.root_id{1}), o);
rows = sa_check(rows, T, 'root and at together, neither, a 2-vector and a NaN are refused (E_ENTRY); a root entry still works', ...
    all(strcmp(ids, 'SA:TRACE:E_ENTRY')) && aR.ran && isempty(fieldnames(aR.at_resolution)), strjoin(ids, ','));

% ---- staleness after a resume that adds roots ------------------------------------------------------------------------
db = fullfile(d0, 'budget');
sa_m3_run(db, [1 1], 'settings', struct('budget', struct('max_calls', 6)));
Qb = seeded_atlas_plot(db, struct('coordinates_only', true));
new = find(~ismember(Q.root_id, Qb.root_id), 1);    % the completed run d holds the roots the resume will add
old = find(ismember(Q.root_id, Qb.root_id), 1);
eN = struct('name', 'N', 'at', Q.coord(new, :));
eO = struct('name', 'O', 'at', Q.coord(old, :));
pN = seeded_atlas_trace(db, eN, o); pO = seeded_atlas_trace(db, eO, o);
sa_m3_run(db, [1 1], 'settings', struct('budget', struct('max_calls', 5000), 'resume', true));
Qa = seeded_atlas_plot(db, struct('coordinates_only', true));
n1 = cnt('n');
aN = seeded_atlas_trace(db, eN, o);
aO = seeded_atlas_trace(db, eO, o);
rows = sa_check(rows, T, ['after a resume adds roots: an at entry whose nearest root changed is integrated again with ' ...
    'the old and the new root in the reason; one whose nearest root is unchanged is reused'], ...
    isequal(Qa.coord, Q.coord) && isequal(Qa.root_id, Q.root_id) && ~isempty(new) && ~strcmp(pN.at_resolution.root_id, Qa.root_id{new}) && aN.ran && ...
    contains(aN.reason, 'nearest root') && contains(aN.reason, pN.at_resolution.root_id) && ...
    contains(aN.reason, Qa.root_id{new}) && strcmp(pO.at_resolution.root_id, Qa.root_id{old}) && ~aO.ran && ...
    cnt('n') == n1 + 1, sprintf('%s | %s', aN.reason, aO.reason));

% ---- data tips in the figure --------------------------------------------------------------------------------------------
f = openfig(fullfile(d, 'plots', 'roots_3d.fig'), 'invisible');
sc = findobj(f, 'Type', 'scatter');
okt = ~isempty(sc); picked = true;
for k = 1:numel(sc)
    rr = sc(k).DataTipTemplate.DataTipRows;
    okt = okt && numel(rr) == 2 && strcmp(rr(1).Label, 'at') && strcmp(rr(2).Label, 'root') && ...
        numel(rr(1).Value) == numel(sc(k).XData);
    for q = 1:numel(rr(1).Value)
        at = str2num(char(rr(1).Value(q))); %#ok<ST2NM>
        a = seeded_atlas_trace(d, struct('name', 'tip', 'at', at), setfield(o, 'plot', false)); %#ok<SFLD>
        picked = picked && strcmp(a.at_resolution.root_id, char(rr(2).Value(q))) && ...
            max(abs(at - [sc(k).XData(q) sc(k).YData(q) sc(k).ZData(q)]) ./ Q.range) < 1e-8;
    end
end
close(f);
fb = dir(fullfile(d, 'plots', 'roots_3d.fig'));
rows = sa_check(rows, T, ['every scatter carries the data tip rows at ([x y z]) and root; pasting a tip''s at picks ' ...
    'that tip''s root; the figure stays small'], okt && picked && fb.bytes < 100e3, sprintf('%d bytes', fb.bytes));
end

% =========================================================================

function C = read_plot_csv(f)
L = strsplit(strtrim(fileread(f)), newline);
L = L(3:end);
C = struct('id', {cell(1, numel(L))}, 'C', zeros(numel(L), 3));
for k = 1:numel(L)
    p = regexp(L{k}, '^"([^"]*)","[^"]*",[^,]*,([^,]*),([^,]*),([^,]*),', 'tokens', 'once');
    C.id{k} = p{1}; C.C(k, :) = str2double(p(2:4));
end
end

function fr = stub_integrate(theta, Y0, entry, cnt)
cnt('n') = cnt('n') + 1;
P = sa_post_model('stub_two_pool');
t = entry.t_out; F = numel(t);
x = linspace(0, 2 * theta.Lh, size(Y0, 2));
O = P.observables(theta, x, Y0);
fr = struct('t', t, 'coord', repmat(O.c, F, 1), 'masses', repmat(O.masses, F, 1), ...
    'umax', repmat(O.umax, F, 1), 'rhs_inf', nan(1, F), 'profile_frames', [1 F], 'stop', 'horizon_reached', ...
    'stop_detail', '', 'wall_s', 0, 'x', x);
fr.profiles = repmat(Y0, 1, 1, 2);
fr.mass_drift_rel = zeros(F, size(O.masses, 2));
fr.Y_final = Y0;
end
