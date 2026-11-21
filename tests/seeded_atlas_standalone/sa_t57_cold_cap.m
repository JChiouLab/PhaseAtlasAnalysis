function rows = sa_t57_cold_cap(E)
%SA_T57_COLD_CAP  The
%   cold-only mesh cap. (1) solver_for: cold capped (nmax = mesh_reject =
%   cap), refinement and continuation unchanged; cap Inf and a struct without a cap
%   leave S.solver. (2) The cut_by_cap rule on planted returns through the
%   real record path (stub model, cold pass under cap 10,000 with
%   solver.nmax = mesh_reject = 20,000): proposed mesh > nmax; returned mesh
%   >= reject; a coarsened 3,601-point RelTolNotMet return (proposed 10,801);
%   the 9,695 case (proposed > nmax, returned mesh below the cap): cut; the
%   case of a positivity refusal at a returned mesh above the cap:
%   cut, outcome unchanged; a positivity-only refusal and a plain no-return:
%   not cut; an accepted return: not cut; the same calls without the cap: no
%   cut at all. The rule on its own: a cap set through mesh_reject alone; a
%   no-return whose reason is the module's mesh blow-up; continuation calls
%   never cut. (3) Journal round trip: live 'end' lines carry ' cut_by_cap';
%   a kill between the attempt file and its 'end' line -> the reconcile
%   writes the same suffix from the stored row; the replayed rows are
%   isequaln. (4) A changed cap is refused on resume at 0 calls.
T = 'T57 cold cap';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
% ---- (1) solver_for ----------------------------------------------------------------------
Mf = sa_model('original_polarity');
S = sa_settings('resolve', struct('model', 'original_polarity', 'run_dir', 'x'), Mf);
c = sa_search('solver_for', S, 'cold'); r = sa_search('solver_for', S, 'refinement');
k = sa_search('solver_for', S, 'continuation'); q = sa_search('solver_for', S, 'rescue');
S2 = S; S2.cold.mesh_cap = Inf; c2 = sa_search('solver_for', S2, 'cold');
S3 = rmfield(S, 'cold'); c3 = sa_search('solver_for', S3, 'cold');
rows = sa_check(rows, T, ['solver_for: cold nmax = mesh_reject = 10000 (model default), refinement nmax = mesh_reject = ' ...
    '20000 at the refinement tolerances, continuation/rescue = S.solver; cap Inf or no cold group: cold = S.solver'], ...
    c.nmax == 10000 && c.mesh_reject == 10000 && c.reltol == S.solver.reltol && r.nmax == 20000 && ...
    r.mesh_reject == 20000 && r.reltol == S.validity.refine_reltol && isequal(k, S.solver) && isequal(q, S.solver) && ...
    isequal(c2, S.solver) && isequal(c3, S.solver), sprintf('cold %d/%d refinement %d/%d', c.nmax, c.mesh_reject, ...
    r.nmax, r.mesh_reject));
% ---- (2) planted returns through the record path --------------------------------------------
addpath(E.stub);
d = sa_m3_testdir(E, 't57');
plan = {'proposed', 'reject', 'coarse3601', 'c9695', 'a100', 'positivity', 'noreturn', 'ok'};
want_cut = [true true true true true false false false];
want_out = {'invalid', 'invalid', 'invalid', 'invalid', 'invalid', 'invalid', 'no_return', 'accepted'};
tl = [0 0; 1 0; 2 0; 3 0; 0 1; 1 1; 2 1; 3 1];
o = struct('targets', struct('list', tl), 'routes', struct('continuation', false, 'rescue', false), ...
    'cold', struct('mesh_cap', 10000));
R = fullfile(d, 'cap');
res = run_planted(R, o, plan);
A = res.attempts;
ok = numel(A) == numel(plan) && isequal([A.cut_by_cap], want_cut) && isequal({A.outcome}, want_out) && ...
    all([A.nmax] == 10000);
rows = sa_check(rows, T, ['cut_by_cap on planted cold returns (cap 10,000): proposed > nmax, mesh >= reject, coarsened ' ...
    '3,601 RelTolNotMet, the 9,695 case, the A000100 case (positivity at 11,743: cut, still invalid) -> cut; ' ...
    'positivity-only, no-return, accepted -> not cut; row nmax = the call''s own 10,000'], ok, ...
    sprintf('cut %s; outcomes %s; failing(A000005) %s', mat2str([A.cut_by_cap]), strjoin({A.outcome}, ','), A(5).failing));
R0 = fullfile(d, 'nocap'); o0 = o; o0.cold.mesh_cap = Inf;
r0 = run_planted(R0, o0, plan);
A0 = r0.attempts;
rows = sa_check(rows, T, 'the same planted calls without the cap: no call cut by cap; row nmax = 20000', ...
    ~any([A0.cut_by_cap]) && all([A0.nmax] == 20000), sprintf('outcomes %s', strjoin({A0.outcome}, ',')));
sv = struct('nmax', 20000, 'mesh_reject', 20000);
t1 = sa_record_call('cut_by_cap', 'cold', struct('nmax', 20000, 'mesh_reject', 10000), sv, {'mesh'}, '', true);
t2 = sa_record_call('cut_by_cap', 'cold', struct('nmax', 20000, 'mesh_reject', 20000), sv, {'mesh'}, '', true);
t3 = sa_record_call('cut_by_cap', 'cold', struct('nmax', 10000, 'mesh_reject', 10000), sv, ...
    {'convergence_maxres', 'convergence_proposed_mesh', 'mesh'}, 'mesh blow-up (12001 points >= 10000)', false);
t4 = sa_record_call('cut_by_cap', 'cold', struct('nmax', 10000, 'mesh_reject', 10000), sv, ...
    {'convergence_maxres', 'convergence_proposed_mesh', 'mesh'}, 'SingJac', false);
t5 = sa_record_call('cut_by_cap', 'continuation', struct('nmax', 10000, 'mesh_reject', 10000), sv, {'mesh'}, '', true);
t6 = sa_record_call('cut_by_cap', 'cold', struct('nmax', 10000, 'mesh_reject', 10000), sv, 'positivity', '', true);
rows = sa_check(rows, T, ['the rule alone: a cap through mesh_reject alone -> cut; no tighter option -> not; a no-return ' ...
    'with the mesh blow-up reason -> cut; another no-return (every check NaN) -> not; continuation never; ' ...
    'positivity only -> not'], t1 && ~t2 && t3 && ~t4 && ~t5 && ~t6, mat2str([t1 t2 t3 t4 t5 t6]));
% ---- (3) journal round trip ------------------------------------------------------------------
J = strsplit(strtrim(fileread(fullfile(sa_layout('state', R), 'journal.txt'))), newline);
ends = J(startsWith(J, 'end '));
suf = cellfun(@(s) endsWith(s, ' cut_by_cap'), ends);
rows = sa_check(rows, T, 'live journal: exactly the cut calls'' end lines carry the suffix " cut_by_cap"', ...
    isequal(suf, want_cut), strjoin(ends, ' | '));
K = fullfile(d, 'fileend');
copyfile(R, K);
Z = load(fullfile(sa_layout('state', K), 'settings_effective.mat'));
Z.settings.run_dir = K;
sa_write('mat', fullfile(sa_layout('state', K), 'settings_effective.mat'), Z);
L = strsplit(strtrim(fileread(fullfile(sa_layout('state', K), 'journal.txt'))), newline);
L = L(~startsWith(L, 'end A000003') & ~startsWith(L, 'end A000006'));   % a cut and an uncut call
fid = fopen(fullfile(sa_layout('state', K), 'journal.txt'), 'w'); fprintf(fid, '%s\n', L{:}); fclose(fid);
oK = o; oK.resume = true;
rK = run_planted(K, oK, plan);
J2 = strsplit(strtrim(fileread(fullfile(sa_layout('state', K), 'journal.txt'))), newline);
AK = load(fullfile(sa_layout('state', K), 'attempts.mat'));
A1 = load(fullfile(sa_layout('state', R), 'attempts.mat'));
rows = sa_check(rows, T, ['kill between the attempt file and its end line: the reconcile appends "end A000003 invalid ' ...
    'cut_by_cap" and "end A000006 invalid" from the stored rows; the resume replays every row isequaln'], ...
    any(strcmp(J2, 'end A000003 invalid cut_by_cap')) && any(strcmp(J2, 'end A000006 invalid')) && ...
    isempty(rK.replay_interrupted) && isequaln(rmfield(AK.attempts, 'wall_s'), rmfield(A1.attempts, 'wall_s')) && ...
    isequal([rK.attempts.cut_by_cap], want_cut), strjoin(J2(end-3:end), ' | '));
% ---- (4) a changed cap refused on resume --------------------------------------------------------
n0 = numel(dir(fullfile(sa_layout('state', R), 'attempts', 'A*.mat')));
o4 = o; o4.resume = true; o4.cold.mesh_cap = 9000;
g = '';
try, run_planted(R, o4, plan); catch ME, g = ME.identifier; end
rows = sa_check(rows, T, 'a resume with another cold.mesh_cap is refused (SA:RESUME:E_SETTINGS) at 0 calls', ...
    strcmp(g, 'SA:RESUME:E_SETTINGS') && n0 == numel(dir(fullfile(sa_layout('state', R), 'attempts', 'A*.mat'))), g);
rmpath(E.stub);
end

% =========================================================================

function res = run_planted(R, o, plan)
s = sa_m3_stub('settings', R, []);
s = merge(s, o);
s.budget = struct('max_calls', numel(plan));
p = sa_m3_stub('params');
st = containers.Map('KeyType', 'char', 'ValueType', 'any'); st('n') = 0;
res = seeded_atlas_run(p, struct(), s, @(th, g, t, op) planted(th, g, t, op, st, plan));
end

function [sol, info] = planted(th, g, t, op, st, plan) %#ok<INUSL>
n = st('n') + 1; st('n') = n;
info = struct('class', 'STUB', 'warning_id', '');
sol = [];
switch plan{n}
    case 'proposed'      % proposed mesh above nmax under the cap; 201 points returned
        sol = mk(th, t.n, 0.1, 201); sol.stats.nmeshpoints = 12000;
        if sol.stats.nmeshpoints > op.nmax, info.warning_id = 'MATLAB:bvp4c:RelTolNotMet'; sol.stats.maxres = 2e-6; end
    case 'reject'        % a returned mesh at 10,000 points: >= reject only under the cap
        sol = mk(th, t.n, 0.1, 10000);
    case 'coarse3601'    % refine-then-coarsen: 3,601 returned, 10,801 proposed, RelTolNotMet
        sol = mk(th, t.n, 0.1, 3601); sol.stats.nmeshpoints = 10801;
        if 10801 > op.nmax, info.warning_id = 'MATLAB:bvp4c:RelTolNotMet'; sol.stats.maxres = 2e-6; end
    case 'c9695'         % NMax-bound at a last solved mesh below the cap
        sol = mk(th, t.n, 0.1, 9695); sol.stats.nmeshpoints = 10554;
        if 10554 > op.nmax, info.warning_id = 'MATLAB:bvp4c:RelTolNotMet'; sol.stats.maxres = 2e-6; end
    case 'a100'          % a converged positivity refusal returned at 11,743 points (above the cap)
        sol = mk(th, t.n, 1.5, 11743);
    case 'positivity'    % positivity only, small mesh
        sol = mk(th, t.n, 1.5, 201);
    case 'noreturn'
        sol = [];
    case 'ok'
        sol = mk(th, t.n, 0.1, 201);
end
end

function sol = mk(th, n, a, m)
k = pi / th.Lh;
x = linspace(0, th.Lh, m);
y = zeros(8, m); yp = y;
for j = 1:4
    y(2*j-1, :) = n(j) * (1 + a * cos(k * x));
    y(2*j, :) = -n(j) * a * k * sin(k * x);
    yp(2*j-1, :) = y(2*j, :);
    yp(2*j, :) = -n(j) * a * k^2 * cos(k * x);
end
sol = struct('solver', 'bvp4c', 'x', x, 'y', y, 'yp', yp, 'parameters', n(:), ...
    'stats', struct('maxres', 1e-9, 'nmeshpoints', m));
end

function a = merge(a, b)
f = fieldnames(b);
for k = 1:numel(f)
    if isstruct(b.(f{k})) && isfield(a, f{k}) && isstruct(a.(f{k}))
        a.(f{k}) = merge(a.(f{k}), b.(f{k}));
    else
        a.(f{k}) = b.(f{k});
    end
end
end
