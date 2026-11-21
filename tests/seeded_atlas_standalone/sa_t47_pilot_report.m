function rows = sa_t47_pilot_report(E)
%SA_T47_PILOT_REPORT  The pilot report
%   is a pure function of the run folder (rebuilt twice and from a copy: the
%   same text and tables, the live free-space reading masked); written at
%   every stop reason (completed, budget_calls, budget_wall, disk_floor);
%   a budget-stopped pilot resumes to the uninterrupted run; its sums equal
%   the attempt rows; the projection arithmetic; a planted repeated guess
%   fingerprint is counted. Stub model, default pilot mode (pilot_sample
%   targets, the cross region).
%   Also: f over the whole planned region; the
%   partial report after a budget stop (PARTIAL, pending nodes unreached, b
%   bounded while no refused node is measured); the process/corrector factor
%   on live calls with replay removed at report.replay_s_per_record; section
%   1 shows the current segment's budget; one row per guard (hit counts
%   summed); all 11 tables and a NaN cap-hours at 0 calls; the replay rate
%   is a model-layer leaf (fission yeast 3.6, noxa 0.65) used in the
%   projection.
%   Multiplicity: MULTIPLICITY NOT
%   TESTED in the pilot report (section 7), on section 8's central_nine
%   line and in summary.txt whenever targets.policy is not all_seedable,
%   absent for all_seedable; the repeated-call count excludes the retry of
%   an interrupted record (a stub pilot killed mid-call, then killed again
%   in the retry, resumed: 0 repeats, 2 retries, equivalent to the
%   uninterrupted run) while a planted genuine repeat with retry_of set
%   still counts; noxa report.replay_s_per_record 0.05, and the model-layer
%   digest (model_version) excludes report.*.
T = 'T47 pilot report';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't47');
fail = @(dd, b, ex) sa_m3_wall(dd, ex);
base = struct('mode', 'pilot', 'targets', struct('policy', 'pilot_sample', 'list', zeros(0, 2)), ...
    'continuation', struct('region', 'cross'));
U = fullfile(d, 'U');
rU = run(U, fail, base);
txt = fileread(fullfile(U, 'pilot_report.txt'));
secs = all(cellfun(@(s) contains(txt, s), {'1. SETTINGS', '2. VALIDITY', '3. COST', '4. RETURNED MESH', ...
    '5. ROUTE YIELD', '6. REACH MAP', '7. IDENTITY', '8. PRODUCTION PROJECTION'}));
csvs = {'settings', 'validity_by_check', 'no_return_by_guard', 'cost', 'mesh_by_depth', 'route_yield', ...
    'reach_map', 'frontier', 'identity', 'large_change', 'projection'};
hasc = all(cellfun(@(c) exist(fullfile(sa_layout('state', U), 'pilot', [c '.csv']), 'file') == 2, csvs));
Q = rU.Q;
rows = sa_check(rows, T, 'pilot run (pilot_sample targets, cross region): report with sections 1-8 and the 11 tables', ...
    strcmp(rU.stop_reason, 'completed') && secs && hasc && sum(Q.box) == 25 && ...
    all(Q.index(Q.box, 1) == 0 | Q.index(Q.box, 2) == 0), sprintf('%d calls; region %d nodes', rU.calls, sum(Q.box)));
% ---- pure function of the folder ------------------------------------------------
b1 = snapshot(U);
sa_report('pilot', U);
b2 = snapshot(U);
C2 = fullfile(d, 'copy'); copyfile(U, C2);
sa_report('pilot', C2);
b3 = snapshot(C2);
b3.text = strrep(b3.text, C2, U);
rows = sa_check(rows, T, 'pure function: rebuilt in place and from a copy, identical text (live free-space line masked) and tables', ...
    isequal(b1, b2) && isequal(b1, b3), '');
% ---- every stop reason ------------------------------------------------------------
B = fullfile(d, 'budget');
rB = run(B, fail, setf(base, 'budget', struct('max_calls', 40)));
W = fullfile(d, 'wall');
rW = run(W, fail, setf(base, 'budget', struct('max_wall_s', 1e-6)));
F = fullfile(d, 'disk');
rF = run(F, fail, setf(base, 'checkpoint', struct('every_calls', 5, 'min_free_gib', 1e9)));
ex = @(x) exist(fullfile(x, 'pilot_report.txt'), 'file') == 2;
% ---- the partial report of a pilot stopped mid-cross ----------------------
% (B above stops inside the cold pass, before any scheduler state; B2 stops
% 4 calls into the continuation, as the S5 control did)
k1 = find(~ismember({rU.attempts.kind}, {'cold', 'refinement'}), 1);
B2 = fullfile(d, 'budget_midcross');
rB2 = run(B2, fail, setf(base, 'budget', struct('max_calls', k1 + 3)));
RB = sa_report('pilot', B2);
PB = RB.projection; tB = fileread(fullfile(B2, 'pilot_report.txt'));
% the projection's cold count is PB.cold (the per_seed seed count, the production
% default), not PB.cold_all (all_seedable)
cB = PB.cold + PB.N * PB.f * PB.a + PB.N * (1 - PB.f) * PB.b;
pend = RB.reach.classes.unattempted + RB.reach.classes.policy_skipped + RB.reach.classes.unfinished;
okb = (PB.n_refused_measured > 0 && strcmp(PB.b_source, 'measured')) || ...
    (PB.n_refused_measured == 0 && PB.b == PB.b_bound && startsWith(PB.b_source, 'bound'));
rows = sa_check(rows, T, 'partial report: PARTIAL printed, f = reached / region with the pending nodes unreached, b measured or bounded, arithmetic', ...
    strcmp(rB2.stop_reason, 'budget_calls') && PB.partial && contains(tB, 'PARTIAL: stop reason budget_calls') && pend > 0 && ...
    abs(PB.f - RB.reach.classes.reached / RB.reach.n_region) < 1e-12 && PB.unreached == RB.reach.n_region - ...
    RB.reach.classes.reached && okb && abs(PB.calls(1) - cB) < 1e-9 * cB, ...
    sprintf('reached %d / region %d, pending %d, refused measured %d, b %.3g (%s), central %.0f', ...
    RB.reach.classes.reached, RB.reach.n_region, pend, PB.n_refused_measured, PB.b, PB.b_source, PB.calls(1)));
rows = sa_check(rows, T, 'written at every stop: completed, budget_calls (mid-run), budget_wall, disk_floor (0 calls)', ...
    ex(U) && ex(B) && ex(W) && ex(F) && strcmp(rB.stop_reason, 'budget_calls') && rB.calls == 40 && ...
    strcmp(rW.stop_reason, 'budget_wall') && strcmp(rF.stop_reason, 'disk_floor') && rF.calls == 0, ...
    sprintf('%s %d / %s %d / %s %d', rB.stop_reason, rB.calls, rW.stop_reason, rW.calls, rF.stop_reason, rF.calls));
rR = run(B, fail, setf(setf(base, 'budget', struct('max_calls', 5000)), 'resume', true));
D = sa_m3_equiv(U, B);
tb = fileread(fullfile(B, 'pilot_report.txt'));
rows = sa_check(rows, T, 'the budget-stopped pilot resumes to the uninterrupted run; the resumed report shows the replay cost', ...
    D.pass && strcmp(rR.stop_reason, 'completed') && contains(tb, 'replay s per record'), strjoin(D.diffs, ' | '));
% ---- live overhead after the resume; section 1 = this segment ---------
RR = sa_report('pilot', B);
ZA = load(fullfile(sa_layout('state', B), 'attempts.mat')); ZP = load(fullfile(sa_layout('state', B), 'process_time.mat'));
wB = [ZA.attempts.wall_s]; wl = sum(wB(41:end), 'omitnan');
rem = min(40 * 0.65, max(0, ZP.process_s - wl));
rows = sa_check(rows, T, 'resumed: replay removed = 40 records x 0.65 s (package rate, stub), factor = (process - replay) / live corrector s', ...
    abs(RR.cost.replay_removed_s - rem) < 1e-12 && abs(RR.cost.overhead - (ZP.process_s - rem) / wl) < 1e-9 && ...
    RR.projection.replay_rate == 0.65, sprintf('process %.2f s, replay removed %.2f s, live corrector %.3f s, factor %.3g', ...
    ZP.process_s, RR.cost.replay_removed_s, wl, RR.cost.overhead));
rows = sa_check(rows, T, 'section 1 of the resumed run: this segment''s budget (5000), the first segment''s (40), leaf sourced resume_1', ...
    contains(tb, 'budget 5000 this segment; segment 2 after resume 1; first segment''s budget 40') && ...
    ~isempty(regexp(tb, 'budget\.max_calls\s+5000\s+resume_1', 'once')), '');
% ---- every table at 0 calls; cap hours not computed ------------------------
pj = rd(fullfile(sa_layout('state', F), 'pilot', 'projection.csv'));
hasF = all(cellfun(@(c) exist(fullfile(sa_layout('state', F), 'pilot', [c '.csv']), 'file') == 2, csvs));
rows = sa_check(rows, T, 'disk_floor at 0 calls: all 11 tables written; MATLAB-h cap NaN (not 0)', ...
    hasF && rF.calls == 0 && strcmp(pj{strcmp(pj(:, 1), 'matlab_h'), 4}, 'NaN'), ...
    sprintf('matlab_h cap %s', pj{strcmp(pj(:, 1), 'matlab_h'), 4}));
% ---- sums equal the attempt rows --------------------------------------------------------
A = rU.attempts;
R = sa_report('pilot', U);
cost = rd(fullfile(sa_layout('state', U), 'pilot', 'cost.csv'));
cl = cost(~startsWith(cost(:, 1), 'nmax_share'), :);
vb = rd(fullfile(sa_layout('state', U), 'pilot', 'validity_by_check.csv'));
nf = 0;
for k = 1:numel(A)
    if ~any(strcmp(A(k).outcome, {'no_return', 'interrupted'})) && ~isempty(A(k).failing)
        nf = nf + numel(strsplit(A(k).failing, ';'));
    end
end
nfirst = sum(~cellfun(@isempty, {A(~strcmp({A.outcome}, 'no_return')).failing}));
ry = rd(fullfile(sa_layout('state', U), 'pilot', 'route_yield.csv'));
ec = sum(str2double(ry(startsWith(ry(:, 1), 'edges '), 3)));
rows = sa_check(rows, T, 'sums: cost classes = attempt rows; all-failing = failing entries; first-failing = failing returns; completed-edge calls', ...
    sum(str2double(cl(:, 2))) == numel(A) && sum(str2double(vb(:, 5))) == nf && sum(str2double(vb(:, 4))) == nfirst ...
    && ec == sum([Q.edges(strcmp({Q.edges.status}, 'completed')).calls]), ...
    sprintf('rows %d; failing %d; first %d', numel(A), nf, nfirst));
P = R.projection;
% cold = P.cold, the per_seed seed count (production default), <= P.cold_all
c = P.cold + P.N * P.f * P.a + P.N * (1 - P.f) * P.b;
h = P.cold + P.N * P.a + P.N * (1 - P.f) * P.b;
rows = sa_check(rows, T, 'projection arithmetic: central = cold + N f a + N (1 - f) b, high = cold + N a + N (1 - f) b, cap >= high', ...
    abs(P.calls(1) - c) < 1e-9 * c && abs(P.calls(2) - h) < 1e-9 * h && P.calls(3) >= P.calls(2) && ...
    P.N == 169 && P.cold_all > 0 && strcmp(P.cold_policy, 'per_seed') && P.cold > 0 && P.cold <= P.cold_all, ...
    sprintf('cold %d (all_seedable %d) N %d f %.3g a %.3g b %.3g -> %s', P.cold, P.cold_all, P.N, P.f, P.a, ...
    P.b, mat2str(round(P.calls))));
rc = R.reach.classes;
rows = sa_check(rows, T, 'f = reached / ALL region nodes (refused and unattempted nodes unreached), not reached / attempted', ...
    ~P.partial && abs(P.f - rc.reached / R.reach.n_region) < 1e-12 && R.reach.n_region == 25 && ...
    P.unreached == rc.refused + rc.unattempted + rc.policy_skipped + rc.unfinished, ...
    sprintf('reached %d refused %d unattempted %d of %d: f %.3g (attempted-only would be %.3g)', rc.reached, ...
    rc.refused, rc.unattempted, R.reach.n_region, P.f, rc.reached / max(1, rc.reached + rc.refused)));
% ---- the replay rate is a model-layer leaf, used in the projection ---------
[Sf, Pf] = sa_settings('resolve', struct('model', 'original_polarity', 'run_dir', 'x'), sa_model('original_polarity'));
[Sn, Pn] = sa_settings('resolve', struct('model', 'guzzo_noxa', 'run_dir', 'x'), sa_model('guzzo_noxa'));
Sx = load(fullfile(sa_layout('state', U), 'settings_effective.mat'));
rows = sa_check(rows, T, 'report.replay_s_per_record: model layer fission yeast 0.35, noxa 0.05; stub 0.65 from the package', ...
    Sf.report.replay_s_per_record == 0.35 && strcmp(Pf.report.replay_s_per_record, 'model') && ...
    Sn.report.replay_s_per_record == 0.05 && strcmp(Pn.report.replay_s_per_record, 'model') && ...
    Sx.settings.report.replay_s_per_record == 0.65 && strcmp(Sx.settings_sources.report.replay_s_per_record, 'package'), '');
Wd = fullfile(d, 'wallseg'); copyfile(U, Wd);
Zw = load(fullfile(sa_layout('state', Wd), 'settings_effective.mat'));
Zw.settings.budget.wall_s = 1e-6; Zw.settings.report.replay_s_per_record = 2;
sa_write('mat', fullfile(sa_layout('state', Wd), 'settings_effective.mat'), Zw);
Pw = sa_report('pilot', Wd); Pw = Pw.projection;
nres = floor(Pw.hours_live(1) * 3600 / 1e-6);
fin = isfinite(Pw.calls);                         % the cap column may be unbounded (Inf)
rows = sa_check(rows, T, 'projection: replay term = resumes (budget.wall_s) x rate x calls / 2 at the planted rate 2 s', ...
    nres > 0 && Pw.resumes == nres && all(abs(Pw.hours(fin) - (Pw.hours_live(fin) + nres * 2 * Pw.calls(fin) / 2 / 3600)) <= ...
    1e-9 * Pw.hours(fin)) && Pw.replay_rate == 2 && all(abs(Pw.hours_live(fin) - P.hours_live(fin)) <= 1e-12 * P.hours_live(fin)) && ...
    isequal(~fin, isinf(Pw.calls)), ...   % the cap column is unbounded (Inf) at max_lineages 10000
    sprintf('resumes %d, hours %s', nres, mat2str(Pw.hours, 4)));
% ---- one row per guard and error identifier, hits summed ------------------
Gd = fullfile(d, 'guard'); copyfile(U, Gd);
Zg = load(fullfile(sa_layout('state', Gd), 'attempts.mat'));
Zg.attempts(3).outcome = 'no_return'; Zg.attempts(5).outcome = 'no_return';
Zg.attempts(3).corrector_reason = 'singularity_guard: 12 hit(s); MATLAB:bvp4c:SingJac Unable to solve the collocation equations';
Zg.attempts(5).corrector_reason = 'singularity_guard: 1345 hit(s); MATLAB:bvp4c:SingJac Unable to solve the collocation equations';
sa_write('mat', fullfile(sa_layout('state', Gd), 'attempts.mat'), Zg);
sa_report('pilot', Gd);
gg = rd(fullfile(sa_layout('state', Gd), 'pilot', 'no_return_by_guard.csv'));
ig = find(strcmp(gg(:, 1), 'singularity_guard MATLAB:bvp4c:SingJac'));
rows = sa_check(rows, T, 'no_return by guard: two planted guard calls with different hit counts are one row (calls 2, hits 1357)', ...
    numel(ig) == 1 && str2double(gg{ig, 2}) == 2 && str2double(gg{ig, 3}) == 1357 && ...
    sum(contains(gg(:, 1), 'singularity_guard')) == 1, strjoin(gg(:, 1).', ' | '));
% ---- planted repeat --------------------------------------------------------------------
Pd = fullfile(d, 'planted'); copyfile(U, Pd);
Z = load(fullfile(sa_layout('state', Pd), 'attempts.mat'));
Z.guess_fingerprint{7} = Z.guess_fingerprint{4};
sa_write('mat', fullfile(sa_layout('state', Pd), 'attempts.mat'), Z);
Rp = sa_report('pilot', Pd);
rows = sa_check(rows, T, 'repeated calls: 0 in the run; a planted repeated guess fingerprint counts 1', ...
    R.yield.repeated == 0 && Rp.yield.repeated == 1 && R.yield.fingerprinted == numel(A), ...
    sprintf('%d / %d of %d', R.yield.repeated, Rp.yield.repeated, R.yield.fingerprinted));
% ---- the retry of an interrupted record is not a repeat ------------------
kk = k1 + 2;
K = fullfile(d, 'killed');
kill(K, fail, base, kk);                                   % kill mid-call (continuation)
kill(K, fail, setf(base, 'resume', true), 1);              % the retry itself is killed
rK = run(K, fail, setf(setf(base, 'resume', true), 'budget', struct('max_calls', 5000)));
DK = sa_m3_equiv(U, K);
RK = sa_report('pilot', K);
tK = fileread(fullfile(K, 'pilot_report.txt'));
ry = rd(fullfile(sa_layout('state', K), 'pilot', 'route_yield.csv'));
nint = sum(strcmp({rK.attempts.outcome}, 'interrupted'));
rows = sa_check(rows, T, ['a pilot killed mid-call and again in the retry, resumed: equivalent to the ' ...
    'uninterrupted run; repeated calls 0, retries of interrupted records 2 (reported separately)'], ...
    DK.pass && nint == 2 && RK.yield.repeated == 0 && RK.yield.retries == 2 && ...
    contains(tK, 'repeated calls (equal guess fingerprint to an earlier call; the retry of an interrupted record excluded): 0 of') && ...
    str2double(ry{strcmp(ry(:, 1), 'retries_of_interrupted'), 2}) == 2, ...
    sprintf('interrupted %d; repeated %d; retries %d; %s', nint, RK.yield.repeated, RK.yield.retries, ...
    strjoin(DK.diffs, ' | ')));
Kp = fullfile(d, 'planted_retry'); copyfile(U, Kp);
Z = load(fullfile(sa_layout('state', Kp), 'attempts.mat'));
Z.guess_fingerprint{9} = Z.guess_fingerprint{6}; Z.attempts(9).retry_of = Z.attempts(6).attempt_id;
sa_write('mat', fullfile(sa_layout('state', Kp), 'attempts.mat'), Z);
Rq = sa_report('pilot', Kp);
rows = sa_check(rows, T, 'retry_of set but the earlier equal fingerprint is not an interrupted record: still a repeat', ...
    Rq.yield.repeated == 1 && Rq.yield.retries == 0, sprintf('repeated %d retries %d', Rq.yield.repeated, Rq.yield.retries));
% ---- multiplicity warning unless the cold targets are all_seedable -------
sU = fileread(fullfile(U, 'summary.txt'));
idc = rd(fullfile(sa_layout('state', U), 'pilot', 'identity.csv'));
txt = fileread(fullfile(U, 'pilot_report.txt'));
s7 = regexp(txt, '7\. IDENTITY(.*)8\. PRODUCTION', 'tokens', 'once');
c9l = regexp(txt, '[^\n]*with central_nine cold instead[^\n]*', 'match', 'once');
Ad = fullfile(d, 'all_seedable'); copyfile(U, Ad);
Za = load(fullfile(sa_layout('state', Ad), 'settings_effective.mat')); Za.settings.targets.policy = 'all_seedable';
sa_write('mat', fullfile(sa_layout('state', Ad), 'settings_effective.mat'), Za);
Ra = sa_report('pilot', Ad);
ta = fileread(fullfile(Ad, 'pilot_report.txt'));
ida = rd(fullfile(sa_layout('state', Ad), 'pilot', 'identity.csv'));
na = sa_report('multiplicity_note', struct('targets', struct('policy', 'all_seedable')));
nc = sa_report('multiplicity_note', struct('targets', struct('policy', 'central_nine')));
rows = sa_check(rows, T, ['MULTIPLICITY NOT TESTED: pilot_sample run -> section 7, summary.txt, identity.csv 0; ' ...
    'section 8 central_nine line carries the caveat; all_seedable -> no warning, identity.csv 1'], ...
    ~isempty(s7) && contains(s7{1}, 'MULTIPLICITY NOT TESTED: cold targets = pilot_sample') && ...
    contains(sU, 'MULTIPLICITY NOT TESTED: cold targets = pilot_sample') && ...
    contains(c9l, 'branch connected to the central nine only; multiplicity not tested') && ...
    strcmp(idc{strcmp(idc(:, 1), 'multiplicity_tested'), 2}, '0') && ~R.identity.multiplicity_tested && ...
    ~contains(ta, 'MULTIPLICITY NOT TESTED') && Ra.identity.multiplicity_tested && ...
    strcmp(ida{strcmp(ida(:, 1), 'multiplicity_tested'), 2}, '1') && isempty(na) && ...
    contains(nc, 'connected to the central nine nodes only'), c9l);
% ---- the model-layer digest EXCLUDES report.* -- a report-only recalibration does not
% change model_version ---------------------------------------------------------------
Mn = sa_model('guzzo_noxa');
% the solver: component follows the setting_defaults tag
tag = regexp(Mn.version, 'setting_defaults:([0-9a-f]+)(?= solver:|$)', 'tokens', 'once');
sd = Mn.setting_defaults; sd2 = rmfield(sd, 'report');
sd3 = sd; sd3.report.replay_s_per_record = 0.65;
rows = sa_check(rows, T, ['noxa model_version digest = SHA-256 of the model layer WITHOUT ' ...
    'report.* (changing report.replay_s_per_record leaves model_version unchanged; USER_GUIDE section 7)'], ...
    ~isempty(tag) && strcmp(tag{1}, digest(sd2)) && ~strcmp(tag{1}, digest(sd)) && ...
    strcmp(digest(rmfield(sd3, 'report')), tag{1}), Mn.version(end-28:end));
rmpath(E.stub);
end

function kill(dd, fail, o, k)
try
    sa_m3_run(dd, [], 'fail', fail, 'kill_at', k, 'settings', o);
catch ME
    assert(strcmp(ME.identifier, 'SA:TEST:KILL'), ME.message);
end
end

function h = digest(sd)
% the model factories' setting_defaults_tag, recomputed independently
txt = leaf_text(sd, '');
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(txt), 'int8'));
b = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(b, 2).', 1, []));
h = h(1:12);
end

function t = leaf_text(s, path)
t = '';
f = fieldnames(s);
for k = 1:numel(f)
    p = [path '.' f{k}];
    v = s.(f{k});
    if isstruct(v), t = [t leaf_text(v, p)]; %#ok<AGROW>
    elseif ischar(v), t = [t p '=' v ';']; %#ok<AGROW>
    else, t = [t p '=' mat2str(v, 17) ';']; %#ok<AGROW>
    end
end
end

function r = run(dd, fail, o)
r = sa_m3_run(dd, [], 'fail', fail, 'settings', o);
end

function s = setf(s, f, v)
s.(f) = v;
end

function b = snapshot(d)
t = fileread(fullfile(d, 'pilot_report.txt'));
t = regexprep(t, 'free above the [^\n]*', '<live free space>');
f = dir(fullfile(sa_layout('state', d), 'pilot', '*.csv'));
c = cellfun(@(n) fileread(fullfile(sa_layout('state', d), 'pilot', n)), {f.name}, 'UniformOutput', false);
b = struct('text', t, 'names', {{f.name}}, 'csv', {c});
end

function C = rd(f)
t = strsplit(strtrim(fileread(f)), newline);
t = t(2:end);
C = cell(numel(t), 0);
for k = 1:numel(t)
    p = regexp(t{k}, ',(?=(?:[^"]*"[^"]*")*[^"]*$)', 'split');
    C(k, 1:numel(p)) = p;
end
end
