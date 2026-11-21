function rows = sa_t102_release_selfcheck(E)
%SA_T102_RELEASE_SELFCHECK  The
%   release self-check on synthetic "other release" toolboxes. Fake toolbox
%   folders are built in scratch from this MATLAB's toolbox files (deleted
%   after; never shipped) and sa_bvp4c_setup('generate', <scratch out>,
%   <fake toolbox>) is called on each. The installed copy is never touched.
%   About 100 tiny bvp4c calls, single thread.
%     1  byte-identical copy: accepted, the label of this toolbox (reference
%        and deval replica on the reference toolbox);
%     2  one changed byte in a comment of bvp4c.m: accepted, self-checked;
%        the solver: component differs, so a (stub) run made under 1 is
%        refused on resume under 2 with the toolbox-change message;
%     3  re-indented anchor lines: accepted through the tolerant path
%        (recorded), inserted lines take the anchor's indentation;
%     4  a duplicated proposal anchor in the main body: E_ANCHOR with counts;
%     5  a removed private helper: E_CLOSURE;
%     6  a numerical change in the residual subfunction: with the fake folder
%        first on the path (so the user's bvp4c is the changed one) accepted
%        (the self-check certifies the patch, not the numerics); without it,
%        refused E_SELFCHECK (the copy differs from the path's bvp4c);
%     7  a changed byte in deval.m: accepted with deval_mode toolbox, and
%        sa_deval in that mode equals the toolbox deval;
%     8  a new_profile call without canRemovePoints (the older signature,
%        the subfunction defaulting it): accepted through the tolerant path
%        (fake folder first on the path);
%     9  a .p sibling of bvp4c.m: E_TOOLBOX;
%    10  a CRLF-converted bvp4c.m: accepted; every line keeps CRLF and the
%        copy equals the LF copy with LF -> CRLF (self-check case g holds).
T = 'T102 release self-check';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
P = sa_bvp4c_setup('pins');
tb = sa_bvp4c_setup('locate');
d = sa_m3_testdir(E, 't102');
sd = fullfile(E.pkg, 'solver');
inst = fullfile(sd, 'sa_bvp4c_traced.m');
inst0 = bytes_or_empty(inst);
R0 = sa_bvp4c_setup('generate', fullfile(d, 'gen_real'));       % this MATLAB's own toolbox (the comparison)
ref = strcmp(R0.toolbox_label, 'reference');
H = R0.helpers;
LF = char(bytes(fullfile(tb, 'bvp4c.m')));
lines = strsplit(LF, newline, 'CollapseDelimiters', false);
i2 = num(R0.anchor_initial); i3 = num(R0.anchor_proposal);
fakes = {};
% ---- 1 byte-identical -------------------------------------------------------------------------
f1 = fake(tb, fullfile(d, 'fake1'), H); fakes{end+1} = f1;
[e1, R1] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen1'), f1));
ok1 = isempty(e1) && strcmp(R1.toolbox_label, R0.toolbox_label) && strcmp(R1.deval_mode, R0.deval_mode) && ...
    (~ref || (strcmp(R1.toolbox_label, 'reference') && strcmp(R1.deval_mode, 'replica'))) && ...
    isequal(bytes(fullfile(d, 'gen1', 'sa_bvp4c_traced.m')), bytes(fullfile(d, 'gen_real', 'sa_bvp4c_traced.m')));
rows = sa_check(rows, T, ['1: a byte-identical fake toolbox is accepted with this toolbox''s label and deval mode ' ...
    '(reference, replica on the reference toolbox); the copy equals the one from the real toolbox'], ok1, ...
    sprintf('"%s" %s %s', e1, field(R1, 'toolbox_label'), field(R1, 'deval_mode')));
% ---- 2 one changed comment byte; resume refused with the toolbox message ----------------------------
f2 = fake(tb, fullfile(d, 'fake2'), H); fakes{end+1} = f2;
L2 = lines;
c = find(startsWith(strtrim(L2), '%') & cellfun(@numel, L2) >= 3, 1);
ch = L2{c}(2);
if isletter(ch), L2{c}(2) = char(bitxor(double(ch), 32)); else, L2{c}(2) = 'x'; end
put(fullfile(f2, 'bvp4c.m'), strjoin(L2, newline));
[e2, R2] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen2'), f2));
addpath(E.stub);
fac = fullfile(d, 'factory'); mkdir(fac);
put(fullfile(fac, 'sa_model_stub_t102.m'), strjoin({
    'function M = sa_model_stub_t102()'
    '% T102 test fixture: the stub model carrying the solver component of a scratch traced copy'
    'M = sa_model_stub_two_pool();'
    'M.id = ''stub_t102'';'
    'M.version = [''stub_t102/1 (test fixture)'' sa_bvp4c_setup(''component'', getappdata(0, ''SA_T102_MODEL''), getappdata(0, ''SA_T102_SOLVER''), getappdata(0, ''SA_T102_TB''))];'
    'end'}, newline));
addpath(fac);
c0 = onCleanup(@rm_apps);
setappdata(0, 'SA_T102_MODEL', fullfile(E.pkg, 'models', 'original_polarity'));
st = struct('model', 'stub_t102', 'budget', struct('max_calls', 3));
K = fullfile(d, 'run1');
setappdata(0, 'SA_T102_SOLVER', fullfile(d, 'gen1')); setappdata(0, 'SA_T102_TB', f1);
c1 = sa_bvp4c_setup('component', fullfile(E.pkg, 'models', 'original_polarity'), fullfile(d, 'gen1'), f1);
[ek, ~] = err(@() sa_m3_run(K, [1 1], 'settings', st));
setappdata(0, 'SA_T102_SOLVER', fullfile(d, 'gen2')); setappdata(0, 'SA_T102_TB', f2);
c2 = sa_bvp4c_setup('component', fullfile(E.pkg, 'models', 'original_polarity'), fullfile(d, 'gen2'), f2);
st2 = st; st2.resume = true;
er = ''; mr = '';
try
    sa_m3_run(K, [1 1], 'settings', st2);
catch ME
    er = ME.identifier; mr = ME.message;
end
b1 = regexp(c1, 'toolbox_bvp4c=([0-9a-f]+)', 'tokens', 'once'); b2 = regexp(c2, 'toolbox_bvp4c=([0-9a-f]+)', 'tokens', 'once');
% the message the wrapper (run_seeded_atlas.m) prints: the resume status of the folder under 2. The model
% version is also a derived settings leaf, so the direct seeded_atlas_run resume stops first on E_SETTINGS
% (model_version); the status lists every refusal, the toolbox change included
U = load(fullfile(sa_layout('state', K), 'settings_user.mat')); U = U.settings;
U.run_dir = K; U.resume = true;
M2 = sa_model('stub_t102');
Z2 = sa_resume('status', K, sa_settings('resolve', U, M2), struct('model_version', M2.version));
want = sprintf('model version: the MATLAB toolbox changed since this run started (bvp4c %s -> %s)', b1{1}, b2{1});
ok2 = isempty(e2) && strcmp(R2.toolbox_label, 'self-checked') && isempty(ek) && ~strcmp(c1, c2) && ...
    ~strcmp(b1{1}, b2{1}) && startsWith(er, 'SA:RESUME:E_') && strcmp(Z2.state, 'refused') && contains(Z2.message, want);
rows = sa_check(rows, T, ['2: one changed comment byte in bvp4c.m is accepted as self-checked; the solver: component ' ...
    'differs, so a stub run made under 1 is refused on resume under 2, and the resume status (the wrapper''s ' ...
    'message) names the toolbox change'], ok2, sprintf('"%s" %s; run "%s"; resume %s: %s || status %s: %s', e2, ...
    field(R2, 'toolbox_label'), ek, er, mr, Z2.id, Z2.message));
% ---- 3 re-indented anchors: tolerant path ------------------------------------------------------------
f3 = fake(tb, fullfile(d, 'fake3'), H); fakes{end+1} = f3;
L3 = lines; L3{i2} = ['  ' L3{i2}]; L3{i3} = ['  ' L3{i3}];
put(fullfile(f3, 'bvp4c.m'), strjoin(L3, newline));
[e3, R3] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen3'), f3));
G3 = lines_of(fullfile(d, 'gen3', 'sa_bvp4c_traced.m'));
ok3 = isempty(e3) && strcmp(R3.anchor_initial, sprintf('tolerant line %d', i2)) && ...
    strcmp(R3.anchor_proposal, sprintf('tolerant line %d', i3)) && ...
    strcmp(G3{i2 + 1}, ['  ' regexp(lines{i2}, '^\s*', 'match', 'once') strtrim(P.initial_insert)]) && ...
    strcmp(G3{i3 + 2}, ['  ' regexp(lines{i3}, '^\s*', 'match', 'once') strtrim(P.note_insert)]);
rows = sa_check(rows, T, ['3: re-indented anchor lines are accepted through the tolerant path (recorded); the inserted ' ...
    'lines take the anchors'' indentation'], ok3, sprintf('"%s" %s / %s', e3, field(R3, 'anchor_initial'), ...
    field(R3, 'anchor_proposal')));
% ---- 4 duplicated proposal anchor ----------------------------------------------------------------------
f4 = fake(tb, fullfile(d, 'fake4'), H); fakes{end+1} = f4;
L4 = [lines(1:i3), lines(i3), lines(i3 + 1:end)];
put(fullfile(f4, 'bvp4c.m'), strjoin(L4, newline));
[e4, ~, m4] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen4'), f4));
ok4 = strcmp(e4, 'SA:SOLVER:E_ANCHOR') && ~isempty(regexp(m4, 'new_profile exact \d+ / tolerant 2', 'once')) && ...
    exist(fullfile(d, 'gen4'), 'dir') ~= 7;
rows = sa_check(rows, T, '4: a duplicated proposal anchor in the main body is E_ANCHOR with the counts; nothing written', ...
    ok4, sprintf('%s: %s', e4, m4));
% ---- 5 removed helper -------------------------------------------------------------------------------------
f5 = fake(tb, fullfile(d, 'fake5'), H); fakes{end+1} = f5;
delete(fullfile(f5, 'private', H{end}));
[e5, ~, m5] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen5'), f5));
rows = sa_check(rows, T, '5: a removed private helper is E_CLOSURE naming it', strcmp(e5, 'SA:SOLVER:E_CLOSURE') && ...
    contains(m5, strrep(H{end}, '.m', '')), sprintf('%s: %s', e5, m5));
% ---- 6 numerical change in residual --------------------------------------------------------------------------
f6 = fake(tb, fullfile(d, 'fake6'), H); fakes{end+1} = f6;
L6 = lines;
k6 = find(~cellfun(@isempty, regexp(L6, '^lobw3\s*=\s*32/45;', 'once')), 1);
if ~isempty(k6)
    L6{k6} = strrep(L6{k6}, '32/45', '32/44');          % a Lobatto weight of the residual
else                                                    % another release: scale the residual on return
    kr = find(~cellfun(@isempty, regexp(L6, '^function .*=\s*residual\s*\(', 'once')), 1);
    kn = kr + find(~cellfun(@isempty, regexp(L6(kr + 1:end), '^function\>', 'once')), 1);
    L6 = [L6(1:kn - 1), {'res = 1.5 * res;'}, L6(kn:end)];
    k6 = kn;
end
put(fullfile(f6, 'bvp4c.m'), strjoin(L6, newline));
[e6n, ~, m6n] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen6_nopath'), f6));
p0 = path;
addpath(f6, '-begin');
w6 = which('bvp4c');
[e6, R6] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen6'), f6));
path(p0);
ok6 = isscalar(k6) && strcmp(w6, fullfile(f6, 'bvp4c.m')) && isempty(e6) && strcmp(R6.toolbox_label, 'self-checked') && ...
    strcmp(e6n, 'SA:SOLVER:E_SELFCHECK') && contains(m6n, 'case (a)');
rows = sa_check(rows, T, ['6: a numerical change in the residual subfunction is accepted when the changed bvp4c is the ' ...
    'user''s (fake folder first on the path): the self-check certifies the patch, not the numerics; against another ' ...
    'bvp4c the copy is refused E_SELFCHECK (case a)'], ok6, sprintf('line %s; which %s; "%s" %s; without path %s: %s', ...
    mat2str(k6), w6, e6, field(R6, 'toolbox_label'), e6n, m6n));
% ---- 7 changed deval byte -----------------------------------------------------------------------------------
f7 = fake(tb, fullfile(d, 'fake7'), H); fakes{end+1} = f7;
fid = fopen(fullfile(f7, 'deval.m'), 'a'); fprintf(fid, '%% changed\n'); fclose(fid);
[e7, R7] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen7'), f7));
C0 = sa_bvp4c_setup('cache');
Ct = struct('recorder', 'x', 'deval_mode', field(R7, 'deval_mode'), 'label', 'x');
sa_bvp4c_setup('cache', Ct);
s = synth();
xs = [s.x, (s.x(1:end-1) + s.x(2:end)) / 2];
xs = xs([2:2:end, 1:2:end]);
[A1, A2] = deval(s, xs); [B1, B2] = sa_deval(s, xs);
[A3] = deval(s, xs, [2 1]); [B3] = sa_deval(s, xs, [2 1]);
if isempty(C0), sa_bvp4c_setup('reset_cache'); else, sa_bvp4c_setup('cache', C0); end
ok7 = isempty(e7) && strcmp(R7.deval_mode, 'toolbox') && strcmp(R7.toolbox_label, 'self-checked') && ...
    contains(R7.selfcheck_cases, 'e=not run') && isequal({A1, A2, A3}, {B1, B2, B3});
rows = sa_check(rows, T, ['7: a changed byte in deval.m is accepted with deval_mode toolbox (case e not run); in that ' ...
    'mode sa_deval equals the toolbox deval'], ok7, sprintf('"%s" %s %s: %s', e7, field(R7, 'toolbox_label'), ...
    field(R7, 'deval_mode'), field(R7, 'deval_reason')));
% ---- 8 the older new_profile signature ------------------------------------------------------------------------
f8 = fake(tb, fullfile(d, 'fake8'), H); fakes{end+1} = f8;
L8 = lines;
L8{i3} = regexprep(L8{i3}, ',\s*canRemovePoints\s*\)', ')');
kf = find(~cellfun(@isempty, regexp(L8, '^function .*=\s*new_profile\(', 'once')), 1);
L8 = [L8(1:kf), {'if nargin < 9, canRemovePoints = true; end'}, L8(kf + 1:end)];
put(fullfile(f8, 'bvp4c.m'), strjoin(L8, newline));
p0 = path;
addpath(f8, '-begin');
[e8, R8] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen8'), f8));
path(p0);
ok8 = ~contains(L8{i3}, 'canRemovePoints') && isscalar(kf) && isempty(e8) && ...
    strcmp(R8.anchor_proposal, sprintf('tolerant line %d', i3)) && strcmp(R8.toolbox_label, 'self-checked');
rows = sa_check(rows, T, ['8: a new_profile call without canRemovePoints (the older signature) is accepted through the ' ...
    'tolerant path'], ok8, sprintf('"%s" %s', e8, field(R8, 'anchor_proposal')));
% ---- 9 .p sibling --------------------------------------------------------------------------------------------
f9 = fake(tb, fullfile(d, 'fake9'), H); fakes{end+1} = f9;
put(fullfile(f9, 'bvp4c.p'), 'not a real p-file');
[e9, ~, m9] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen9'), f9));
rows = sa_check(rows, T, '9: a .p sibling of bvp4c.m is E_TOOLBOX (an unpatchable bvp4c); nothing written', ...
    strcmp(e9, 'SA:SOLVER:E_TOOLBOX') && exist(fullfile(d, 'gen9'), 'dir') ~= 7, sprintf('%s: %s', e9, m9));
% ---- 10 CRLF ---------------------------------------------------------------------------------------------------
f10 = fake(tb, fullfile(d, 'fake10'), H); fakes{end+1} = f10;
put(fullfile(f10, 'bvp4c.m'), strrep(LF, newline, [char(13) newline]));
[e10, R10] = err(@() sa_bvp4c_setup('generate', fullfile(d, 'gen10'), f10));
g10 = char(bytes_or_empty(fullfile(d, 'gen10', 'sa_bvp4c_traced.m')));
g1 = char(bytes(fullfile(d, 'gen1', 'sa_bvp4c_traced.m')));
ok10 = isempty(e10) && contains(R10.selfcheck_cases, 'g=pass') && strcmp(g10, strrep(g1, newline, [char(13) newline])) && ...
    startsWith(R10.line_endings, 'CR kept');
rows = sa_check(rows, T, ['10: a CRLF-converted bvp4c.m is accepted; the copy keeps CRLF on every line, inserted ones ' ...
    'included (it equals the LF copy with LF -> CRLF); self-check case g holds'], ok10, sprintf('"%s" %s; %s', e10, ...
    field(R10, 'line_endings'), field(R10, 'selfcheck_cases')));
% ---- the installed copy is untouched; fakes removed --------------------------------------------------------------
for k = 1:numel(fakes), rmdir(fakes{k}, 's'); end
rows = sa_check(rows, T, 'the installed traced copy is untouched; the fake toolboxes are deleted', ...
    isequal(bytes_or_empty(inst), inst0) && ~any(cellfun(@(f) exist(f, 'dir') == 7, fakes)), inst);
rmpath(fac); rmpath(E.stub);
end

% =========================================================================

function s = synth()
x = linspace(0, 2, 9);
s = struct('solver', 'bvp4c', 'x', x, 'y', [sin(x); cos(3 * x)], 'yp', [cos(x); -3 * sin(3 * x)]);
end

function rm_apps()
for k = {'SA_T102_MODEL', 'SA_T102_SOLVER', 'SA_T102_TB'}
    if isappdata(0, k{1}), rmappdata(0, k{1}); end
end
end

function v = field(R, f)
if isstruct(R) && isfield(R, f), v = R.(f); else, v = ''; end
end

function fk = fake(tb, fk, helpers)
mkdir(fullfile(fk, 'private'));
copyfile(fullfile(tb, 'bvp4c.m'), fk); copyfile(fullfile(tb, 'deval.m'), fk);
for k = 1:numel(helpers), copyfile(fullfile(tb, 'private', helpers{k}), fullfile(fk, 'private')); end
fileattrib(fk, '+w', 'u', 's');
end

function n = num(s)
n = str2double(regexp(s, '\d+$', 'match', 'once'));
end

function [e, out, msg] = err(fn)
e = ''; out = []; msg = '';
try
    if nargout > 1, out = fn(); else, fn(); end
catch ME
    e = ME.identifier; msg = ME.message;
end
end

function put(f, txt)
fid = fopen(f, 'w'); fwrite(fid, uint8(txt), 'uint8'); fclose(fid);
end

function b = bytes(f)
fid = fopen(f, 'r'); b = fread(fid, Inf, '*uint8').'; fclose(fid);
end

function b = bytes_or_empty(f)
if isfile(f), b = bytes(f); else, b = []; end
end

function L = lines_of(f)
L = strsplit(char(bytes(f)), newline, 'CollapseDelimiters', false);
end
