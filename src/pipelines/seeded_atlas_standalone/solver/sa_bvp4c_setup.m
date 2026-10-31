function varargout = sa_bvp4c_setup(action, varargin)
%SA_BVP4C_SETUP  Generate, self-check and verify the traced bvp4c copy
%   (M6 packet v0.2 section 2.3; milestone 9 packet v0.2 section 5, E1).
%
%   R = SA_BVP4C_SETUP('generate')
%   R = SA_BVP4C_SETUP('generate', out_dir, toolbox_dir)
%       Makes out_dir/sa_bvp4c_traced.m (default out_dir: this solver/
%       folder) from the user's own toolbox bvp4c.m (default: the bvp4c.m
%       under matlabroot that which('bvp4c','-all') lists), checks it, and
%       writes the generation record out_dir/sa_bvp4c_generated.txt. Steps:
%         1. locate the toolbox; a .p or MEX sibling of bvp4c.m (or no .m
%            under matlabroot) is SA:SOLVER:E_TOOLBOX (an unpatchable bvp4c);
%         2. read the bytes; split on LF, each line keeping its own ending
%            (a CR before the LF is kept, anchors match with or without it,
%            inserted lines take the anchor line's ending);
%         3. three edits at text anchors, each matched exactly first and by
%            a tolerant pattern second, unique in the main body (the lines
%            before the first subfunction), else SA:SOLVER:E_ANCHOR with the
%            counts:
%              - the function line (the first code line) is replaced by the
%                generated header, templated on the release, so that on
%                R2025b the bytes are those of milestone 6;
%              - sa_mesh_trace('initial', N); after N = length(x);
%              - sa_mesh_trace('note', N); after the [N,x,y,mbcidx] =
%                new_profile(...) line (after the first anchor);
%            the local solver_name = 'bvp4c' line is unchanged (recorded,
%            not required);
%         4. the private closure, computed: every identifier of the copy that
%            names a file of <toolbox>/private (or one of the five reference
%            helpers), transitively over the copied helpers, is copied byte
%            for byte into out_dir/private/; a named helper that cannot be
%            copied is SA:SOLVER:E_CLOSURE; on the reference toolbox the
%            closure must be exactly the five reference helpers;
%         5. write, rehash, self-check (below), then the record. A failed
%            self-check deletes the generated files and is
%            SA:SOLVER:E_SELFCHECK naming the case. A folder that cannot be
%            written is SA:SOLVER:E_WRITE.
%       The generated files are NOT committed (git-ignored): no MathWorks
%       source is redistributed; anyone regenerates them from their own
%       MATLAB.
%
%   The self-check (about 12 tiny bvp4c calls, single thread, warnings left
%   on): on fixed toy problems the traced copy and the user's own toolbox
%   bvp4c (resolved from this folder, which has no forwarder) must give
%   isequal outputs, with a same-process control (the toolbox bvp4c twice);
%   every case also meets the per-call trace consistency rule of
%   sa_record_call (one initial mesh; P = max of the proposals; the last
%   proposal equals stats.nmeshpoints when a solution returned; else
%   nmeshpoints = initial = numel(sol.x); the returned mesh <= max(initial,
%   P)). Cases: (a) a nonlinear problem refining from 5 points (>= 1
%   proposal, stats.maxres present); (b) the same with NMax too small (both
%   warn MATLAB:bvp4c:RelTolNotMet; the last proposal > NMax, over_nmax);
%   (c) an exactly singular problem (a repeated boundary row; both throw
%   MATLAB:bvp4c:SingJac); (d) one unknown parameter (sol.parameters too);
%   (e) deval, only on the reference deval.m and private/ntrp3h.m: the
%   replica sa_deval against the toolbox deval -> deval_mode 'replica',
%   else 'toolbox' (never a refusal); (f) information: sa_validity's
%   bvp4c_residual over stats.maxres; (g) the text diff: the copy equals the
%   toolbox file except the function line and the two inserted lines; (h)
%   the identifier catalog: RelTolNotMet, IllCondJac and SingJac resolve
%   (getString(message(id, ...)), zero calls). The self-check certifies the
%   patch, not the numerics of a release.
%
%   R = SA_BVP4C_SETUP('verify')
%   R = SA_BVP4C_SETUP('verify', solver_dir, toolbox_dir)
%       The run-start check (zero calls): the record exists and names this
%       MATLAB (release, version, architecture) and this toolbox folder, the
%       toolbox files have the recorded SHA-256 (else E_TOOLBOX: the MATLAB
%       toolbox changed since setup, run setup_atlas), and every generated
%       file exists with its recorded SHA-256 (else E_TRACED).
%   st = SA_BVP4C_SETUP('status')   verify without throwing: .ok .id .message .record
%   [txt, A] = SA_BVP4C_SETUP('edit_text', txt [, release])   the three edits on a text
%   t = SA_BVP4C_SETUP('component', model_dir [, solver_dir, toolbox_dir])
%       The model-version 'solver:' component: verify, then SHA-256 prefixes
%       (12 hex) of the traced copy, its helpers, sa_deval.m, the toolbox
%       bvp4c.m and deval.m (measured, from the record), the forwarders, and
%       the tags deval=<mode> and tb=<label>. On the installed solver folder
%       it also loads the per-process cache (the recorder string of
%       sa_mesh_trace and the deval mode of sa_deval), so the per-call paths
%       never read a file. E_MODEL is never waived.
%   s = SA_BVP4C_SETUP('recorder')     'sa_bvp4c_traced/1 (bvp4c <release> <sha12>)'
%   m = SA_BVP4C_SETUP('deval_mode')   'replica' or 'toolbox'
%       Both from the cache; when component() has not run in this process
%       (stub runs, tests) the cache is filled once from the verified
%       installed record, or 'unverified' / 'toolbox' without one.
%   C = SA_BVP4C_SETUP('cache'); SA_BVP4C_SETUP('cache', C);
%   SA_BVP4C_SETUP('reset_cache')      read, set (tests) or empty the cache
%   R = SA_BVP4C_SETUP('record', dir)  the parsed record of a folder
%   tb = SA_BVP4C_SETUP('locate')      the toolbox folder of this MATLAB
%   P = SA_BVP4C_SETUP('pins')         the R2025b (Update 4) reference values:
%       a toolbox matching them is labelled 'reference' (the commissioned
%       scope), any other 'self-checked'. They are a record, not a gate.

switch action
    case 'generate',    [varargout{1:nargout}] = generate(varargin{:});
    case 'verify',      [varargout{1:nargout}] = verify(varargin{:});
    case 'status',      varargout{1} = status(varargin{:});
    case 'edit_text',   [varargout{1:max(1, nargout)}] = edit_text(varargin{:});
    case 'component',   varargout{1} = component(varargin{:});
    case 'recorder',    varargout{1} = cached('recorder');
    case 'deval_mode',  varargout{1} = cached('deval_mode');
    case 'reset_cache', cache_store('clear');
    case 'cache'
        if isempty(varargin), varargout{1} = cache_store('get'); else, cache_store('set', varargin{1}); end
    case 'record',      varargout{1} = read_record(varargin{:});
    case 'locate',      varargout{1} = locate_toolbox();
    case 'pins',        varargout{1} = pins();
    otherwise, error('SA:SOLVER', 'Unknown sa_bvp4c_setup action %s.', action);
end
end

% =========================================================================

function P = pins()
% The reference record: R2025b (Update 4) toolbox files and the copy made from them.
P = struct();
P.toolbox_dir = fullfile(matlabroot, 'toolbox', 'matlab', 'funfun');
P.bvp4c = '597839cef613775894e9f869fd8de7f95ef711a9ae0dee8f9dea855497d60d09';
P.deval = 'c8920dd29039a0b78b66460edaa13c9291d26e2db15f7128ef3e449e389d2319';
P.helpers = struct('name', {'bvparguments.m', 'bvpfunctions.m', 'bvpsingular.m', 'ntrp3h.m', 'odenumjac.m'}, ...
    'sha256', {'fe5fbbefa1387f0d099ec3c012279ba3b2c6506a5e0fcfb5953b8d5855027343', ...
               'd8666d0c9fef066ce2b14a475c31b60728250ccb28ed066929a0dd931bfa630a', ...
               '3fbe1e17e32bcb2b16485c82b97f6c4e04445ad9b3bfc502323ef22d56f4307f', ...
               '5c44de3580a2d15f53ef7551b541b17cb086212a07820c783f31350c3732af31', ...
               'e57bed8ff1a4697c87bc327fd8c12c35dbbfdead9c74bb82ca87afb633273118'});
P.traced = '56e8727d1b1e117dac59f5509521dada12e1d352526cfa7ec9587414aa2e9889';
P.release = 'R2025b';
P.function_line = 'function sol = bvp4c(ode, bc, solinit, options, varargin)';
P.function_new = header('R2025b');
P.initial_anchor = 'N = length(x);      % number of mesh points';
P.initial_insert = 'sa_mesh_trace(''initial'', N);   % sa_bvp4c_setup: inserted (M6 mesh trace)';
P.note_anchor = '        [N,x,y,mbcidx] = new_profile(n,x,y,yp,res,mbcidx,rtol,Nmax,canRemovePoints);';
P.note_insert = '        sa_mesh_trace(''note'', N);   % sa_bvp4c_setup: inserted (M6 mesh trace)';
P.solver_name_line = 'solver_name = ''bvp4c'';';
P.function_pattern = ['^function\s+sol\s*=\s*bvp4c\s*\(ode\s*,\s*bc\s*,\s*solinit\s*,\s*options\s*,' ...
    '\s*varargin\s*\)'];
P.initial_pattern = '^\s*N\s*=\s*length\(x\);';
P.note_pattern = '^\s*\[N,x,y,mbcidx\]\s*=\s*new_profile\(';
P.record_file = 'sa_bvp4c_generated.txt';
end

function h = header(release)
h = sprintf(['function sol = sa_bvp4c_traced(ode, bc, solinit, options, varargin) ' ...
    '%% GENERATED by sa_bvp4c_setup from %s bvp4c.m (M6): do not edit'], release);
end

function r = current_release()
r = ['R' version('-release')];
end

function d = here()
d = fileparts(mfilename('fullpath'));
end

% ---- the per-process cache -------------------------------------------------------------------

function out = cache_store(op, val)
persistent C
switch op
    case 'get',   out = C;
    case 'set',   C = val; out = C;
    case 'clear', C = []; out = [];
end
end

function v = cached(field)
C = cache_store('get');
if isempty(C)
    C = cache_store('set', fill_cache());
end
v = C.(field);
end

function C = fill_cache()
% Only when component() has not run in this process: the installed record, verified.
try
    C = cache_of(verify());
catch
    C = struct('recorder', 'sa_bvp4c_traced/1 (bvp4c unverified)', 'deval_mode', 'toolbox', ...
        'label', 'unverified');
end
end

function C = cache_of(R)
C = struct('recorder', sprintf('sa_bvp4c_traced/1 (bvp4c %s %s)', R.release, R.toolbox_bvp4c(1:12)), ...
    'deval_mode', R.deval_mode, 'label', R.toolbox_label);
end

% ---- locate, verify, status -------------------------------------------------------------------

function tb = locate_toolbox()
w = which('bvp4c', '-all');
if ischar(w), w = {w}; end
root = [matlabroot filesep];
hit = w(startsWith(w, root));
if isempty(hit)
    error('SA:SOLVER:E_TOOLBOX', 'No toolbox bvp4c under matlabroot (%s); which lists: %s.', matlabroot, ...
        strjoin(w, ', '));
end
[tb, ~, ext] = fileparts(hit{1});
if ~strcmpi(ext, '.m')
    error('SA:SOLVER:E_TOOLBOX', ['The toolbox bvp4c is %s, not a text .m file: it cannot be patched with the ' ...
        'mesh trace, so the pipeline cannot run on this MATLAB.'], hit{1});
end
end

function check_shadow(tb)
s = dir(fullfile(tb, 'bvp4c.*'));
s = s(~[s.isdir]);
if ~(isscalar(s) && strcmp(s(1).name, 'bvp4c.m'))
    error('SA:SOLVER:E_TOOLBOX', ['The toolbox folder %s has no bvp4c.m or a .p/MEX file beside it (%s): ' ...
        'the bvp4c MATLAB runs cannot be patched with the mesh trace.'], tb, strjoin({s.name}, ', '));
end
if ~isfile(fullfile(tb, 'deval.m'))
    error('SA:SOLVER:E_TOOLBOX', 'The toolbox folder %s has no deval.m.', tb);
end
end

function R = verify(solver_dir, tb)
if nargin < 1 || isempty(solver_dir), solver_dir = here(); end
if nargin < 2 || isempty(tb), tb = locate_toolbox(); end
P = pins();
f = fullfile(solver_dir, 'sa_bvp4c_traced.m');
if ~isfile(f)
    error('SA:SOLVER:E_TRACED', ['The generated solver copy %s is missing: run setup_atlas (or ' ...
        'sa_bvp4c_setup(''generate'')) (USER_GUIDE section 2).'], f);
end
if ~isfile(fullfile(solver_dir, P.record_file))
    error('SA:SOLVER:E_TRACED', ['The generation record %s is missing: regenerate with setup_atlas (or ' ...
        'sa_bvp4c_setup(''generate'')).'], fullfile(solver_dir, P.record_file));
end
R = read_record(solver_dir);
now_id = sprintf('%s %s %s %s', current_release(), version, computer('arch'), tb);
was_id = sprintf('%s %s %s %s', R.release, R.version, R.arch, R.toolbox_dir);
if ~strcmp(now_id, was_id)
    error('SA:SOLVER:E_TOOLBOX', ['The MATLAB toolbox changed since setup (recorded: %s; now: %s): run ' ...
        'setup_atlas (USER_GUIDE section 2).'], was_id, now_id);
end
check_shadow(tb);
tf = [{fullfile(tb, 'bvp4c.m'), fullfile(tb, 'deval.m')}, ...
    cellfun(@(h) fullfile(tb, 'private', h), R.helpers, 'UniformOutput', false)];
th = [{R.toolbox_bvp4c, R.toolbox_deval}, cellfun(@(h) R.(key('helper', h)), R.helpers, 'UniformOutput', false)];
for k = 1:numel(tf)
    if ~isfile(tf{k}) || ~strcmp(sha(tf{k}), th{k})
        error('SA:SOLVER:E_TOOLBOX', ['The MATLAB toolbox changed since setup: %s is missing or differs from ' ...
            'its recorded SHA-256: run setup_atlas (USER_GUIDE section 2).'], tf{k});
    end
end
if ~strcmp(sha(f), R.generated_traced)
    error('SA:SOLVER:E_TRACED', ['The generated solver copy %s differs from its recorded SHA-256 (edited?): ' ...
        'regenerate with setup_atlas (or sa_bvp4c_setup(''generate'')).'], f);
end
for k = 1:numel(R.helpers)
    g = fullfile(solver_dir, 'private', R.helpers{k});
    if ~isfile(g) || ~strcmp(sha(g), R.(key('generated', R.helpers{k})))
        error('SA:SOLVER:E_TRACED', ['The copied solver helper %s is missing or differs from its recorded ' ...
            'SHA-256: regenerate with setup_atlas (or sa_bvp4c_setup(''generate'')).'], g);
    end
end
end

function st = status(varargin)
st = struct('ok', false, 'id', '', 'message', '', 'record', []);
try
    st.record = verify(varargin{:});
    st.ok = true;
catch ME
    st.id = ME.identifier; st.message = ME.message;
end
end

% ---- generate ---------------------------------------------------------------------------------

function R = generate(out_dir, tb)
P = pins();
if nargin < 1 || isempty(out_dir), out_dir = here(); end
if nargin < 2 || isempty(tb), tb = locate_toolbox(); end
installed = strcmp(out_dir, here());
check_shadow(tb);
src = read_text(fullfile(tb, 'bvp4c.m'));
[new, A] = edit_text(src, current_release());
[names, missing] = closure(tb, new);
if ~isempty(missing)
    error('SA:SOLVER:E_CLOSURE', ['The traced copy names private helper(s) %s that %s/private does not hold ' ...
        'as a copyable .m file.'], strjoin(missing, ', '), tb);
end
H = struct('bvp4c', sha(fullfile(tb, 'bvp4c.m')), 'deval', sha(fullfile(tb, 'deval.m')));
if strcmp(H.bvp4c, P.bvp4c) && ~isequal(sort(names), sort(strrep({P.helpers.name}, '.m', '')))
    error('SA:SOLVER:E_CLOSURE', ['On the reference bvp4c.m the private closure must be the five reference ' ...
        'helpers; it is: %s.'], strjoin(names, ', '));
end
% ---- write
try
    if exist(fullfile(out_dir, 'private'), 'dir') ~= 7, mkdir(fullfile(out_dir, 'private')); end
    old = dir(fullfile(out_dir, 'private', '*.m'));
    for k = 1:numel(old)
        if ~any(strcmp(old(k).name, strcat(names, '.m'))), delete(fullfile(old(k).folder, old(k).name)); end
    end
    delete_if(fullfile(out_dir, P.record_file));
    write_text(fullfile(out_dir, 'sa_bvp4c_traced.m'), new);
    for k = 1:numel(names)
        [ok, msg] = copyfile(fullfile(tb, 'private', [names{k} '.m']), ...
            fullfile(out_dir, 'private', [names{k} '.m']), 'f');
        if ~ok, error('SA:SOLVER:E_COPY', 'Cannot copy %s: %s', names{k}, msg); end
    end
catch ME
    error('SA:SOLVER:E_WRITE', ['Cannot write the traced bvp4c copy into %s (%s). The package folder must be ' ...
        'writable: copy the package to a writable folder, or generate elsewhere with ' ...
        'sa_bvp4c_setup(''generate'', out_dir, toolbox_dir) to inspect the copy.'], out_dir, ME.message);
end
rehash;
% ---- hashes of what was written
hp = cellfun(@(n) sha(fullfile(tb, 'private', [n '.m'])), names, 'UniformOutput', false);
gp = cellfun(@(n) sha(fullfile(out_dir, 'private', [n '.m'])), names, 'UniformOutput', false);
gt = sha(fullfile(out_dir, 'sa_bvp4c_traced.m'));
nt = fullfile(tb, 'private', 'ntrp3h.m');
if isfile(nt), H.ntrp3h = sha(nt); else, H.ntrp3h = 'missing'; end
ref = strcmp(H.bvp4c, P.bvp4c) && strcmp(H.deval, P.deval) && numel(names) == numel(P.helpers);
for k = 1:numel(P.helpers)
    j = strcmp(names, strrep(P.helpers(k).name, '.m', ''));
    ref = ref && any(j) && strcmp(hp{find(j, 1)}, P.helpers(k).sha256);
end
ntref = P.helpers(strcmp({P.helpers.name}, 'ntrp3h.m')).sha256;
% ---- self-check
try
    S = selfcheck(out_dir, tb, src, new, A, strcmp(H.deval, P.deval) && strcmp(H.ntrp3h, ntref));
catch ME
    S = struct('pass', false, 'failed', 'run', 'reason', sprintf('%s: %s', ME.identifier, ME.message));
end
if ~S.pass
    delete_if(fullfile(out_dir, 'sa_bvp4c_traced.m'));
    for k = 1:numel(names), delete_if(fullfile(out_dir, 'private', [names{k} '.m'])); end
    rehash;
    if installed, cache_store('clear'); end
    error('SA:SOLVER:E_SELFCHECK', 'The traced bvp4c copy failed the self-check, case (%s): %s', S.failed, S.reason);
end
% ---- the record
labels = {'self-checked', 'reference'};
L = {
    'record', 'sa_bvp4c_generated/1'
    'release', current_release()
    'version', version
    'arch', computer('arch')
    'toolbox_dir', tb
    'toolbox_bvp4c', H.bvp4c
    'toolbox_deval', H.deval
    'toolbox_ntrp3h', H.ntrp3h
    'helpers', strjoin(strcat(names, '.m'), ',')};
for k = 1:numel(names), L(end+1, :) = {key('helper', [names{k} '.m']), hp{k}}; end %#ok<AGROW>
L(end+1, :) = {'generated_traced', gt};
for k = 1:numel(names), L(end+1, :) = {key('generated', [names{k} '.m']), gp{k}}; end %#ok<AGROW>
L = [L; {
    'anchor_function', sprintf('%s line %d', A.function_path, A.function_line)
    'anchor_initial', sprintf('%s line %d', A.initial_path, A.initial_line)
    'anchor_proposal', sprintf('%s line %d', A.note_path, A.note_line)
    'solver_name_lines', sprintf('%d', A.solver_name_count)
    'main_body_end', sprintf('%d', A.body_end)
    'line_endings', A.endings
    'toolbox_label', labels{1 + ref}
    'deval_mode', S.deval_mode
    'deval_rule', 'reference-hash+case-e'
    'deval_reason', S.deval_reason
    'selfcheck', sprintf('pass(%d/%d)', S.n_pass, S.n_cases)
    'selfcheck_cases', S.summary
    'selfcheck_control', S.control
    'selfcheck_threads', sprintf('%d', S.threads)
    'selfcheck_calls', sprintf('%d', S.calls)
    'time', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'))}];
txt = strjoin(cellfun(@(a, b) sprintf('%s: %s', a, b), L(:, 1), L(:, 2), 'UniformOutput', false), newline);
write_text(fullfile(out_dir, P.record_file), [txt newline]);
R = read_record(out_dir);
if installed, cache_store('clear'); end        % the next use reads the new record
end

function k = key(prefix, file)
k = [prefix '_' regexprep(strrep(file, '.m', ''), '\W', '_')];
end

function [new, A] = edit_text(txt, release)
P = pins();
if nargin < 2 || isempty(release), release = current_release(); end
L = strsplit(txt, newline, 'CollapseDelimiters', false);
B = regexprep(L, '\r$', '');
A = struct();
isc = ~cellfun(@(s) isempty(strtrim(s)) || startsWith(strtrim(s), '%'), B);
i1 = find(isc, 1);
if isempty(i1), i1 = 1; end
ex1 = strcmp(B{i1}, P.function_line);
if ~(ex1 || ~isempty(regexp(B{i1}, P.function_pattern, 'once')))
    error('SA:SOLVER:E_ANCHOR', ['bvp4c.m anchors: the first code line (line %d) is not the bvp4c function ' ...
        'line (exact 0, tolerant 0).'], i1);
end
sub = find(~cellfun(@isempty, regexp(B, '^function\>', 'once')));
sub = sub(sub > i1);
if isempty(sub), be = numel(B); else, be = sub(1) - 1; end
body = (i1 + 1):be;
e2 = body(strcmp(B(body), P.initial_anchor));
t2 = body(~cellfun(@isempty, regexp(B(body), P.initial_pattern, 'once')));
e3 = body(strcmp(B(body), P.note_anchor));
t3 = body(~cellfun(@isempty, regexp(B(body), P.note_pattern, 'once')));
if isscalar(e2), i2 = e2; p2 = 'exact'; elseif isscalar(t2), i2 = t2; p2 = 'tolerant'; else, i2 = []; p2 = ''; end
i3 = []; p3 = '';
if ~isempty(i2)
    if isscalar(e3) && e3 > i2, i3 = e3; p3 = 'exact'; elseif isscalar(t3) && t3 > i2, i3 = t3; p3 = 'tolerant'; end
end
if isempty(i2) || isempty(i3)
    error('SA:SOLVER:E_ANCHOR', ['bvp4c.m anchors not found exactly once in the main body (lines %d-%d): ' ...
        'initial mesh exact %d / tolerant %d, new_profile exact %d / tolerant %d (after the initial-mesh ' ...
        'anchor: %d).'], i1 + 1, be, numel(e2), numel(t2), numel(e3), numel(t3), ...
        sum(t3 > max([i2, 0])) * ~isempty(i2));
end
cr = @(k) L{k}(numel(B{k}) + 1:end);
ind = @(k) regexp(B{k}, '^\s*', 'match', 'once');
ins2 = [ind(i2) strtrim(P.initial_insert) cr(i2)];
ins3 = [ind(i3) strtrim(P.note_insert) cr(i3)];
L{i1} = [header(release) cr(i1)];
L = [L(1:i2), {ins2}, L(i2+1:i3), {ins3}, L(i3+1:end)];
new = strjoin(L, newline);
ncr = sum(~cellfun(@isempty, regexp(L, '\r$', 'once')));
if ncr == 0, A.endings = 'LF'; else, A.endings = sprintf('CR kept on %d lines', ncr); end
A.function_line = i1; A.function_path = tern(ex1, 'exact', 'tolerant');
A.initial_line = i2; A.initial_path = p2;
A.note_line = i3; A.note_path = p3;
A.solver_name_count = sum(strcmp(strtrim(B), P.solver_name_line));
A.body_end = be;
end

function [names, missing] = closure(tb, txt)
P = pins();
pv = dir(fullfile(tb, 'private', '*'));
pv = pv(~[pv.isdir]);
base = regexprep({pv.name}, '\.[^.]*$', '');
cand = unique([base, strrep({P.helpers.name}, '.m', '')]);
names = {}; missing = {};
queue = {txt};
while ~isempty(queue)
    ids = unique(regexp(queue{1}, '[A-Za-z]\w*', 'match'));
    queue(1) = [];
    hs = setdiff(intersect(ids, cand), [names, missing]);
    for j = 1:numel(hs)
        f = fullfile(tb, 'private', [hs{j} '.m']);
        if isfile(f)
            names{end+1} = hs{j}; %#ok<AGROW>
            queue{end+1} = read_text(f); %#ok<AGROW>
        else
            missing{end+1} = hs{j}; %#ok<AGROW>
        end
    end
end
names = sort(names);
end

% ---- the self-check ---------------------------------------------------------------------------

function S = selfcheck(out_dir, tb, src, new, A, deval_ref)
S = struct('pass', false, 'failed', '', 'reason', '', 'deval_mode', 'toolbox', 'deval_reason', '', ...
    'n_pass', 0, 'n_cases', 0, 'summary', '', 'control', 'deterministic', 'threads', NaN, 'calls', 0);
n0 = maxNumCompThreads;
c1 = onCleanup(@() maxNumCompThreads(n0));
maxNumCompThreads(1);
S.threads = maxNumCompThreads;
p0 = path;
c2 = onCleanup(@() path(p0));
addpath(out_dir, '-begin');
rehash;
want = fullfile(out_dir, 'sa_bvp4c_traced.m');
if ~strcmp(which('sa_bvp4c_traced'), want)
    S.failed = 'resolution'; S.reason = sprintf('sa_bvp4c_traced resolves to %s, not %s', which('sa_bvp4c_traced'), want);
    return
end
traced = @sa_bvp4c_traced;
plain = @bvp4c;
% the toy problems (neutral names and constants)
q = 60;
f_a = @(x, y) [y(2, :); q * (y(1, :) + y(1, :).^3 / 3)];
g_a = @(ya, yb) [ya(1) - 1; yb(1) - 0.3];
m_a = bvpinit(linspace(0, 1, 5), [1; 0]);
g_c = @(ya, yb) [ya(1) - 1; ya(1) - 1];
f_d = @(x, y, p) [y(2, :); -p * y(1, :); y(1, :).^2];
g_d = @(ya, yb, p) [ya(1); ya(2) - 1; yb(1); ya(3)];
m_d = bvpinit(linspace(0, 1, 7), @(x) [sin(3 * x); 3 * cos(3 * x); 0 * x], 9);
nm_b = 12;
spec = {
    'a', f_a, g_a, m_a, [], floor(10000 / 2)
    'b', f_a, g_a, m_a, bvpset('NMax', nm_b), nm_b
    'c', f_a, g_c, m_a, [], floor(10000 / 2)
    'd', f_d, g_d, m_d, [], floor(10000 / 3)};
out = struct();
summ = {};
for k = 1:size(spec, 1)
    c = spec{k, 1};
    U1 = run_one(plain, spec(k, 2:5), false, spec{k, 6});
    U2 = run_one(plain, spec(k, 2:5), false, spec{k, 6});
    T = run_one(traced, spec(k, 2:5), true, spec{k, 6});
    S.calls = S.calls + 3;
    out.(c) = U1;
    if ~same(U1, U2)
        S.control = sprintf('non-deterministic in this session (case %s)', c);
        S.failed = 'control'; S.reason = sprintf(['the toolbox bvp4c gave two different results on the same ' ...
            'inputs in this session (case %s): not a defect of the patch, but nothing can be certified'], c);
        return
    end
    why = '';
    if ~same(U1, T), why = 'traced and toolbox outputs differ'; end
    if isempty(why), why = trace_rule(T); end
    switch c
        case 'a'
            if isempty(why) && ~isempty(T.err), why = ['threw ' T.err]; end
            if isempty(why) && T.trace.initial ~= 5, why = sprintf('initial mesh %g, not 5', T.trace.initial); end
            if isempty(why) && T.trace.n_proposals < 1, why = 'no mesh proposal'; end
            if isempty(why) && ~isfield(T.sol.stats, 'maxres'), why = 'sol.stats.maxres missing'; end
            tag = sprintf('a=%s(%d proposals)', tern(isempty(why), 'pass', 'fail'), T.trace.n_proposals);
        case 'b'
            if isempty(why) && ~(strcmp(U1.wid, 'MATLAB:bvp4c:RelTolNotMet') && strcmp(T.wid, U1.wid))
                why = sprintf('warnings %s / %s, not MATLAB:bvp4c:RelTolNotMet', U1.wid, T.wid);
            end
            if isempty(why) && isempty(T.sol), why = 'no solution returned'; end
            if isempty(why) && ~(~isempty(T.trace.proposals) && T.trace.proposals(end) > nm_b && ...
                    T.trace.over_nmax && T.trace.proposals(end) == T.sol.stats.nmeshpoints)
                why = sprintf('last proposal %s, over_nmax %d, nmeshpoints %g', mat2str(T.trace.proposals), ...
                    T.trace.over_nmax, T.sol.stats.nmeshpoints);
            end
            tag = sprintf('b=%s', tern(isempty(why), 'pass', 'fail'));
        case 'c'
            if isempty(why) && ~(strcmp(U1.err, 'MATLAB:bvp4c:SingJac') && strcmp(T.err, U1.err))
                why = sprintf('errors "%s" / "%s", not MATLAB:bvp4c:SingJac', U1.err, T.err);
            end
            tag = sprintf('c=%s', tern(isempty(why), 'pass', 'fail'));
        case 'd'
            if isempty(why) && (~isempty(T.err) || ~isequal(T.sol.parameters, U1.sol.parameters))
                why = 'parameters differ or the call threw';
            end
            tag = sprintf('d=%s', tern(isempty(why), 'pass', 'fail'));
    end
    summ{end+1} = tag; %#ok<AGROW>
    S.n_cases = S.n_cases + 1;
    if ~isempty(why)
        S.failed = c; S.reason = why; S.summary = strjoin(summ, ' ');
        return
    end
    S.n_pass = S.n_pass + 1;
end
% (e) deval
[S.deval_mode, S.deval_reason, ran, ne] = deval_case(out.a.sol, out.d.sol, deval_ref);
if ran
    S.n_cases = S.n_cases + 1;
    if strcmp(S.deval_mode, 'replica'), S.n_pass = S.n_pass + 1; end
    summ{end+1} = sprintf('e=%s(%d evaluations)', tern(strcmp(S.deval_mode, 'replica'), 'pass', ...
        'differs:deval=toolbox'), ne);
else
    summ{end+1} = 'e=not run'; %#ok<AGROW>
end
% (f) information only
r = NaN;
try
    if exist('sa_validity', 'file') == 2
        M = struct('ode', @(x, Y, sol, theta) f_a(x, Y));
        r = sa_validity('bvp4c_residual', M, out.a.sol, struct(), struct('abstol', 1e-6, 'reltol', 1e-3)) / ...
            out.a.sol.stats.maxres;
    end
catch
end
summ{end+1} = sprintf('f=info(residual ratio %.6g)', r);
% (g) the text diff
why = text_diff(src, new, A, fullfile(out_dir, 'sa_bvp4c_traced.m'));
S.n_cases = S.n_cases + 1;
summ{end+1} = sprintf('g=%s', tern(isempty(why), 'pass', 'fail'));
if ~isempty(why)
    S.failed = 'g'; S.reason = why; S.summary = strjoin(summ, ' ');
    return
end
S.n_pass = S.n_pass + 1;
% (h) the identifier catalog
bad = {};
for id = {'MATLAB:bvp4c:RelTolNotMet', 'MATLAB:bvp4c:IllCondJac', 'MATLAB:bvp4c:SingJac'}
    if ~resolves(id{1}), bad{end+1} = id{1}; end %#ok<AGROW>
end
S.n_cases = S.n_cases + 1;
summ{end+1} = sprintf('h=%s', tern(isempty(bad), 'pass', 'fail'));
S.summary = strjoin(summ, ' ');
if ~isempty(bad)
    S.failed = 'h'; S.reason = sprintf('message identifier(s) %s do not resolve in this MATLAB', strjoin(bad, ', '));
    return
end
S.n_pass = S.n_pass + 1;
S.pass = true;
end

function U = run_one(fn, a, traced, nmax)
U = struct('sol', [], 'err', '', 'wid', '', 'trace', []);
lastwarn('', '');
if traced, tok = sa_mesh_trace('reset'); end
try
    if isempty(a{4}), U.sol = fn(a{1}, a{2}, a{3}); else, U.sol = fn(a{1}, a{2}, a{3}, a{4}); end
catch ME
    U.err = ME.identifier;
end
[~, U.wid] = lastwarn;
if traced, U.trace = sa_mesh_trace('take', tok, nmax); end
end

function tf = same(U, V)
tf = strcmp(U.err, V.err) && strcmp(U.wid, V.wid) && isempty(U.sol) == isempty(V.sol);
if tf && ~isempty(U.sol)
    f = {'x', 'y', 'yp', 'stats', 'solver'};
    if isfield(U.sol, 'parameters') || isfield(V.sol, 'parameters'), f{end+1} = 'parameters'; end
    for k = 1:numel(f)
        tf = tf && isfield(U.sol, f{k}) && isfield(V.sol, f{k}) && isequal(U.sol.(f{k}), V.sol.(f{k}));
    end
end
end

function why = trace_rule(T)
% the production per-call consistency rule (sa_record_call trace_problem)
why = '';
tr = T.trace;
if tr.n_initial ~= 1, why = sprintf('%d initial meshes recorded', tr.n_initial); return; end
P = tr.max_proposed;
if ~isequaln(P, sa_mesh_trace('max_of', tr)), why = 'P differs from its proposals'; return; end
if ~isempty(T.sol)
    m = numel(T.sol.x);
    ns = NaN;
    if isfield(T.sol, 'stats') && isfield(T.sol.stats, 'nmeshpoints'), ns = double(T.sol.stats.nmeshpoints); end
    if ~isempty(tr.proposals)
        if ~(P >= ns && tr.proposals(end) == ns)
            why = sprintf('P %g, last proposal %g, stats.nmeshpoints %g', P, tr.proposals(end), ns); return
        end
    elseif ~(ns == tr.initial && tr.initial == m)
        why = sprintf('no proposal but nmeshpoints %g, initial %g, mesh %d', ns, tr.initial, m); return
    end
    if ~(m <= max(tr.initial, P))
        why = sprintf('returned mesh %d above max(initial %g, P %g)', m, tr.initial, P); return
    end
end
end

function [mode, reason, ran, ne] = deval_case(sa, sd, deval_ref)
mode = 'toolbox'; ran = false; ne = 0;
if ~deval_ref
    reason = 'the toolbox deval.m or private/ntrp3h.m is not the reference: the toolbox deval is used';
    return
end
try
    matlab.internal.math.discretize([0.25 0.75], [0 0.5 1], false);
catch
    reason = 'matlab.internal.math.discretize is not callable: the toolbox deval is used';
    return
end
ran = true;
C0 = cache_store('get');
base = C0;
if isempty(base), base = struct('recorder', '', 'deval_mode', '', 'label', ''); end
base.deval_mode = 'replica';
cache_store('set', base);
c = onCleanup(@() cache_store('set', C0));
bad = {};
sols = {sa, sd};
for s = 1:2
    sol = sols{s};
    x = sol.x;
    n = size(sol.y, 1);
    mid = (x(1:end-1) + x(2:end)) / 2;
    perm = [mid(2:2:end), fliplr(mid(1:2:end))];
    pts = {x, x(end), mid, perm, sort([x mid]), mid(:), [x(3) x(3) mid(2) mid(2) x(end)]};
    idxs = {[], 1, n, [n 1], [1 1 n], (1:n).'};
    for p = 1:numel(pts)
        for i = 1:numel(idxs)
            for nout = 1:2
                if isempty(idxs{i}), args = {sol, pts{p}}; else, args = {sol, pts{p}, idxs{i}}; end
                A = cell(1, nout); B = cell(1, nout);
                [A{:}] = deval(args{:});
                [B{:}] = sa_deval(args{:});
                ne = ne + 1;
                if ~isequal(A, B) || ~isequal(cellfun(@class, A, 'UniformOutput', false), ...
                        cellfun(@class, B, 'UniformOutput', false))
                    bad{end+1} = sprintf('solution %d points %d idx %d nout %d', s, p, i, nout); %#ok<AGROW>
                end
            end
        end
    end
end
if isempty(bad)
    mode = 'replica';
    reason = sprintf('reference deval.m and ntrp3h.m; case e equal on %d evaluations', ne);
else
    reason = sprintf('case e: the replica differs (%s): the toolbox deval is used', strjoin(bad(1:min(3, end)), '; '));
end
end

function why = text_diff(src, new, A, gen_file)
why = '';
G = strsplit(read_text(gen_file), newline, 'CollapseDelimiters', false);
if ~isequal(strjoin(G, newline), new), why = 'the written copy differs from the edited text'; return; end
L = strsplit(src, newline, 'CollapseDelimiters', false);
k2 = A.initial_line + 1; k3 = A.note_line + 2;
if numel(G) ~= numel(L) + 2, why = sprintf('%d lines, toolbox %d', numel(G), numel(L)); return; end
rest = G; rest([k2 k3]) = [];
d = find(~cellfun(@isequal, rest, L));
if ~isequal(d, A.function_line)
    why = sprintf('the copy differs from the toolbox file at lines %s, not only at the function line', mat2str(d));
    return
end
P = pins();
ok = startsWith(strtrim(G{k2}), strtrim(P.initial_insert)) && startsWith(strtrim(G{k3}), strtrim(P.note_insert)) && ...
    strcmp(G{k2 - 1}, L{A.initial_line}) && strcmp(G{k3 - 1}, L{A.note_line});
if ~ok, why = 'an inserted line is not at its anchor'; end
end

function tf = resolves(id)
tf = false;
for n = 0:6
    a = num2cell(ones(1, n));
    try
        s = getString(message(id, a{:}));
        tf = ischar(s) && ~isempty(s);
        return
    catch ME
        if strcmp(ME.identifier, 'MATLAB:builtins:MessageNotFound'), return; end
    end
end
end

% ---- component ----------------------------------------------------------------------------------

function t = component(model_dir, solver_dir, tb)
installed = nargin < 2 || isempty(solver_dir);
if installed, solver_dir = here(); end
if nargin < 3, tb = []; end
R = verify(solver_dir, tb);
pkg = fileparts(here());
parts = {};
add = @(name, f) sprintf('%s=%s', name, prefix(f));
parts{end+1} = add('sa_bvp4c_traced', fullfile(solver_dir, 'sa_bvp4c_traced.m'));
for k = 1:numel(R.helpers)
    parts{end+1} = add(strrep(R.helpers{k}, '.m', ''), fullfile(solver_dir, 'private', R.helpers{k})); %#ok<AGROW>
end
parts{end+1} = add('sa_deval', fullfile(here(), 'sa_deval.m'));
parts{end+1} = sprintf('toolbox_bvp4c=%s', R.toolbox_bvp4c(1:12));
parts{end+1} = sprintf('toolbox_deval=%s', R.toolbox_deval(1:12));
parts{end+1} = add('fwd_pkg_deval', fullfile(pkg, 'private', 'deval.m'));
parts{end+1} = add('fwd_bvp4c', fullfile(model_dir, 'private', 'bvp4c.m'));
parts{end+1} = add('fwd_deval', fullfile(model_dir, 'private', 'deval.m'));
parts{end+1} = sprintf('deval=%s', R.deval_mode);
parts{end+1} = sprintf('tb=%s', R.toolbox_label);
t = sprintf(' solver:%s', strjoin(parts, ','));
if installed, cache_store('set', cache_of(R)); end
end

function p = prefix(f)
if ~isfile(f)
    error('SA:SOLVER:E_TRACED', 'solver file %s is missing.', f);
end
h = sha(f);
p = h(1:12);
end

% ---- files ----------------------------------------------------------------------------------------

function R = read_record(d)
f = fullfile(d, pins().record_file);
L = strsplit(read_text(f), newline);
R = struct();
for k = 1:numel(L)
    m = regexp(L{k}, '^([A-Za-z]\w*): (.*)$', 'tokens', 'once');
    if ~isempty(m), R.(m{1}) = m{2}; end
end
need = {'release', 'version', 'arch', 'toolbox_dir', 'toolbox_bvp4c', 'toolbox_deval', 'helpers', ...
    'generated_traced', 'toolbox_label', 'deval_mode'};
if ~all(isfield(R, need))
    error('SA:SOLVER:E_TRACED', 'The generation record %s is incomplete: regenerate with setup_atlas.', f);
end
R.helpers = strsplit(R.helpers, ',');
if isequal(R.helpers, {''}), R.helpers = {}; end
end

function txt = read_text(f)
% bytes 1:1 as character codes 0-255 (ASCII, UTF-8 or any other bytes kept exactly)
fid = fopen(f, 'r');
if fid < 0, error('SA:SOLVER', 'Cannot open %s.', f); end
c = onCleanup(@() fclose(fid));
b = fread(fid, Inf, '*uint8');
txt = char(b(:).');
end

function write_text(f, txt)
tmp = [f '.tmp_gen'];
fid = fopen(tmp, 'w');
if fid < 0, error('SA:SOLVER', 'Cannot write %s.', tmp); end
fwrite(fid, uint8(txt), 'uint8');
fclose(fid);
[ok, msg] = movefile(tmp, f, 'f');
if ~ok, error('SA:SOLVER', 'Cannot move %s into place: %s', tmp, msg); end
end

function delete_if(f)
if isfile(f), delete(f); end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function h = sha(file)
fid = fopen(file, 'r');
if fid < 0, error('SA:SOLVER', 'Cannot open %s.', file); end
c = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('SHA-256');
b = fread(fid, Inf, '*uint8');
if ~isempty(b), md.update(typecast(b, 'int8')); end
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
end
