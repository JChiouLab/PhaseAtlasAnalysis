function out = setup_atlas(mode)
%SETUP_ATLAS  One command for a fresh copy of the package (milestone 9, E7).
%
%   SETUP_ATLAS            puts the seeded-atlas package folder
%                          (src/pipelines/seeded_atlas_standalone) on the
%                          path; if the traced bvp4c copy is missing, or its
%                          generation record does not match this MATLAB and
%                          its toolbox, makes it from your own MATLAB's
%                          bvp4c.m with the self-check
%                          (sa_bvp4c_setup('generate')); prints the one-line
%                          verdict.
%   SETUP_ATLAS('quiet')   the same, silent unless the copy is (re)made;
%                          run_seeded_atlas.m calls this before every run.
%   SETUP_ATLAS('check')   verifies only (no generation, no solve) and prints
%                          the one-line verdict, e.g.
%       setup_atlas: OK  R2025b (25.2.0.3150157) maca64  bvp4c 597839cef613  toolbox=reference  deval=replica  selfcheck=pass(7/7)
%                          or "setup_atlas: FAIL <check>: <reason>". This is
%                          the line to report when testing another MATLAB.
%   out = SETUP_ATLAS(...) also returns .ok, .line and .record.
%
%   toolbox=reference: the R2025b Update 4 files the package was commissioned
%   on; toolbox=self-checked: another bvp4c.m, accepted because the traced
%   copy passed the self-check on it (tested only on R2025b Update 4, macOS
%   arm64). deval=replica: the package's fast deval replica (only on the
%   reference deval.m); deval=toolbox: your own deval. A setup or generation
%   failure names the failed check (E_TOOLBOX, E_ANCHOR, E_CLOSURE,
%   E_SELFCHECK, E_WRITE) and is raised as an error, except in 'check'.
%   Only the package folder stays on the path: building a library directly
%   (atlas_build_libraries) needs the uncoupled pipeline folders too (see the
%   README).

if nargin < 1 || isempty(mode), mode = 'setup'; end
if ~any(strcmp(mode, {'setup', 'quiet', 'check'}))
    error('SA:SETUP', 'setup_atlas mode must be ''setup'', ''quiet'' or ''check'', not %s.', mode);
end
here = fileparts(mfilename('fullpath'));
pkg = fullfile(here, 'src', 'pipelines', 'seeded_atlas_standalone');
solver = fullfile(pkg, 'solver');
if ~contains([pathsep path pathsep], [pathsep pkg pathsep]), addpath(pkg); end
p0 = path;
c = onCleanup(@() path(p0));                % solver/ only while setting up (sa_model adds it per run)
if ~contains([pathsep path pathsep], [pathsep solver pathsep]), addpath(solver); end
st = sa_bvp4c_setup('status');
made = false;
if ~st.ok && ~strcmp(mode, 'check')
    fprintf('setup_atlas: making the traced bvp4c copy from %s (%s)\n', toolbox_text(), st.message);
    try
        sa_bvp4c_setup('generate');
        made = true;
    catch ME
        fprintf('%s\n', fail_line(ME.identifier, ME.message));
        rethrow(ME);
    end
    st = sa_bvp4c_setup('status');
end
if st.ok
    R = st.record;
    line = sprintf(['setup_atlas: OK  %s (%s) %s  bvp4c %s  toolbox=%s  deval=%s  selfcheck=%s'], R.release, ...
        regexp(R.version, '^\S+', 'match', 'once'), R.arch, R.toolbox_bvp4c(1:12), R.toolbox_label, ...
        R.deval_mode, field_or(R, 'selfcheck', '?'));
else
    line = fail_line(st.id, st.message);
end
if ~strcmp(mode, 'quiet') || made || ~st.ok
    fprintf('%s\n', line);
end
if ~st.ok && ~strcmp(mode, 'check')
    error(st.id, '%s', st.message);
end
if nargout > 0
    out = struct('ok', st.ok, 'line', line, 'record', st.record);
end
end

function s = fail_line(id, msg)
k = regexp(id, '[^:]+$', 'match', 'once');
if isempty(k), k = 'ERROR'; end
s = sprintf('setup_atlas: FAIL %s: %s', k, msg);
end

function s = toolbox_text()
try
    s = sa_bvp4c_setup('locate');
catch
    s = 'the MATLAB toolbox';
end
end

function v = field_or(R, f, d)
if isfield(R, f), v = R.(f); else, v = d; end
end
