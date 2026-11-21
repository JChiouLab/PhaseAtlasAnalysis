function varargout = sa_progress(action, varargin)
%SA_PROGRESS  Progress lines and solver output on the screen (display only).
%
%   Nothing here enters a solve, a stored record, a key or the order of work:
%   it only prints. The state is per MATLAB process (persistent) and is cleared
%   when the outermost caller closes it.
%
%   SA_PROGRESS('begin', interval_s)
%       Open one level. The first open starts the elapsed-time clock and the
%       interval timer; every open sets the interval (seconds between status
%       lines; Inf = stage lines only).
%   SA_PROGRESS('end')
%       Close one level; closing the last level clears the state.
%   SA_PROGRESS('stage', name, info)
%       Print a stage line now: "[progress] stage: <name> | <info>".
%   tf = SA_PROGRESS('due')
%       True when a status line is due (the interval has passed since the
%       last progress line); false when no level is open.
%   SA_PROGRESS('status', info)
%       Print a status line for the current stage and restart the interval.
%   k = SA_PROGRESS('kind')
%       The current stage's kind (the code the caller passed in info.kind).
%   SA_PROGRESS('solver_text', txt, show_all)
%       Print the text one solve wrote to the command window (captured by
%       the search driver sa_search). show_all = true (report.solver_warnings):
%       everything, verbatim. Otherwise every bvp4c mesh-limit warning block
%       ("Unable to meet the tolerance without using more than N mesh
%       points ...") is left out and anything else is printed unchanged; the
%       first time a block is left out a one-line note says where the
%       warning is recorded. The warning itself is recorded per call by the
%       coupled bvp4c solve in every case (warning_id, warning_msg).
%   n = SA_PROGRESS('hidden')
%       Mesh-limit warning blocks left out since the outermost open.
%   txt = SA_PROGRESS('strip', txt)
%       The text without its mesh-limit warning blocks (no printing).
%
%   info (all fields optional): kind, within (text: progress inside the
%   stage), calls, max_calls, roots, replay (text). Fields left out keep
%   their last values (a plot or trace stage still shows the run's calls and
%   roots).

persistent P
if isempty(P), P = blank(); end
switch action
    case 'begin'
        if P.depth == 0
            P = blank();
            P.t0 = tic; P.last = tic;
        end
        P.depth = P.depth + 1;
        if ~isempty(varargin) && ~isempty(varargin{1}), P.interval = double(varargin{1}); end
    case 'end'
        P.depth = max(P.depth - 1, 0);
        if P.depth == 0, P = blank(); end
    case 'stage'
        P = merge(P, info_of(varargin, 2));
        P.stage = varargin{1};
        if isfield(info_of(varargin, 2), 'kind'), P.kind = varargin{2}.kind; else, P.kind = ''; end
        fprintf('%s\n', line_text(P, 'stage'));
        if P.depth > 0, P.last = tic; end
    case 'due'
        varargout{1} = P.depth > 0 && isfinite(P.interval) && toc(P.last) >= P.interval;
    case 'status'
        P = merge(P, info_of(varargin, 1));
        fprintf('%s\n', line_text(P, 'status'));
        if P.depth > 0, P.last = tic; end
    case 'kind'
        varargout{1} = P.kind;
    case 'solver_text'
        txt = varargin{1};
        show_all = numel(varargin) >= 2 && isequal(varargin{2}, true);
        if isempty(txt), return; end
        if show_all
            fprintf('%s', txt);
            return
        end
        [rest, n] = strip(txt);
        if n > 0
            if ~P.noted
                fprintf(['[progress] note: bvp4c''s mesh-limit warnings are not printed; each call records ' ...
                    'its warning (attempts.csv: the convergence_warning check). report.solver_warnings = ' ...
                    'true prints them.\n']);
                P.noted = true;
            end
            P.hidden = P.hidden + n;
        end
        if ~isempty(strtrim(rest)), fprintf('%s', rest); end
    case 'hidden'
        varargout{1} = P.hidden;
    case 'strip'
        [varargout{1}, varargout{2}] = strip(varargin{1});
    otherwise
        error('SA:PROGRESS', 'Unknown sa_progress action %s.', action);
end
end

% =========================================================================

function P = blank()
P = struct('depth', 0, 't0', [], 'last', [], 'interval', 300, 'stage', '', 'kind', '', 'within', '', ...
    'calls', NaN, 'max_calls', NaN, 'roots', NaN, 'replay', '', 'noted', false, 'hidden', 0);
end

function I = info_of(args, k)
I = struct();
if numel(args) >= k && isstruct(args{k}), I = args{k}; end
end

function P = merge(P, I)
for f = {'within', 'calls', 'max_calls', 'roots', 'replay'}
    if isfield(I, f{1}), P.(f{1}) = I.(f{1}); end
end
if ~isfield(I, 'within'), P.within = ''; end
if ~isfield(I, 'replay'), P.replay = ''; end
end

function t = line_text(P, what)
parts = {sprintf('[progress] %s: %s', what, P.stage)};
if ~isempty(P.within), parts{end+1} = P.within; end
if ~isempty(P.replay), parts{end+1} = P.replay; end
if isfinite(P.calls)
    if isfinite(P.max_calls)
        parts{end+1} = sprintf('calls %d of %d', P.calls, P.max_calls);
    else
        parts{end+1} = sprintf('calls %d', P.calls);
    end
end
if isfinite(P.roots), parts{end+1} = sprintf('roots %d', P.roots); end
if ~isempty(P.t0), parts{end+1} = ['elapsed ' hms(toc(P.t0))]; end
parts{end+1} = datestr(now, 'HH:MM:SS'); %#ok<DATST,TNOW1>
t = strjoin(parts, ' | ');
end

function s = hms(x)
x = floor(x);
s = sprintf('%d:%02d:%02d', floor(x / 3600), floor(mod(x, 3600) / 60), mod(x, 60));
end

function [rest, n] = strip(txt)
% One block: 'Warning: Unable to meet the tolerance without using more
% than <N> mesh points.' ... 'requested accuracy is <tol>.', with the
% command window's warning markers around it ('[' and ']', each may be
% followed by a backspace character), trailing blanks and the line end.
% Words may be wrapped onto new lines (MATLAB wraps long warnings), so word
% gaps match any white space.
bs = char(8);
pat = ['[\[' bs ']*Warning:\s+Unable\s+to\s+meet\s+the\s+tolerance\s+without\s+using\s+more\s+than\s+\d+\s+' ...
    'mesh\s+points\..*?requested\s+accuracy\s+is\s+\S+?\.[\]' bs ']*[ \t]*(\r?\n|$)'];
[s, e] = regexp(txt, pat, 'start', 'end');
n = numel(s);
if n == 0, rest = txt; return; end
keep = true(1, numel(txt));
for k = 1:n, keep(s(k):e(k)) = false; end
rest = txt(keep);
end
