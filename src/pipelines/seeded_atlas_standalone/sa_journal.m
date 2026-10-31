function varargout = sa_journal(action, varargin)
%SA_JOURNAL  Append-only run journal.
%
%   SA_JOURNAL('start', run_dir, id, unit_key)   before a call
%   SA_JOURNAL('end',   run_dir, id, outcome)    after the attempt file is in place
%   SA_JOURNAL('note',  run_dir, text)           resume / override / stop lines
%   J = SA_JOURNAL('read', run_dir)
%       J.starts, J.ends: struct arrays (id, text) in file order; J.lines.
%
%   One line per event; each line is one open / write / close, so a kill
%   loses at most the line in flight (a partial last line is ignored by
%   'read' and reported in J.partial).

switch action
    case 'start', append(varargin{1}, sprintf('start %s %s', varargin{2}, varargin{3}));
    case 'end',   append(varargin{1}, sprintf('end %s %s', varargin{2}, varargin{3}));
    case 'note',  append(varargin{1}, sprintf('note %s', varargin{2}));
    case 'read',  varargout{1} = read_journal(varargin{1});
    otherwise, error('SA:JOURNAL', 'Unknown sa_journal action %s.', action);
end
end

function append(run_dir, line)
line = strrep(line, newline, ' ');
fid = fopen(fullfile(run_dir, 'journal.txt'), 'a');
if fid < 0, error('SA:JOURNAL', 'Cannot open the journal in %s.', run_dir); end
fprintf(fid, '%s\n', line);
fclose(fid);
end

function J = read_journal(run_dir)
f = fullfile(run_dir, 'journal.txt');
J = struct('starts', struct('id', {}, 'text', {}), 'ends', struct('id', {}, 'text', {}), ...
    'lines', {{}}, 'partial', '');
if exist(f, 'file') ~= 2, return; end
txt = fileread(f);
if ~isempty(txt) && txt(end) ~= newline
    k = find(txt == newline, 1, 'last');
    if isempty(k), J.partial = txt; txt = ''; else, J.partial = txt(k+1:end); txt = txt(1:k); end
end
L = strsplit(txt, newline);
L = L(~cellfun(@isempty, L));
J.lines = L;
for k = 1:numel(L)
    t = regexp(L{k}, '^(start|end) (A\d+) (.*)$', 'tokens', 'once');
    if isempty(t), continue; end
    if strcmp(t{1}, 'start')
        J.starts(end+1) = struct('id', t{2}, 'text', t{3});
    else
        J.ends(end+1) = struct('id', t{2}, 'text', t{3});
    end
end
end
