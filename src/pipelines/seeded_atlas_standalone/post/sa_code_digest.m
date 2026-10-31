function varargout = sa_code_digest(action, varargin)
%SA_CODE_DIGEST  Comment-free token digest of MATLAB source files.
%
%   T = SA_CODE_DIGEST('tokens', file)
%       The parse tree of the file (MATLAB's mtree) as node kinds and
%       strings, with the comment nodes (COMMENT, BLKCOM, CELLMARK) dropped.
%       T.kinds, T.strings (cell rows), T.n, T.has_err (true when the tree
%       holds an ERR node, i.e. the file does not parse), T.digest (SHA-256
%       hex of the stream). Two files with equal streams compile to the same
%       code: a comment edit leaves the stream unchanged, any code edit
%       (operator, literal, identifier, string, line structure inside
%       brackets) changes it.
%   T = SA_CODE_DIGEST('text', txt)        the same for source text
%   h = SA_CODE_DIGEST('files', files)
%       One digest over several files: SHA-256 of the lines
%       "<file name>:<token digest>" in the given order. A file that does not
%       parse is an error (SA:CODE_DIGEST:E_PARSE).
%   [same, why] = SA_CODE_DIGEST('equal', fileA, fileB)
%       Token-stream equality of two files (and neither holds an ERR node).

switch action
    case 'tokens', varargout{1} = tokens(mtree(varargin{1}, '-file'));
    case 'text',   varargout{1} = tokens(mtree(varargin{1}));
    case 'files',  varargout{1} = files_digest(varargin{1});
    case 'equal',  [varargout{1}, varargout{2}] = equal(varargin{:});
    otherwise, error('SA:CODE_DIGEST', 'Unknown sa_code_digest action %s.', action);
end
end

function T = tokens(tree)
k = reshape(tree.kinds, 1, []);
s = reshape(tree.strings, 1, []);
drop = ismember(k, {'COMMENT', 'BLKCOM', 'CELLMARK'});
k = k(~drop); s = s(~drop);
s(cellfun(@isempty, s)) = {''};
T = struct('kinds', {k}, 'strings', {s}, 'n', numel(k), 'has_err', any(strcmp(k, 'ERR')));
T.digest = sha_text(strjoin(strcat(k, char(31), s), char(30)));
end

function h = files_digest(files)
if ischar(files), files = {files}; end
lines = cell(1, numel(files));
for f = 1:numel(files)
    T = tokens(mtree(files{f}, '-file'));
    if T.has_err
        error('SA:CODE_DIGEST:E_PARSE', '%s does not parse (mtree ERR node).', files{f});
    end
    [~, nm, ex] = fileparts(files{f});
    lines{f} = sprintf('%s%s:%s', nm, ex, T.digest);
end
h = sha_text(strjoin(lines, newline));
end

function [same, why] = equal(a, b)
A = tokens(mtree(a, '-file'));
B = tokens(mtree(b, '-file'));
why = '';
if A.has_err || B.has_err
    same = false; why = 'parse error (ERR node)'; return
end
same = isequal(A.kinds, B.kinds) && isequal(A.strings, B.strings);
if ~same
    n = min(A.n, B.n);
    d = find(~strcmp(A.kinds(1:n), B.kinds(1:n)) | ~strcmp(A.strings(1:n), B.strings(1:n)), 1);
    if isempty(d), d = n + 1; end
    why = sprintf('token streams differ at token %d of %d / %d', d, A.n, B.n);
end
end

function h = sha_text(t)
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(unicode2native(t, 'UTF-8'), 'int8'));
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
end
