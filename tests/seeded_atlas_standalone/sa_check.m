function rows = sa_check(rows, test, assertion, pass, detail)
%SA_CHECK  Append one assertion row (test, assertion, pass, detail).
if nargin < 5, detail = ''; end
rows(end+1) = struct('test', test, 'assertion', assertion, ...
    'pass', logical(pass), 'detail', detail);
end
