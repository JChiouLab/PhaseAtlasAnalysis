function sa_write(kind, file, varargin)
%SA_WRITE  Atomic output writes: write a temporary sibling, then rename.
%
%   SA_WRITE('mat',  file, S)      save the fields of struct S (-v7.3 if large)
%   SA_WRITE('text', file, text)   write char text
%
%   A reader sees either the previous complete file or the new complete file,
%   never a partial one; an interruption loses at most the write in flight.

tmp = [file '.tmp_' sprintf('%d', feature('getpid'))];
switch kind
    case 'mat'
        S = varargin{1};
        w = whos('S');
        if w.bytes > 1.5e9
            save(tmp, '-struct', 'S', '-v7.3');
        else
            save(tmp, '-struct', 'S', '-v7');
        end
        % save appends .mat when the name lacks it
        if exist(tmp, 'file') ~= 2 && exist([tmp '.mat'], 'file') == 2
            tmp = [tmp '.mat'];
        end
    case 'text'
        fid = fopen(tmp, 'w');
        if fid < 0, error('SA:WRITE', 'Cannot open %s for writing.', tmp); end
        fwrite(fid, varargin{1}, 'char');
        fclose(fid);
    otherwise
        error('SA:WRITE', 'Unknown write kind %s.', kind);
end
[ok, msg] = movefile(tmp, file, 'f');
if ~ok
    error('SA:WRITE', 'Atomic rename to %s failed: %s', file, msg);
end
end
