function h = sa_sha256(file)
%SA_SHA256  Lower-case hex SHA-256 of a file's bytes (input file identity).
%   Streams the file through java.security.MessageDigest in 16 MB chunks.

fid = fopen(file, 'r');
if fid < 0
    error('SA:SHA256', 'Cannot open %s.', file);
end
cleaner = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('SHA-256');
while true
    b = fread(fid, 16 * 2^20, '*uint8');
    if isempty(b), break; end
    md.update(typecast(b, 'int8'));
end
d = typecast(md.digest(), 'uint8');
h = lower(reshape(dec2hex(d, 2).', 1, []));
end
