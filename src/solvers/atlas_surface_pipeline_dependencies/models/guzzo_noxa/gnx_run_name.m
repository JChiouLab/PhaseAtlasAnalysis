function [name, sha8, canon] = gnx_run_name(P)
%GNX_RUN_NAME  Run-folder name for an M. xanthus (Guzzo_noXA) parameter set.
%
%   [name, sha8, canon] = GNX_RUN_NAME(P)     P = gnx_model('parameters', ...)
%   name = guzzo_noxa_L<L>_s<s>_dba<f>_rR<rho_R>_rB<rho_B>_h<h>_res<res>_<sha8>,
%   numbers printed with %g, '.' -> 'p', '-' -> 'm'; sha8 = first 8 hex of the
%   SHA-256 of the canonical normalized struct (every field, in order,
%   'name=%.17g' joined by ';'), so an unlisted parameter change changes the
%   name. The fission yeast folder pattern starts with 'a' (atlas_run_name) and cannot
%   collide.

f = fieldnames(P);
parts = cell(1, numel(f));
for k = 1:numel(f)
    parts{k} = sprintf('%s=%.17g', f{k}, P.(f{k}));
end
canon = strjoin(parts, ';');
h = gnx_model('sha256_bytes', uint8(canon));
sha8 = h(1:8);
name = sprintf('guzzo_noxa_L%s_s%s_dba%s_rR%s_rB%s_h%s_res%d_%s', num(P.L), num(P.r_speed), ...
    num(P.dba_factor), num(P.rho_R), num(P.rho_B), num(P.dn_half_range), P.res, sha8);
end

function s = num(v)
s = strrep(strrep(sprintf('%g', v), '.', 'p'), '-', 'm');
end
