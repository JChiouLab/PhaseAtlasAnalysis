function C = sa_seedability(M, lib, theta, L)
%SA_SEEDABILITY  Library seedability census on the target lattice.
%
%   C = SA_SEEDABILITY(M, lib, theta, L) for the lattice L of sa_lattice.
%   Rules (as the campaign census sca_pr0b_l14_mask):
%     * a node maps to the library's nearest native cell (module 'cell'
%       action); a node outside the native swept range is not seedable and
%       asserts nothing about existence;
%     * the nearest cell must lie within half a native spacing on every axis,
%       otherwise the census fails closed;
%     * family A_B is seedable at a node iff half type A exists at the node's
%       cell and half type B exists at the exact discrete mirror cell
%       (res+1-row, res+1-col): no interpolation;
%     * families are every ordered pair of the module's library types, and the
%       seedable ones are DERIVED from the masks, never typed.
%
%   Returned fields
%     names            1 x nf family names 'A_B' (library type order, A major)
%     counts           1 x nf seedable node counts
%     seedable         n_nodes x nf logical
%     in_range         n_nodes x 1 logical
%     cell, mirror     n_nodes x 2 native cells (NaN outside the range)
%     offset           n_nodes x 2 |node - cell| in native spacings per axis
%                      (NaN outside the range; used by targets.policy per_seed)
%     outside_reason   n_nodes x 1 cellstr
%     worst_offset     largest |node - cell| in native spacings, per axis

types = M.library_types;
n = L.n_nodes;
C = struct();
C.in_range = false(n, 1);
C.cell = nan(n, 2);
C.mirror = nan(n, 2);
C.offset = nan(n, 2);
C.outside_reason = repmat({''}, n, 1);
worst = [0 0];
for k = 1:n
    c = M.library('cell', lib, theta, L.coordinate(k, :));
    if c.in_range
        C.in_range(k) = true;
        C.cell(k, :) = c.cell;
        C.mirror(k, :) = lib.res + 1 - c.cell;
        C.offset(k, :) = c.offset;
        worst = max(worst, c.offset);
    else
        C.outside_reason{k} = c.reason;
    end
end
if any(worst > 0.5 * (1 + 1e-9))
    error('SA:SEEDABILITY:E_SPACING', ['A lattice node''s nearest native cell is more ' ...
        'than half a native spacing away (worst %g, %g spacings).'], worst(1), worst(2));
end
C.worst_offset = worst;

names = {};
for a = 1:numel(types)
    for b = 1:numel(types)
        names{end+1} = sprintf('%s_%s', types{a}, types{b}); %#ok<AGROW>
    end
end
nf = numel(names);
C.names = names;
C.seedable = false(n, nf);
f = 0;
for a = 1:numel(types)
    for b = 1:numel(types)
        f = f + 1;
        mA = lib.mask.(types{a});
        mB = lib.mask.(types{b});
        for k = find(C.in_range).'
            C.seedable(k, f) = mA(C.cell(k, 1), C.cell(k, 2)) && ...
                               mB(C.mirror(k, 1), C.mirror(k, 2));
        end
    end
end
C.counts = sum(C.seedable, 1);
C.family_types = cell(nf, 2);
f = 0;
for a = 1:numel(types)
    for b = 1:numel(types)
        f = f + 1;
        C.family_types(f, :) = {types{a}, types{b}};
    end
end
C.rule = ['nearest native cell (within half a spacing); right half at the exact ' ...
    'discrete mirror (res+1-row, res+1-col); outside the native range = not seedable'];
end
