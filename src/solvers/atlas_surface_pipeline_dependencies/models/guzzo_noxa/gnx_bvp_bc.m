function res = gnx_bvp_bc(ya, yb, nR, nB)
%GNX_BVP_BC  12 boundary conditions of the Guzzo_noXA steady BVP (packet §4a).
%
%   res = GNX_BVP_BC(ya, yb, nR, nB)
%   No flux of the four fields at both ends (eight conditions; two are
%   redundant by conservation and absorbed by the shifts sigmaR, sigmaB),
%   accumulators mR(0) = mB(0) = 0, mR(ell) = nR, mB(ell) = nB.
%   Same conditions as bvp10_noxa.m:28-33 (bc12) with the accumulators in
%   rows 11-12.
res = [ya(2); ya(4); ya(6); ya(8);
       yb(2); yb(4); yb(6); yb(8);
       ya(11); ya(12);
       yb(11) - nR; yb(12) - nB];
end
