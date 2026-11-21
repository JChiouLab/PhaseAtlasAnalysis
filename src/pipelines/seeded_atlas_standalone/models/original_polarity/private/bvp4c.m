function sol = bvp4c(varargin)
%BVP4C  Private forwarder: bvp4c calls of this model folder go to the traced copy.
%   Every bvp4c call of the functions in this model folder (including the
%   verbatim copy of the coupled-BVP solve) reaches sa_bvp4c_traced: a
%   generated, renamed copy of your MATLAB's bvp4c.m whose only edits record
%   the initial mesh and every mesh proposal (see sa_mesh_trace). The calling
%   code reads as plain bvp4c; this forwarder is what runs.
sol = sa_bvp4c_traced(varargin{:});
end
