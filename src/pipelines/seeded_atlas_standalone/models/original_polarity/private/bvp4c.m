function sol = bvp4c(varargin)
%BVP4C  M6 private forwarder (packet v0.2 section 2.3).
%   Every bvp4c call of the functions in this model folder (the verbatim
%   corrector copy included, with no byte change) reaches sa_bvp4c_traced:
%   the generated, renamed copy of the pinned R2025b bvp4c.m whose only
%   edits record the initial mesh and every mesh proposal (sa_mesh_trace).
%   A reader of the corrector sees bvp4c; this forwarder runs. Pinned by T61
%   and recorded in the runtime evidence (seeded_atlas_run).
sol = sa_bvp4c_traced(varargin{:});
end
