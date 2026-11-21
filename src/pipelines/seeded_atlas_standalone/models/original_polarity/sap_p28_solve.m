function [sol, info] = sap_p28_solve(th, guess, tgt, opts)
%SAP_P28_SOLVE  ONE call of the coupled S. pombe BVP (28 states, both halves).
%
%   Exactly one bvp4c call per invocation; the caller counts calls against its
%   budget. Applies the checks recorded in the info.flag_G* fields: G1 the
%   solver returned; G2 positivity; G3a source antisymmetry; G4 concentration
%   continuity at the middle; G5 flux match at the middle; G6a accumulator mass;
%   G6b dense mass; G9 velocity consistency.
%
%   G3b (all sources zero), G7 (profile distance) and G8 (full-domain steady
%   residual) are ROOT-ONLY and are applied by the caller: opts.zero_forcing_root
%   = true adds G3b here; G7/G8 stay with the caller because they need a target
%   profile. A forced off-root solution is a constrained quasi-steady state, not
%   an equilibrium, so G8 must never be applied to it.

if nargin < 4 || isempty(opts), opts = struct(); end
opts = defaults(opts);
info = init_info(opts, tgt);

t0 = tic; sol = [];
% Capture, never suppress, the bvp4c convergence warnings.
prevw = warning;
warning('on','MATLAB:bvp4c:RelTolNotMet'); warning('on','MATLAB:bvp4c:IllCondJac');
warning('off','backtrace'); lastwarn('');
cl = onCleanup(@() warning(prevw)); %#ok<NASGU>
try
    s = bvp4c(@(x,y) sap_forced28('ode', x, y, th), ...
              @(ya,yb) sap_forced28('bc', ya, yb, th, tgt), guess, ...
              bvpset('RelTol', opts.reltol, 'AbsTol', opts.abstol, 'NMax', opts.nmax));
    if numel(s.x) >= opts.mesh_reject
        info.reason = sprintf('mesh blow-up (%d points >= %d)', numel(s.x), opts.mesh_reject);
        try
            info.rejected_diag = compact(sap_forced28('diagnostics', s, th, tgt));
            info.rejected_n_mesh = numel(s.x);
        catch
        end
        % Recording only. Placed OUTSIDE the try above so the statistics land even
        % when the diagnostics call throws.
        info = ca0r3_record(info, s, opts);
    elseif any(~isfinite(s.y(:)))
        info.reason = 'non-finite solution';
        info = ca0r3_record(info, s, opts);      % recording only
    else
        sol = s;
    end
catch ME
    info.reason = ME.message;
end
info.wall_s = toc(t0);
if isempty(sol), info.class = 'FAILED'; return; end

info.converged = true;
D = sap_forced28('diagnostics', sol, th, tgt);
info.diag = D;
info.n_mesh = D.n_mesh;

% ---- validity checks ----------------------------------------------------
info.flag_G1  = true;
info.flag_G2  = D.all_positive;
info.antisym_gate = max(opts.antisym_abs, opts.antisym_rel*D.s_max);
info.flag_G3a = D.antisym_max < info.antisym_gate;
info.flag_G3b = D.s_max < opts.s_tol;                       % evaluated always, gating only at roots
info.flag_G4  = max(abs(D.seam_conc_mismatch)) < opts.seam_tol;
info.flag_G5  = D.max_seam_flux < opts.seam_tol;
info.flag_G6a = D.G6a_residual < opts.g6a_tol;
info.flag_G6b = (D.G6b_residual < opts.g6b_tol) && (D.G6b_refine < opts.g6b_refine_tol);
info.vel_gate = max(opts.vel_abs, opts.vel_rel*D.vel_max);
info.flag_G9  = D.vel_discrepancy_max < info.vel_gate;

% ---- scalar summary -----------------------------------------------------
info.minimum = D.minimum;
info.s_GL = D.s_GL; info.s_PL = D.s_PL; info.s_GR = D.s_GR; info.s_PR = D.s_PR;
info.s_max = D.s_max;
info.antisym_max = D.antisym_max; info.antisym_rel = D.antisym_rel;
info.seam_conc_max = max(abs(D.seam_conc_mismatch));
info.seam_flux_max = D.max_seam_flux;
info.G6a_residual = D.G6a_residual;
info.G6b_residual = D.G6b_residual; info.G6b_refine = D.G6b_refine;
info.G6c_residual = D.G6c_residual;
info.identity_max = D.identity_max;
info.vel_flux = D.vel_flux; info.vel_shift = D.vel_shift;
info.vel_discrepancy_max = D.vel_discrepancy_max;
info.vel_discrepancy_rel = D.vel_discrepancy_rel;
info.J_G = D.J_G; info.J_P = D.J_P;
info.d_G = D.d_G; info.d_P = D.d_P;        % authoritative: accumulator endpoints
info.M_G = D.M_G; info.M_P = D.M_P;
info.bc_max_abs = D.bc_max_abs;

core = info.flag_G1 && info.flag_G2 && info.flag_G3a && info.flag_G4 && ...
       info.flag_G5 && info.flag_G6a && info.flag_G6b && info.flag_G9;
if opts.zero_forcing_root, core = core && info.flag_G3b; end
% MATLAB's `lastwarn` returns [MESSAGE, IDENTIFIER].
[wmsg_, wid_] = lastwarn;
info.convergence = sap_bvp_convergence(sol, opts, wid_, wmsg_);
info.bvp_converged = info.convergence.converged;
info.returned_n_mesh = info.convergence.returned_n_mesh;
info.stats_nmeshpoints = info.convergence.stats_nmeshpoints;
info.stats_maxres = info.convergence.stats_maxres;
core = core && info.bvp_converged;   % convergence is required
info.accepted_core = core;
info.zero_forcing_root = opts.zero_forcing_root;

if core
    info.class = 'ACCEPTED_CORE'; info.reason = '';
else
    r = {};
    if ~info.flag_G2,  r{end+1} = sprintf('G2 min %.3g', D.minimum); end
    if ~info.flag_G3a, r{end+1} = sprintf('G3a antisym %.3g (gate %.3g)', D.antisym_max, info.antisym_gate); end
    if opts.zero_forcing_root && ~info.flag_G3b, r{end+1} = sprintf('G3b |s| %.3g', D.s_max); end
    if ~info.flag_G4,  r{end+1} = sprintf('G4 seam conc %.3g', info.seam_conc_max); end
    if ~info.flag_G5,  r{end+1} = sprintf('G5 seam flux %.3g', info.seam_flux_max); end
    if ~info.flag_G6a, r{end+1} = sprintf('G6a accumulator %.3g', D.G6a_residual); end
    if ~info.flag_G6b, r{end+1} = sprintf('G6b dense %.3g / refine %.3g', D.G6b_residual, D.G6b_refine); end
    if ~info.flag_G9,  r{end+1} = sprintf('G9 velocity %.3g (gate %.3g)', D.vel_discrepancy_max, info.vel_gate); end
    info.reason = strjoin(r, '; ');
    info.class = 'REJECTED';
    % Recording only. The returned-but-rejected path also ships its statistics
    % and its iterate, so that every refusal, not only the mesh_reject one, can
    % be inspected afterwards.
    info = ca0r3_record(info, sol, opts);
end
end

% -------------------------------------------------------------------------
% Recording-only helper: it stores the solver statistics (stats.maxres and the
% proposed mesh size) and the rejected iterate of a refused call, so that a
% refusal can be inspected afterwards.
%
% It must be a part of this function, not a wrapper: on the mesh_reject path
% this function sets sol = [] and the rejected iterate `s` is discarded inside
% the function, so it is unreachable from any caller.
%
% It changes no equation, boundary condition, check, tolerance, NMax,
% mesh_reject, acceptance expression or `defaults` value, and nothing it
% writes is ever read by an acceptance decision.
% -------------------------------------------------------------------------
function info = ca0r3_record(info, s, opts)
D = struct('has_stats',false,'stats_maxres',NaN,'stats_nmeshpoints',NaN, ...
    'returned_n_mesh',NaN,'reltol',NaN,'nmax',NaN,'mesh_reject',NaN);
try
    D.returned_n_mesh = numel(s.x);
    if isfield(s,'stats') && isstruct(s.stats)
        if isfield(s.stats,'maxres'),      D.stats_maxres      = s.stats.maxres;      end
        if isfield(s.stats,'nmeshpoints'), D.stats_nmeshpoints = s.stats.nmeshpoints; end
        D.has_stats = isfield(s.stats,'maxres') && isfield(s.stats,'nmeshpoints');
    end
    D.reltol = opts.reltol; D.nmax = opts.nmax; D.mesh_reject = opts.mesh_reject;
catch
end
if isfield(info,'rejected_diag') && isstruct(info.rejected_diag)
    fn = fieldnames(D);
    for k = 1:numel(fn), info.rejected_diag.(fn{k}) = D.(fn{k}); end
else
    info.rejected_diag = D;
end
info.rejected_sol = s;                 % the iterate is retained
info.rejected_retained = true;
info.ca0r3_recording = ['Recording only: stats.maxres and the proposed mesh ' ...
    'persisted in info.rejected_diag, and the rejected iterate retained in info.rejected_sol. ' ...
    'No solver policy, tolerance, gate or acceptance is affected.'];
end

% =========================================================================

function o = defaults(o)
d = struct('reltol',1e-6,'abstol',1e-8,'nmax',20000,'mesh_reject',3000, ...
    's_tol',1e-7,'seam_tol',1e-6, ...
    'g6a_tol',1e-8,'g6b_tol',1e-7,'g6b_refine_tol',1e-9, ...
    'antisym_abs',1e-9,'antisym_rel',1e-6, ...
    'vel_abs',1e-9,'vel_rel',1e-6, ...
    'zero_forcing_root',false);
fn = fieldnames(d);
for k = 1:numel(fn)
    if ~isfield(o,fn{k}) || isempty(o.(fn{k})), o.(fn{k}) = d.(fn{k}); end
end
end

function c = compact(D)
c = struct('n_mesh',D.n_mesh,'s_GL',D.s_GL,'s_PL',D.s_PL,'s_GR',D.s_GR,'s_PR',D.s_PR, ...
    's_max',D.s_max,'antisym_max',D.antisym_max, ...
    'max_seam_conc',max(abs(D.seam_conc_mismatch)),'max_seam_flux',D.max_seam_flux, ...
    'G6a_residual',D.G6a_residual,'G6b_residual',D.G6b_residual, ...
    'identity_max',D.identity_max,'vel_flux',D.vel_flux,'vel_shift',D.vel_shift, ...
    'minimum',D.minimum,'d_G',D.d_G,'d_P',D.d_P);
end

function i = init_info(o, tgt)
i = struct('converged',false,'accepted_core',false,'class','FAILED', ...
    'reason','not solved','wall_s',NaN,'n_mesh',NaN,'minimum',NaN, ...
    's_GL',NaN,'s_PL',NaN,'s_GR',NaN,'s_PR',NaN,'s_max',NaN, ...
    'antisym_max',NaN,'antisym_rel',NaN,'antisym_gate',NaN, ...
    'seam_conc_max',NaN,'seam_flux_max',NaN, ...
    'G6a_residual',NaN,'G6b_residual',NaN,'G6b_refine',NaN,'G6c_residual',NaN, ...
    'identity_max',NaN,'vel_flux',[NaN NaN],'vel_shift',[NaN NaN], ...
    'vel_discrepancy_max',NaN,'vel_discrepancy_rel',NaN,'vel_gate',NaN, ...
    'J_G',NaN,'J_P',NaN,'d_G',NaN,'d_P',NaN,'M_G',NaN,'M_P',NaN,'bc_max_abs',NaN, ...
    'flag_G1',false,'flag_G2',false,'flag_G3a',false,'flag_G3b',false, ...
    'flag_G4',false,'flag_G5',false,'flag_G6a',false,'flag_G6b',false,'flag_G9',false, ...
    'zero_forcing_root',o.zero_forcing_root,'targets',tgt, ...
    'reltol',o.reltol,'abstol',o.abstol,'nmax',o.nmax, ...
    'gate_values',struct('s_tol',o.s_tol,'seam_tol',o.seam_tol,'g6a_tol',o.g6a_tol, ...
        'g6b_tol',o.g6b_tol,'g6b_refine_tol',o.g6b_refine_tol, ...
        'antisym_abs',o.antisym_abs,'antisym_rel',o.antisym_rel, ...
        'vel_abs',o.vel_abs,'vel_rel',o.vel_rel), ...
    'diag',struct());
end
