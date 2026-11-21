function out = sa_mesh_trace(action, varargin)
%SA_MESH_TRACE  Per-call recorder of the mesh sizes bvp4c proposes.
%
%   t = SA_MESH_TRACE('reset')
%       A new call token and an empty, armed store. The adapters (sap_solve,
%       sag_solve) call it right before the coupled-BVP solve.
%   SA_MESH_TRACE('initial', N)
%       Record the initial mesh N0 (sa_bvp4c_traced, after N = length(x)).
%   SA_MESH_TRACE('note', N)
%       Append one proposal: the count returned by new_profile
%       (sa_bvp4c_traced, after the new_profile line).
%       'initial' and 'note' NEVER error, NEVER warn, never touch lastwarn and
%       need no prior 'reset'. Misuse is caught at 'take' (a take without a
%       reset, or with another token, is SA:TRACE:E_STALE) and by the
%       consistency checks of sa_record_call. Only sa_bvp4c_traced and test
%       stubs call them.
%   tr = SA_MESH_TRACE('take', t, nmax)
%       The recorded trace of the call with token t, and clears/disarms the
%       store:
%         initial      N0, or NaN when not recorded
%         proposals    row of the proposals, in order
%         n_proposals  number of proposals
%         max_proposed P = max(proposals); 0 when n_proposals = 0 and N0 was
%                      recorded (stub solver, or a throw before N0)
%         over_nmax    the last proposal > nmax (the call stopped on a
%                      proposal above its own NMax)
%         n_initial    how many 'initial' notes arrived (1 for one bvp4c run)
%         recorder     'sa_bvp4c_traced/1 (bvp4c <release> <sha12>)': the
%                      release (['R' version('-release')]) and the toolbox
%                      bvp4c.m SHA-256 prefix of the generation record,
%                      cached at run start by sa_bvp4c_setup('component')
%                      (never read from a file per call); on the reference
%                      toolbox 'sa_bvp4c_traced/1 (bvp4c R2025b 597839cef613)'
%         token        t
%       A stale or mismatched token, or a take without a reset, is
%       SA:TRACE:E_STALE (raised by the adapter, outside the solver's own
%       try/catch, so it is never turned into a no-return).
%   tr = SA_MESH_TRACE('blank')    the "not recorded" trace (NaN), e.g. for
%       a record written without a trace
%   P = SA_MESH_TRACE('max_of', tr)   P recomputed from a stored trace's
%       recorded or tr is not a trace). Replay uses it.

persistent S
if isempty(S)
    S = empty_store(0, false);
end
switch action
    case 'note'
        % never error, never warn: everything inside a bare try/catch
        try
            S.proposals(end+1) = varargin{1};
            if ~S.armed, S.stray = S.stray + 1; end
        catch
        end
    case 'initial'
        try
            S.initial = varargin{1};
            S.n_initial = S.n_initial + 1;
            if ~S.armed, S.stray = S.stray + 1; end
        catch
        end
    case 'reset'
        S = empty_store(S.token + 1, true);
        out = S.token;
    case 'take'
        t = varargin{1};
        nmax = Inf; if numel(varargin) >= 2, nmax = varargin{2}; end
        if ~S.armed || ~isequal(t, S.token)
            S = empty_store(S.token, false);
            error('SA:TRACE:E_STALE', ['mesh trace taken with a stale or mismatched token (or without a ' ...
                'reset): the recorder and the adapter are out of step.']);
        end
        out = make_trace(S.initial, S.proposals, S.n_initial, nmax, S.token);
        S = empty_store(S.token, false);
    case 'blank'
        out = make_trace(NaN, zeros(1, 0), 0, Inf, NaN);
    case 'max_of'
        out = max_of(varargin{1});
    otherwise
        error('SA:TRACE', 'Unknown sa_mesh_trace action %s.', action);
end
end

function S = empty_store(token, armed)
S = struct('token', token, 'armed', armed, 'initial', NaN, 'proposals', zeros(1, 0), ...
    'n_initial', 0, 'stray', 0);
end

function tr = make_trace(initial, proposals, n_initial, nmax, token)
k = numel(proposals);
if k > 0
    P = max(proposals);
    over = proposals(end) > nmax;
elseif isfinite(initial)
    P = 0; over = false;
else
    P = NaN; over = false;
end
tr = struct('initial', double(initial), 'proposals', double(proposals(:).'), 'n_proposals', k, ...
    'max_proposed', double(P), 'over_nmax', logical(over), 'n_initial', n_initial, ...
    'recorder', sa_bvp4c_setup('recorder'), 'token', token);
end

function P = max_of(tr)
P = NaN;
if ~(isstruct(tr) && isscalar(tr) && isfield(tr, 'proposals') && isfield(tr, 'initial')), return; end
if ~isempty(tr.proposals)
    P = max(double(tr.proposals));
elseif isfinite(tr.initial)
    P = 0;
end
end
