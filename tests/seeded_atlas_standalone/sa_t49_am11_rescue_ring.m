function rows = sa_t49_am11_rescue_ring(E)
%SA_T49_AM11_RESCUE_RING  Rescue is on in every mode and
%   rescue.ring defaults to true; rescue.ring = false gives continuation only in
%   the ring band while the interior is unchanged; with rescue.ring = true the
%   ring behaviour is that of sa_t40_ring, which runs under it
%   unchanged.
T = 'T49 AM-11 rescue scope';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
M = sa_model('stub_two_pool');
ok = true; det = {};
for m = sa_modes('modes')
    u = struct('model', 'stub_two_pool', 'run_dir', 'x', 'mode', m{1});
    if any(strcmp(m{1}, {'production', 'focus'})), u.budget = struct('max_calls', 10, 'max_wall_s', 10); end
    if strcmp(m{1}, 'focus'), u.import = struct('source_dir', 'x'); end   % focus mode needs a source
    [S, P] = sa_settings('resolve', u, M);
    % every mode layer sets rescue.ring
    % false and rescue.rungs {secant}; the package default keeps ring true and
    % the full ladder
    ok = ok && S.routes.rescue && ~S.rescue.ring && S.routes.continuation && strcmp(P.routes.rescue, 'mode') && ...
        strcmp(P.rescue.ring, 'mode') && isequal(S.rescue.rungs, {'secant'});
    det{end+1} = sprintf('%s rescue %d ring %d rungs {%s}', m{1}, S.routes.rescue, S.rescue.ring, ...
        strjoin(S.rescue.rungs, ',')); %#ok<AGROW>
end
D0 = sa_settings('defaults', M);
rows = sa_check(rows, T, ['defaults: rescue on, rescue.ring = false and rungs {secant} (mode) in smoke, ' ...
    'pilot and production; package default ring true, ladder diagonal-secant-deep'], ...
    ok && D0.rescue.ring && isequal(D0.rescue.rungs, {'diagonal', 'secant', 'deep'}), strjoin(det, '; '));
d = sa_m3_testdir(E, 't49');
fail = @(dd, b, ex) ringwall(dd, ex);
rT = sa_m3_run(fullfile(d, 'ring_true'), [], 'fail', fail);
rF = sa_m3_run(fullfile(d, 'ring_false'), [], 'fail', fail, 'settings', struct('rescue', struct('ring', false)));
cnt = @(r, kind) sum(strcmp({r.Q.edges.kind}, kind) & r.Q.ring([r.Q.edges.q]).');
resc = @(r) cnt(r, 'rescue_diagonal') + cnt(r, 'rescue_secant') + cnt(r, 'rescue_deep');
iT = rT.attempts(~strcmp({rT.attempts.phase}, 'ring')); iF = rF.attempts(~strcmp({rF.attempts.phase}, 'ring'));
same_int = numel(iT) == numel(iF) && isequal({iT.unit_key}, {iF.unit_key});
rows = sa_check(rows, T, ['rescue.ring = false: no rescue edge into a ring node, continuation edges into the ring remain, ' ...
    'interior work identical; rescue.ring = true rescues in the ring'], resc(rT) > 0 && resc(rF) == 0 && ...
    cnt(rF, 'continuation') > 0 && same_int && sum(rF.Q.rescue_done & rF.Q.ring) == 0, ...
    sprintf('ring rescue edges %d / %d; ring continuation edges %d / %d', resc(rT), resc(rF), ...
    cnt(rT, 'continuation'), cnt(rF, 'continuation')));
rmpath(E.stub);
r40 = sa_t40_ring(E);
rows = sa_check(rows, T, 'T40 (ring, v0.4 behaviour) passes under the rescue.ring = true default', all([r40.pass]), ...
    sprintf('%d of %d', sum([r40.pass]), numel(r40)));
end

function c = ringwall(dd, ex)
% ring wall: ring nodes with k1 = 6 and |k2| <= 1 refuse every non-exact guess
c = 'S';
if abs(dd(1) - 3.9) < 1e-9 && abs(dd(2)) < 0.7 && ~ex, c = 'R'; end
end
