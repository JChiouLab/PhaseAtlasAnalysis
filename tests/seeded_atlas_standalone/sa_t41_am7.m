function rows = sa_t41_am7(E)
%SA_T41_AM7  identity.candidate_mesh =
%   'comparison' stores candidate profiles evaluated by deval on the
%   comparison mesh (and computes pair distances there). Against the default
%   'dense' run of the same stub search: identical memberships and outcomes;
%   comparison profiles and representative distances agree to the retired
%   interpolation error; the stored returned profile is ~10x smaller.
T = 'T41 AM-7 comparison mesh';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't41');
rD = sa_m3_mc_run(fullfile(d, 'dense'));
rC = sa_m3_mc_run(fullfile(d, 'comparison'), struct('identity', struct('candidate_mesh', 'comparison')));
rows = sa_check(rows, T, 'same work units and outcomes in both modes', ...
    isequal({rD.attempts.unit_key}, {rC.attempts.unit_key}) && isequal({rD.attempts.outcome}, {rC.attempts.outcome}) && ...
    isequal({rD.attempts.root_id}, {rC.attempts.root_id}), sprintf('%d calls', rD.calls));
mx = 0; dmx = 0; same = true; nd = 0; nc = 0; nn = 0;
FD = sa_store_io('full_map', rD); FC = sa_store_io('full_map', rC);   % result.shards are light
for k = FD.keys
    a = FD(k{1}); b = FC(k{1});
    same = same && numel(a.roots) == numel(b.roots);
    for r = 1:numel(a.roots)
        mx = max(mx, max(abs(a.roots(r).representative_profile - b.roots(r).representative_profile), [], 'all'));
        da = [a.roots(r).discoveries.representative_distance]; db = [b.roots(r).discoveries.representative_distance];
        dmx = max(dmx, max(abs(da - db)));
        same = same && isequal({a.roots(r).discoveries.attempt_id}, {b.roots(r).discoveries.attempt_id});
        nd = nd + numel(a.roots(r).discoveries(1).returned_mesh); nc = nc + numel(b.roots(r).discoveries(1).returned_mesh);
        nn = nn + 1;
    end
end
rows = sa_check(rows, T, 'same memberships; comparison profiles and distances agree to <= 1e-10', ...
    same && mx <= 1e-10 && dmx <= 1e-10, sprintf('max profile diff %.3g, max distance diff %.3g', mx, dmx));
Z = load(fullfile(sa_layout('state', fullfile(d, 'comparison')), 'settings_effective.mat'));
rows = sa_check(rows, T, 'stored candidate mesh = the comparison mesh (2001 nodes) instead of n_dense (20001)', ...
    nc / nn == 2001 && nd / nn == 20001 && strcmp(Z.settings.identity.candidate_mesh, 'comparison'), ...
    sprintf('nodes per root %g vs %g (%d roots)', nc / nn, nd / nn, nn));
M = sa_model('stub_two_pool');
s1 = sa_m3_stub('sol', sa_m3_stub('params'), [2 2 2 2], 'A'); s2 = sa_m3_stub('sol', sa_m3_stub('params'), [2 2 2 2], 'B');
S = Z.settings; K = Z.store_settings;
Pc = sa_identity('pair', M, s1, s2, K, S);
S.identity.candidate_mesh = 'dense';
Pd = sa_identity('pair', M, s1, s2, K, S);
rows = sa_check(rows, T, 'pair distance on the comparison mesh equals the dense-grid interpolation to <= 3.3e-6 relative', ...
    abs(Pc.scaled_rms - Pd.scaled_rms) <= 3.3e-6 * Pd.scaled_rms, sprintf('%.12g vs %.12g', Pc.scaled_rms, Pd.scaled_rms));
rmpath(E.stub);
end
