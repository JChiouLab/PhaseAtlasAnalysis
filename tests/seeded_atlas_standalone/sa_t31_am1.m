function rows = sa_t31_am1(E)
%SA_T31_AM1  A record written under the
%   current software version, MATLAB version and thread mode is re-derived
%   and must equal its stored checks: a tampered value, tolerance,
%   comparison or opts is an integrity failure. An older-version record
%   (no or different provenance) is re-derived under the current rules and
%   reported BESIDE its stored verdict: no integrity failure, no relabel.
T = 'T31 AM-1 re-derivation';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't31');
res = sa_m3_run(fullfile(d, 'a'), [1 0]);
M = sa_model('stub_two_pool');
Ef = load(fullfile(sa_layout('state', fullfile(d, 'a')), 'settings_effective.mat'));
S = Ef.settings; th = Ef.theta;
rec = load(fullfile(sa_layout('state', fullfile(d, 'a')), 'attempts', 'A000002.mat'));
prov = sa_settings('provenance');
ctx = struct('M', M, 'theta', th, 'S', S);
o = sa_search('opts_for', ctx, sa_search('solver_for', S, rec.unit.kind));
[V, info] = sa_validity('rederive', M, rec, th, S, o, prov);
rows = sa_check(rows, T, 'untampered same-version record: integrity ok, verdict and checks reproduced', ...
    info.integrity_ok && info.same_version && isequaln(V, rec.validity) && V.valid, info.reason);
k = find(strcmp({rec.validity.checks.name}, 'stub_noflux'));
tam = {'value', 5e-9; 'tolerance', 1; 'comparison', '<='};
for t = 1:size(tam, 1)
    r2 = rec; r2.validity.checks(k).(tam{t, 1}) = tam{t, 2};
    [~, i2] = sa_validity('rederive', M, r2, th, S, o, prov);
    rows = sa_check(rows, T, sprintf('tampered %s -> integrity failure', tam{t, 1}), ...
        ~i2.integrity_ok && contains(i2.reason, ['stub_noflux.' tam{t, 1}]), i2.reason);
end
r2 = rec; r2.opts.reltol = 1e-5;
[~, i2] = sa_validity('rederive', M, r2, th, S, o, prov);
rows = sa_check(rows, T, 'tampered opts -> integrity failure (opts rebuilt from the effective settings)', ...
    ~i2.integrity_ok && contains(i2.reason, 'opts'), i2.reason);
r3 = rmfield(rec, 'provenance'); r3.validity.checks(k).value = 1; r3.validity.valid = false;
[V3, i3] = sa_validity('rederive', M, r3, th, S, o, prov);
rows = sa_check(rows, T, 'older-version record (no provenance), stored value changed: no integrity failure; recomputed verdict beside the stored one', ...
    i3.integrity_ok && ~i3.same_version && i3.valid_stored == 0 && i3.valid_recomputed == 1 && i3.flip && V3.valid, ...
    sprintf('stored %d recomputed %d', i3.valid_stored, i3.valid_recomputed));
r4 = rec; r4.provenance.matlab_version = 'older';
[~, i4] = sa_validity('rederive', M, r4, th, S, o, prov);
rows = sa_check(rows, T, 'record from another MATLAB version: treated as older (side by side)', ...
    i4.integrity_ok && ~i4.same_version, '');
rmpath(E.stub);
end
