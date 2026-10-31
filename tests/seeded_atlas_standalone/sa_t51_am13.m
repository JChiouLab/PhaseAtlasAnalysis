function rows = sa_t51_am13(E)
%SA_T51_AM13  Milestone-4 T51 (packet v0.3 AM-13): the guess-mesh defaults.
%   Model layer: 'subsample' for both models, N = 401 for fission yeast (the
%   commissioned value; 802 only as the documented pilot lever) and the S3a
%   decision for Guzzo_noXA (H1: see the S3 record); the package default stays
%   'parent' (packet v0.3 exactly, so T42 is unchanged; the stub model has no
%   model layer). T42 itself is run unchanged by the M3 suite.
T = 'T51 AM-13 guess-mesh defaults';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
n_noxa = 482;                                       % H1 S3a result under the AM-13 rule (S3 record)
Mf = sa_model('original_polarity'); Mn = sa_model('guzzo_noxa');
addpath(E.stub); cl = onCleanup(@() rmpath(E.stub));
Ms = sa_model('stub_two_pool');
[Sf, Pf] = sa_settings('resolve', struct('model', 'original_polarity', 'run_dir', 'x'), Mf);
[Sn, Pn] = sa_settings('resolve', struct('model', 'guzzo_noxa', 'run_dir', 'x'), Mn);
Ss = sa_settings('resolve', struct('model', 'stub_two_pool', 'run_dir', 'x'), Ms);
rows = sa_check(rows, T, sprintf('fission yeast subsample 401; Guzzo_noXA subsample %d; both from the model layer', n_noxa), ...
    strcmp(Sf.continuation.guess_mesh, 'subsample') && Sf.continuation.guess_nodes == 401 && ...
    strcmp(Sn.continuation.guess_mesh, 'subsample') && Sn.continuation.guess_nodes == n_noxa && ...
    strcmp(Pf.continuation.guess_mesh, 'model') && strcmp(Pn.continuation.guess_nodes, 'model'), ...
    sprintf('%s %d / %s %d', Sf.continuation.guess_mesh, Sf.continuation.guess_nodes, ...
    Sn.continuation.guess_mesh, Sn.continuation.guess_nodes));
rows = sa_check(rows, T, 'package default unchanged: guess_mesh ''parent'', guess_nodes 1000 (a model without a model layer)', ...
    strcmp(Ss.continuation.guess_mesh, 'parent') && Ss.continuation.guess_nodes == 1000);
U = struct('model', 'original_polarity', 'run_dir', 'x', 'continuation', struct('guess_nodes', 802));
Su = sa_settings('resolve', U, Mf);
U2 = struct('model', 'guzzo_noxa', 'run_dir', 'x', 'continuation', struct('guess_nodes', 241));
Su2 = sa_settings('resolve', U2, Mn);
rows = sa_check(rows, T, 'the user may set N (fission yeast 802 = the documented pilot lever; noxa 241 = the M3 commissioned value)', ...
    Su.continuation.guess_nodes == 802 && Su2.continuation.guess_nodes == 241);
end
