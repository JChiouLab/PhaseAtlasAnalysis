function rows = sa_t24_determinism(E)
%SA_T24_DETERMINISM  The scheduler is a pure function of
%   settings and completed outcomes. Two identical stub runs (continuation +
%   rescue, a refusal wall forcing rescue) give identical work units,
%   iterates, validity records, stores and scheduler state; the route files
%   read no clock (source scan); sa_search reads the clock only in its
%   budget check (budget_wall, live calls only).
T = 'T24 determinism';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't24');
fail = @(dd, b, ex) sa_m3_wall(dd, ex);
r1 = sa_m3_run(fullfile(d, 'a'), [3 1], 'fail', fail);
r2 = sa_m3_run(fullfile(d, 'b'), [3 1], 'fail', fail);
D = sa_m3_equiv(fullfile(d, 'a'), fullfile(d, 'b'));
rows = sa_check(rows, T, 'two identical runs: same units, iterates, validity, rows and shards', ...
    D.pass && D.n_units == r1.calls && D.n_interrupted == 0 && r1.calls > 20, ...
    sprintf('%d units; %s', D.n_units, strjoin(D.diffs, ' | ')));
rows = sa_check(rows, T, 'identical scheduler state Q (lineages, edges, rescue table)', ...
    isequaln(r1.Q, r2.Q), sprintf('%d edges', numel(r1.Q.edges)));
pkg = E.pkg; hits = {};
for f = {'sa_continuation.m', 'sa_rescue.m'}
    c = code_lines(fullfile(pkg, f{1}));
    h = c(~cellfun(@isempty, regexp(c, '\<(tic|toc|now|clock|datetime|datestr|rand|randi|randn|cputime|rng)\>', 'once')));
    hits = [hits, h]; %#ok<AGROW>
end
rows = sa_check(rows, T, 'sa_continuation / sa_rescue read no clock and no random stream', isempty(hits), strjoin(hits, ' || '));
txt = fileread(fullfile(pkg, 'sa_search.m'));
fb = regexp(txt, 'function \[ok, why\] = budget_ok.*?\nend', 'match', 'once');
rest = strrep(txt, fb, '');
c = code_lines_text(rest);
h = c(~cellfun(@isempty, regexp(c, '\<(toc|now|clock|datetime|rand|randi|randn|cputime)\>', 'once')));
rows = sa_check(rows, T, 'sa_search reads the clock only inside budget_ok (budget_wall, skipped while replaying)', ...
    ~isempty(fb) && contains(fb, 'toc(') && contains(fb, '~ctx.replay.active') && isempty(h), strjoin(h, ' || '));
rmpath(E.stub);
end

function c = code_lines(file)
c = code_lines_text(fileread(file));
end

function c = code_lines_text(t)
L = strsplit(t, newline, 'CollapseDelimiters', false);
c = L(~cellfun(@(s) startsWith(strtrim(s), '%'), L));
end
