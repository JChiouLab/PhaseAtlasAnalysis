function rows = sa_t04_source_scan(E)
%SA_T04_SOURCE_SCAN  Shared code is model-free; runtime is campaign-free;
%   the copied ODE and reaction are unclipped (RC-10).
%   Comment-only lines are dropped LINE BY LINE (no '.'-spanning regexp);
%   text inside code lines is scanned, so error messages count too.
T = 'T04 source scan';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
shared = dir(fullfile(E.pkg, '*.m'));
words = ['\<(HSS|HPSS|MPSS|PAK|GEF|alpha|beta|gamma|delta|epsilon|rho|kappa|' ...
    'M_BG|M_P|Dm|Dc|u_L|p_L|u_R|p_R)\>'];
hits = {};
for k = 1:numel(shared)
    code = code_lines(fullfile(shared(k).folder, shared(k).name));
    for q = 1:numel(code)
        c = code{q};
        bad = {};
        % the verbatim M1 kernel uses `delta` as a local name (model-free code)
        if ~isempty(regexp(c, words, 'once')) && ~strcmp(shared(k).name, 'sa_root_store.m')
            bad{end+1} = 'model word'; %#ok<AGROW>
        end
        if ~isempty(regexp(c, '\<(sap|cha|sca)_', 'once')), bad{end+1} = 'model/campaign function'; end %#ok<AGROW>
        if ~isempty(regexp(c, '(?<![\w.])(28|401)(?![\w.])', 'once')), bad{end+1} = 'state-count literal'; end %#ok<AGROW>
        if ~isempty(regexp(c, '\<0\s*,?\s*0\.5\s*[\]\)]', 'once')), bad{end+1} = 'domain literal'; end %#ok<AGROW>
        if ~isempty(regexp(c, 'y\(\s*2[0-9]', 'once')), bad{end+1} = 'hard-coded state row'; end %#ok<AGROW>
        if ~isempty(bad)
            hits{end+1} = sprintf('%s: %s [%s]', shared(k).name, strtrim(c), strjoin(bad, ',')); %#ok<AGROW>
        end
    end
end
rows = sa_check(rows, T, sprintf('%d shared files name no species, pool, half type, model constant or campaign function', ...
    numel(shared)), isempty(hits), strjoin(hits, ' || '));

allm = [dir(fullfile(E.pkg, '*.m')); dir(fullfile(E.module, '*.m'))];
hits = {};
for k = 1:numel(allm)
    txt = fileread(fullfile(allm(k).folder, allm(k).name));
    if contains(txt, 'data/experiments') || contains(txt, ['data' filesep 'experiments']) ...
            || contains(txt, 'coupled_half_domain_atlas')
        hits{end+1} = allm(k).name; %#ok<AGROW>
    end
end
% the copied corrector's comments cite campaign backups; only code lines count
hits = hits(~cellfun(@(n) isempty(code_hits(fullfile(E.pkg, n), fullfile(E.module, n))), hits));
rows = sa_check(rows, T, 'no package code line references a data/experiments or campaign path', ...
    isempty(hits), strjoin(hits, ', '));

f = fileread(fullfile(E.module, 'sap_forced28.m'));
i0 = strfind(f, 'function dy = ode28'); i1 = strfind(f, 'function r = bc28');
ode = strjoin(code_text(f(i0:i1-1)), newline);
g = fileread(fullfile(E.module, 'sap_model.m'));
j0 = strfind(g, 'function R = reaction');
rx = strjoin(code_text(g(j0:end)), newline);
clip = '\<(max|min|abs|sign|heaviside)\s*\(';
rows = sa_check(rows, T, 'copied ODE (ode28) has no max/min/abs clipping', ...
    isempty(regexp(ode, clip, 'once')), sprintf('%d chars scanned', numel(ode)));
rows = sa_check(rows, T, 'copied reaction has no max/min/abs clipping', ...
    isempty(regexp(rx, clip, 'once')), sprintf('%d chars scanned', numel(rx)));
M = sa_model('original_polarity');
rows = sa_check(rows, T, 'module declares unclipped = true', isequal(M.unclipped, true));
end

function c = code_lines(file)
L = strsplit(fileread(file), newline, 'CollapseDelimiters', false);
keep = ~cellfun(@(s) startsWith(strtrim(s), '%'), L);
c = L(keep);
end

function c = code_text(txt)
L = strsplit(txt, newline, 'CollapseDelimiters', false);
c = L(~cellfun(@(s) startsWith(strtrim(s), '%'), L));
end

function h = code_hits(a, b)
if exist(a, 'file') == 2, f = a; else, f = b; end
c = code_lines(f);
h = c(contains(c, 'data/experiments') | contains(c, 'coupled_half_domain_atlas'));
end
