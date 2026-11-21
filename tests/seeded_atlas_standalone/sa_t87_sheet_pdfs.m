function rows = sa_t87_sheet_pdfs(E)
%SA_T87_SHEET_PDFS  The 2-D sheets as PDF files (post/sa_plot_sheets.py).
%   Zero bvp4c calls (stub model, the multi-seed run of sa_t28_lineages).
%   Without Python (opts.python names a file that does not exist): no
%   error, sheets_2d_roots.csv (one row per root) and sheets_2d.json are
%   written, no PDF and no key, out.sheets empty. With the system's python3
%   (only when it imports numpy, matplotlib and contourpy; otherwise those
%   assertions are reported as skipped): one PDF per observed family plus
%   all_families.pdf, no .fig sheet, and a second call redraws nothing.
T = 'T87 sheet PDFs';
rows = struct('test', {}, 'assertion', {}, 'pass', {}, 'detail', {});
addpath(E.stub);
d = sa_m3_testdir(E, 't87');
run_dir = fullfile(d, 'U');
sa_m3_mc_run(run_dir, struct());
n = numel(seeded_atlas_roots(run_dir, struct('verify', false)));

% ---- no Python: the inputs only -------------------------------------------------------
o1 = fullfile(d, 'nopy');
out = seeded_atlas_plot(run_dir, struct('out', o1, 'python', fullfile(d, 'no_such_python')));
sd = fullfile(o1, 'sheets_2d');
L = splitlines(strtrim(fileread(fullfile(sd, 'sheets_2d_roots.csv'))));
rows = sa_check(rows, T, 'without Python: sheets_2d_roots.csv has a header and one row per root', ...
    numel(L) == n + 1 && strcmp(L{1}, 'root_id,family,ix,iy,x,y,dn2,vx,vy'), sprintf('%d rows, %d roots', numel(L) - 1, n));
J = jsondecode(fileread(fullfile(sd, 'sheets_2d.json')));
rows = sa_check(rows, T, 'without Python: sheets_2d.json lists the families and the lattice', ...
    isfield(J, 'order') && isfield(J, 'step') && numel(J.step) == 2, '');
rows = sa_check(rows, T, 'without Python: no PDF, no key, out.sheets empty, sheets_redrawn false', ...
    isempty(dir(fullfile(sd, '*.pdf'))) && exist(fullfile(sd, 'sheets_2d_key.mat'), 'file') ~= 2 && ...
    isempty(out.sheets) && ~out.sheets_redrawn, '');

% ---- with Python ------------------------------------------------------------------------
[st, ~] = system('python3 -c "import numpy, matplotlib, contourpy"');
if st ~= 0
    rows = sa_check(rows, T, 'with Python: skipped (python3 with numpy, matplotlib, contourpy not found)', true, 'SKIPPED');
    rmpath(E.stub);
    return
end
o2 = fullfile(d, 'py');
out = seeded_atlas_plot(run_dir, struct('out', o2, 'python', 'python3'));
sd = fullfile(o2, 'sheets_2d');
J = jsondecode(fileread(fullfile(sd, 'sheets_2d.json')));
fams = cellstr(J.order);
want = [cellfun(@(f) [sa_mode_name('file_tag', f) '.pdf'], fams(:).', 'UniformOutput', false), {'all_families.pdf'}];
have = {dir(fullfile(sd, '*.pdf')).name};
rows = sa_check(rows, T, 'with Python: one PDF per family plus all_families.pdf, and no .fig sheet', ...
    isempty(setxor(want, have)) && isempty(dir(fullfile(sd, '*.fig'))) && out.sheets_redrawn && ...
    numel(out.sheets) == numel(want), strjoin(have, ' '));
out2 = seeded_atlas_plot(run_dir, struct('out', o2, 'python', 'python3'));
rows = sa_check(rows, T, 'with Python: a second call redraws nothing (keys equal)', ...
    ~out2.redrawn && ~out2.sheets_redrawn && numel(out2.sheets) == numel(want), '');
rmpath(E.stub);
end
