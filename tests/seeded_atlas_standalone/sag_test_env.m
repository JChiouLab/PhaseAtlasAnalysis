function E = sag_test_env(scratch)
%SAG_TEST_ENV  Paths and reference inputs of the Guzzo_noXA zero-solve tests.
%   Read-only references only; every test output goes under E.scratch.
%   E.myxo_code, E.myxo_sheets, E.libs, E.lib14 and E.sheet_dirs point to the authors'
%   development data (reference M. xanthus sheets, uncoupled libraries and
%   stored runs), which is not shipped. Tests that need it SKIP without it.
here = fileparts(mfilename('fullpath'));
E = struct();
E.root = fileparts(fileparts(here));
E.tests = here;
E.pkg = fullfile(E.root, 'src', 'pipelines', 'seeded_atlas_standalone');
E.module = fullfile(E.pkg, 'models', 'guzzo_noxa');
E.core = fullfile(E.root, 'src', 'models', 'guzzo_noxa');
E.solver_noxa = fullfile(E.root, 'src', 'solvers', 'atlas_surface_pipeline_dependencies', ...
    'models', 'guzzo_noxa');
E.original_module = fullfile(E.pkg, 'models', 'original_polarity');
E.myxo_code = fullfile(E.root, 'myxo', 'codes', 'atlas_noxa');
E.myxo_sheets = fullfile(E.root, 'myxo', 'simulations', 'atlas_noxa', 's0p35', 'l6');
E.libs = fullfile(E.root, 'data', 'experiments', 'standalone_seeded_atlas_20261012', 'm2m', ...
    'phaseB2_20261012_234352', 'libraries');
E.lib14 = fullfile(E.root, 'data', ...
    'a0p5945_b0p1852_c75p785_d0p552_e0p8278_f1_Dm0p05_Dc5_L14_MBG28_MP6p89388_res95', 'outputs');
E.model_version = '04a0838273bebf2a';             % model_version prefix of the shipped gnx_* code
E.signals = [0 0.02 0.05 0.07 0.35 1];
E.lib_dirs = {'guzzo_noxa_L6_s0p35_dba0_rR1_rB1_h1p8_res25_c07275cf', ...
    'guzzo_noxa_L6_s0p35_dba0p02_rR1_rB1_h1p8_res25_2903fe00', ...
    'guzzo_noxa_L6_s0p35_dba0p05_rR1_rB1_h1p8_res25_8cac20b8', ...
    'guzzo_noxa_L6_s0p35_dba0p07_rR1_rB1_h1p8_res25_e9f0298a', ...
    'guzzo_noxa_L6_s0p35_dba0p35_rR1_rB1_h1p8_res25_b8820988', ...
    'guzzo_noxa_L6_s0p35_dba1_rR1_rB1_h1p8_res25_e0dd357e'};
E.sheet_dirs = {'stage3c_coupled_dba0', 'stage3c_coupled_dba0p02', 'stage3c_coupled_dba0p05', ...
    'stage3c_coupled_dba0p07', 'stage3c_coupled_dba0p35', 'stage3c_coupled'};
if nargin < 1 || isempty(scratch)
    scratch = fullfile(E.root, 'data', 'experiments', 'standalone_seeded_atlas_20261012', 'm2c', 'tests_scratch');
end
E.scratch = scratch;
end
