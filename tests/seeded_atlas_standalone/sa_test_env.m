function E = sa_test_env()
%SA_TEST_ENV  Paths and reference inputs shared by the zero-solve tests.
%   Read-only references only; every test output goes under E.scratch.
%   E.exp, E.campaign_matlab, E.raw, E.lattice_csv, E.inherit, E.lib14 and E.lib79
%   point to the authors' development data (fission-yeast libraries at L = 14 and
%   L = 7.9, stored solver returns), which is not shipped. Tests that need it SKIP
%   without it; the portable tests do not use it.

here = fileparts(mfilename('fullpath'));
E = struct();
E.root = fileparts(fileparts(here));
E.tests = here;
E.pkg = fullfile(E.root, 'src', 'pipelines', 'seeded_atlas_standalone');
E.module = fullfile(E.pkg, 'models', 'original_polarity');
E.stub = fullfile(here, 'fixtures', 'stub_two_pool');
E.exp = fullfile(E.root, 'data', 'experiments', 'coupled_half_domain_atlas_20260811');
E.campaign_matlab = fullfile(E.exp, 'matlab');
E.raw = fullfile(E.exp, 'outputs', 'sca_prb', 'raw');
E.lattice_csv = fullfile(E.exp, 'outputs', 'sca_n0b', 'sca_n0b_lattice_20260904_N0B2.csv');
E.inherit = fullfile(E.exp, 'outputs', 'sca_drv', 'configs', ...
    'inherit_l14_stage_b_frozen_20260906_SCAPR0B.mat');      % read only
E.m1_root_store = fullfile(E.root, 'src', 'pipelines', 'seeded_coupled_atlas', ...
    'seeded_atlas_root_store.m');
E.lib14 = fullfile(E.root, 'data', ...
    'a0p5945_b0p1852_c75p785_d0p552_e0p8278_f1_Dm0p05_Dc5_L14_MBG28_MP6p89388_res95', 'outputs');
E.lib79 = fullfile(E.root, 'data', ...
    'a0p5945_b0p1852_c75p785_d0p552_e0p8278_f1_Dm0p05_Dc5_L7p9_MBG15p8_MP3p89012_res95', 'outputs');
E.scratch = fullfile(E.root, 'data', 'experiments', 'standalone_seeded_atlas_20261012');

% parameter sets of the authors' development runs, not the manuscript's Table 2 setting (L = 7)
base = struct('a', 0.5945, 'b', 0.1852, 'c', 75.785, 'd', 0.552, 'e', 0.8278, ...
    'f', 1, 'Dm', 0.05, 'Dc', 5);
E.p14 = base; E.p14.L = 14;  E.p14.M_BG = 28;   E.p14.M_P = 6.89388;
E.p79 = base; E.p79.L = 7.9; E.p79.M_BG = 15.8; E.p79.M_P = 3.890118;   % library par (folder name rounds it)

% stored solver returns (development data) used by the re-diagnosis and replay tests (read only)
E.raw_accepted = 'SCA_PRB_W1B_RAW_call1093_HSS_HSS_cold_i_0_j_0.mat';
E.raw_g2       = 'SCA_PRB_CLOSUREB1_RAW_call0167_HSS_MPSS_subdiv1_s2_i_3_j_28_s2.mat';
E.raw_hpss     = 'SCA_PRB_W1B_RAW_call5668_HPSS_HPSS_cold_i_0_j_0.mat';
end
