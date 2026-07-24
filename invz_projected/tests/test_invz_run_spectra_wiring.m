function tests = test_invz_run_spectra_wiring
% Step-5 Task 8: static source-contract test for the invz_run_spectra opt-in Ewald wiring.
% The driver script is interactive and its default parameters are expensive, so this test does
% NOT execute it -- it asserts the wiring CONTRACT on the source text. The callable backend
% behaviour is covered by test_invz_spectra_forward_ewald and the Gate-C7 end-to-end matrix.
tests = functiontests(localfunctions);
end

function src = read_src()
here = fileparts(mfilename('fullpath'));
src  = fileread(fullfile(here, '..', 'invz_run_spectra.m'));
end

function test_default_backend_is_bruteforce(testCase)
% The default driver must reproduce the legacy bruteforce path (Step 7 flips the production
% default, not this script): the knob defaults to bruteforce and is never hardcoded to ewald.
src = read_src();
verifyNotEmpty(testCase, regexp(src, "dipoleBackend\s*=\s*'bruteforce'\s*;", 'once'), ...
    'run_spectra must default dipoleBackend to bruteforce.');
verifyEmpty(testCase, regexp(src, "dipoleBackend\s*=\s*'ewald'\s*;", 'once'), ...
    'run_spectra default backend must not be hardcoded to ewald.');
end

function test_backend_knobs_declared(testCase)
src = read_src();
for kn = {'dipoleBackend', 'ewaldOpts', 'gridConvention', 'gridOffset', 'gammaPolicy'}
    verifyNotEmpty(testCase, regexp(src, [kn{1} '\s*='], 'once'), ...
        sprintf('run_spectra is missing the backend knob %s.', kn{1}));
end
end

function test_coupling_opts_built_and_conditional(testCase)
% coupling_opts is built from the knob; Ewald controls only for ewald; grid-policy fields only
% when their knob is nonempty -- so the default script activates NO new grid route.
src = read_src();
verifyNotEmpty(testCase, regexp(src, "coupling_opts\s*=\s*struct\('dipole',\s*dipoleBackend\)", 'once'), ...
    'coupling_opts must be built from the dipoleBackend knob.');
verifyNotEmpty(testCase, regexp(src, "if\s+strcmp\(dipoleBackend,\s*'ewald'\)[^\n]*coupling_opts\.ewald", 'once'), ...
    'Ewald controls must be attached only when the ewald backend is selected.');
verifyNotEmpty(testCase, regexp(src, "if\s*~isempty\(gridConvention\)[^\n]*coupling_opts\.gridConvention", 'once'), ...
    'gridConvention must be forwarded only when nonempty.');
verifyNotEmpty(testCase, regexp(src, "if\s*~isempty\(gridOffset\)[^\n]*coupling_opts\.gridOffset", 'once'), ...
    'gridOffset must be forwarded only when nonempty.');
verifyNotEmpty(testCase, regexp(src, "if\s*~isempty\(gammaPolicy\)[^\n]*coupling_opts\.gammaPolicy", 'once'), ...
    'gammaPolicy must be forwarded only when nonempty.');
end

function test_all_three_calls_use_merged_opts(testCase)
% Every spectra_map / spectra_qpath call must receive the merged coupling opts; the pre-Task-8
% un-merged inline opts struct must be gone from the call sites.
src = read_src();
% anchor on the final argument + closing paren (robust against inner parens like Bq(k)*dhat)
verifyEqual(testCase, numel(regexp(src, "wqMeV,\s*qpOpts\)")), 2, ...
    'both q-path calls (scalar + loop) must pass the merged qpOpts.');
verifyNotEmpty(testCase, regexp(src, "wMeV,\s*mapOpts\)", 'once'), ...
    'the field-map call must pass the merged mapOpts.');
verifyEmpty(testCase, regexp(src, "wqMeV,\s*struct\('eta',\s*eta,\s*'solve_opts',\s*solve_opts\)", 'once'), ...
    'q-path calls must no longer use the un-merged inline opts struct.');
end
