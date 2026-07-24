function tests = test_invz_anchor_couplings
%TEST_INVZ_ANCHOR_COUPLINGS Ewald Step-5 Task 5: parity + injection tests for the thin,
% presence-preserving invz_anchor_couplings helper (docs/superpowers/plans/2026-07-24-ewald-
% step5-integration.md Task 5; authority docs/invzp_ewald_prereg.md Sec.7 "Frozen Gate-E",
% docs/invzp_ewald_design.md item 6, docs/invzp_ewald_integration_map.md Sec.5.4(B)).
%
% Covers: (1) default-call bit-for-bit parity against BOTH an inline legacy qVec_generator +
% Gamma-drop + invz_jq_modes construction AND a direct invz_bz_couplings call at the same
% explicit dpRng/cache (proving the helper adds nothing beyond dpRng/cache defaulting); (2)
% dpRng-default-applies-only-when-absent presence semantics; (3) an Ewald + halfopen + [0 0 0]
% + P_complete injection call verifying COMPLETE backend and grid provenance; (4) a second,
% mechanical P_drop call verifying the retained row count differs from P_complete as expected;
% (5) that opts.dipole/opts.ewald alone (no grid-policy field), forwarded through the helper,
% does not activate the new grid route -- exactly mirroring invz_bz_couplings' own Task-4
% routing rule, since the helper adds no routing logic of its own.
%
% Mechanical provenance ONLY throughout: no Sigma_c, Tc, Bc, Jensen field, or anchor inequality is
% computed anywhere in this file (Step 7's job, out of this task's scope).
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                        % invz_projected: invz_bz_couplings, invz_jq_modes, invz_phase1_qgrid (invz_anchor_couplings itself lives right here in tests/, already on path when this file runs)
addpath(fullfile(here, '..', '..'));                  % repo root: qVec_generator, invz_dipole_ewald, MF_dipole, exchange
addpath(fullfile(here, '..', '..', 'invz_common'));   % invz_ion, invz_const, getf
end

% =====================================================================
% shared helpers (names deliberately free of "test" so functiontests skips them)
% =====================================================================
function a0 = alpha0_of(a)
a0 = sqrt(pi)/abs(det(a))^(1/3);
end

function eo = default_eopts(a)
% Frozen production Ewald defaults (docs/invzp_ewald_prereg.md Sec.2), built from the live
% alpha0_of -- never a rounded displayed constant.
a0 = alpha0_of(a);
eo = struct('alpha', a0, 'r_cut', 5.5/a0, 'g_cut', 11*a0, 'boundary', 'conducting_k0_omitted');
end

function qc = legacy_qc(ion, grid)
% Independent, byte-for-byte reproduction of invz_bz_couplings.m's own absent-grid-field
% qVec_generator + Gamma-drop construction -- the actual inline legacy construction the helper's
% default path must keep reproducing exactly.
[qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5], 'verbose', false);
qc = qc(any(abs(qc) > 1e-12, 2), :);
end

% =====================================================================
% Parity: default helper call == inline legacy construction, bit-for-bit
% =====================================================================
function test_default_helper_matches_inline_legacy_construction(testCase)
ion = invz_ion();
grid = [4 4 4];

qcRef = legacy_qc(ion, grid);
[JnuRef, infoRef] = invz_jq_modes(ion, qcRef, struct('dpRng', 30, 'cache', false));
JnuRef = JnuRef(:);
Jaa0Ref = ion.Jxx0; if isfield(infoRef, 'Jaa0'), Jaa0Ref = infoRef.Jaa0; end

[JnuHelp, infoHelp, Jaa0Help] = invz_anchor_couplings(ion, struct('grid', grid, 'cache', false));

verifyEqual(testCase, JnuHelp, JnuRef);
verifyEqual(testCase, infoHelp, infoRef);
verifyEqual(testCase, Jaa0Help, Jaa0Ref);
verifyFalse(testCase, isfield(infoHelp, 'grid'), 'the default helper route must keep info.grid ABSENT.');
end

function test_default_helper_matches_direct_bz_couplings_call(testCase)
% Second, independent parity anchor at a different cheap grid size: the helper at its
% dpRng/cache defaults must be bit-identical to a direct invz_bz_couplings call given the SAME
% explicit dpRng=30/cache=false -- proving the helper adds nothing beyond defaulting those two
% fields (every other field, here none, forwarded unchanged).
ion = invz_ion();
grid = [6 6 6];

[JnuBz, infoBz, Jaa0Bz] = invz_bz_couplings(ion, struct('grid', grid, 'dpRng', 30, 'cache', false));
[JnuHelp, infoHelp, Jaa0Help] = invz_anchor_couplings(ion, struct('grid', grid, 'cache', false));

verifyEqual(testCase, JnuHelp, JnuBz);
verifyEqual(testCase, infoHelp, infoBz);
verifyEqual(testCase, Jaa0Help, Jaa0Bz);
verifyFalse(testCase, isfield(infoHelp, 'grid'));
end

% =====================================================================
% Presence semantics: dpRng defaults ONLY when the caller's opts omits it
% =====================================================================
function test_dpRng_default_applies_only_when_absent(testCase)
ion = invz_ion();
grid = [4 4 4];

[~, infoDefault] = invz_anchor_couplings(ion, struct('grid', grid, 'cache', false));
verifyEqual(testCase, infoDefault.dpRng, 30);

[~, infoExplicit] = invz_anchor_couplings(ion, struct('grid', grid, 'cache', false, 'dpRng', 12));
verifyEqual(testCase, infoExplicit.dpRng, 12);
end

% =====================================================================
% Ewald injection: complete backend + grid provenance, mechanical only
% =====================================================================
function test_ewald_halfopen_pcomplete_injection_complete_provenance(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
grid = [4 4 4];
opts = struct('grid', grid, 'cache', false, 'dipole', 'ewald', 'ewald', eo, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0], 'gammaPolicy', 'P_complete');

[Jnu, info, Jaa0] = invz_anchor_couplings(ion, opts);

% -- complete backend provenance --
verifyEqual(testCase, sort(fieldnames(info.dipole)), sort({'backend'; 'ewald'; 'q_reduction'; 'primitive_schema'}));
verifyEqual(testCase, info.dipole.backend, 'ewald');
verifyEqual(testCase, info.dipole.ewald, eo);
verifyTrue(testCase, all(isfinite(Jnu(:))));
verifyTrue(testCase, isfinite(Jaa0));

% -- complete grid provenance --
verifyTrue(testCase, isfield(info, 'grid'));
gr = info.grid;
expectedFields = sort({'schema', 'convention', 'offset', 'gammaPolicy', 'requested', 'nominal', ...
                       'retained', 'n_gamma', 'qhash'});
verifyEqual(testCase, sort(fieldnames(gr)), expectedFields(:));
verifyEqual(testCase, gr.convention, 'halfopen');
verifyEqual(testCase, gr.offset, logical([0 0 0]));
verifyEqual(testCase, gr.gammaPolicy, 'P_complete');
verifyEqual(testCase, gr.requested, grid);
verifyEqual(testCase, gr.nominal, prod(grid));
verifyEqual(testCase, gr.retained, prod(grid));   % P_complete: nothing dropped
verifyGreaterThanOrEqual(testCase, gr.n_gamma, 1, ...
    'the [0 0 0]-offset halfopen grid must contain the exact Gamma row.');

% -- dpRng is still helper-defaulted/threaded through to invz_bz_couplings/invz_jq_modes, which
% then (correctly, per invz_jq_modes.m's own documented Ewald contract) reports the canonical
% NaN sentinel: "info.dpRng is NaN (dpRng does not affect the Ewald calculation or its cache
% identity)" -- the helper does not, and must not, override that backend-owned sentinel.
verifyTrue(testCase, isscalar(info.dpRng) && isnan(info.dpRng));
end

function test_ewald_halfopen_pdrop_retained_rowcount_differs(testCase)
% Second, mechanical call: P_drop must retain STRICTLY FEWER rows than P_complete (exact-Gamma
% row(s) removed), same nominal cardinality/backend. No physics computed -- mechanical row-count
% provenance only.
ion = invz_ion();
eo = default_eopts(ion.a);
grid = [4 4 4];
common = struct('grid', grid, 'cache', false, 'dipole', 'ewald', 'ewald', eo, ...
    'gridConvention', 'halfopen', 'gridOffset', [0 0 0]);
optsComplete = common; optsComplete.gammaPolicy = 'P_complete';
optsDrop     = common; optsDrop.gammaPolicy     = 'P_drop';

[~, infoComplete] = invz_anchor_couplings(ion, optsComplete);
[~, infoDrop]      = invz_anchor_couplings(ion, optsDrop);

verifyEqual(testCase, infoDrop.dipole.backend, 'ewald');
verifyEqual(testCase, infoComplete.grid.nominal, infoDrop.grid.nominal);
verifyEqual(testCase, infoComplete.grid.n_gamma, infoDrop.grid.n_gamma);
verifyEqual(testCase, infoComplete.grid.retained, infoComplete.grid.nominal);
verifyEqual(testCase, infoDrop.grid.retained, infoDrop.grid.nominal - infoDrop.grid.n_gamma);
verifyLessThan(testCase, infoDrop.grid.retained, infoComplete.grid.retained);
end

% =====================================================================
% Dipole/Ewald alone (no grid-policy field) must not activate the new grid
% route THROUGH the wrapper either -- the helper adds no routing logic of
% its own beyond invz_bz_couplings' own Task-4 presence rule.
% =====================================================================
function test_ewald_alone_through_helper_does_not_activate_grid_route(testCase)
ion = invz_ion();
eo = default_eopts(ion.a);
[Jnu, info] = invz_anchor_couplings(ion, struct('grid', [4 4 4], 'cache', false, ...
    'dipole', 'ewald', 'ewald', eo));
verifyFalse(testCase, isfield(info, 'grid'), ...
    'opts.dipole=''ewald'' alone (no grid-policy field), forwarded through the helper, must NOT activate the new grid route.');
verifyEqual(testCase, info.dipole.backend, 'ewald');
verifyTrue(testCase, all(isfinite(Jnu(:))));
end
