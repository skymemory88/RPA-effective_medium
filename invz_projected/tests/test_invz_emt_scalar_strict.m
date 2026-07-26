function tests = test_invz_emt_scalar_strict
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function [G0, Sigma, Jnu] = fixture()
nw = 8;
G0    = -(0.5 ./ (1:nw)').^0 * 0.5;      % [nw x 1], O(1) negative (G = -chi)
G0(1) = -0.9;
Sigma = 0.05*ones(nw, 1);
Jnu   = linspace(-2e-3, 6e-3, 24).';
end

% Absent field => legacy path, BIT-IDENTICAL. This is the G9 guard at leaf level.
function test_absent_scheme_is_bit_identical(testCase)
[G0, Sigma, Jnu] = fixture();
a = invz_emt_scalar(G0, Sigma, Jnu, struct());
b = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'resummed'));
verifyEqual(testCase, a.K, b.K, 'AbsTol', 0);
verifyEqual(testCase, a.G, b.G, 'AbsTol', 0);
verifyEqual(testCase, a.medium_status, 'not_applicable');
end

% Strict mode changes slot 1 ONLY.
function test_strict_touches_only_slot_one(testCase)
[G0, Sigma, Jnu] = fixture();
leg = invz_emt_scalar(G0, Sigma, Jnu, struct());
st  = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, st.K(2:end), leg.K(2:end), 'AbsTol', 0);
verifyEqual(testCase, st.G(2:end), leg.G(2:end), 'AbsTol', 0);
verifyTrue(testCase, st.dynamic_converged);
verifyNotEqual(testCase, st.K(1), leg.K(1));
verifyEqual(testCase, st.medium_status, 'ok');
verifyEqual(testCase, st.medium.scheme, 'strict_1z_dyson_ref');
end

% The strict slot is exactly the primitive composition -- no re-derivation inside the leaf.
function test_strict_slot_equals_the_primitives(testCase)
[G0, Sigma, Jnu] = fixture();
st = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
mom  = invz_coupling_moments(Jnu);
Gref = invz_static_medium_reference(G0(1), Sigma(1), 'strict_1z_dyson_ref');
K0   = invz_medium_moment_closure(Gref, mom, 'strict_1z_dyson_ref');
verifyEqual(testCase, st.K(1), K0, 'AbsTol', 0);
verifyEqual(testCase, st.G(1), G0(1)/(1 + Sigma(1) + K0*G0(1)), 'AbsTol', 0);
end

% Supplying Jmom must give the identical answer to deriving it (Task 3 hot-path optimisation
% must not be a semantic change).
function test_supplied_jmom_matches_derived(testCase)
[G0, Sigma, Jnu] = fixture();
o = struct('static_medium', 'strict_1z_dyson_ref');
a = invz_emt_scalar(G0, Sigma, Jnu, o);
o.Jmom = invz_coupling_moments(Jnu);
b = invz_emt_scalar(G0, Sigma, Jnu, o);
verifyEqual(testCase, a.K, b.K, 'AbsTol', 0);
verifyEqual(testCase, a.G, b.G, 'AbsTol', 0);
end

% An out-of-domain reference is a STATUS, not a throw, and it is not reported as
% non-convergence-without-explanation (spec §5.2).
function test_out_of_domain_reference_is_a_status(testCase)
[G0, Sigma, Jnu] = fixture();
Sigma(1) = -1;                                   % 1 + Sigma0 = 0 exactly
med = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
verifyEqual(testCase, med.medium_status, 'ref_denom_nonpositive');
verifyTrue(testCase, isnan(med.K(1)));
verifyFalse(testCase, med.converged);
verifyEqual(testCase, med.medium.ref.status, 'ref_denom_nonpositive');
end

% The bare-reference comparator is selectable and differs at finite Sigma0.
% The two conventions differ only through Gref, so their K(1) separation has a CLOSED FORM:
%   K_dyson - K_bare = -mu2*(Gref_dyson - Gref_bare)
%                    = -mu2*(G0/(1+Sigma0) - G0) = mu2*G0*Sigma0/(1+Sigma0)
% Asserting that identity, rather than "the difference exceeds some number", pins the actual
% O(1/z^2) scheme-choice physics and cannot silently drift if the fixture changes.
%
% The brief originally asserted a >1e-3 relative separation here. That threshold is unreachable
% with this fixture: measured separation is 1.239e-4 (dK = -2.48447205e-07), and clearing 1e-3
% would need Sigma0 >= 1 -- outside the regime, since Sigma0 is O(1/z). The fixture is right and
% the threshold was wrong; replaced on the user's ruling of 2026-07-26.
function test_bare_ref_comparator_differs(testCase)
[G0, Sigma, Jnu] = fixture();
d = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_dyson_ref'));
b = invz_emt_scalar(G0, Sigma, Jnu, struct('static_medium', 'strict_1z_bare_ref'));
mom = invz_coupling_moments(Jnu);
dK   = d.K(1) - b.K(1);
pred = mom.mu2*G0(1)*Sigma(1)/(1 + Sigma(1));
% Tolerance is ABSOLUTE and scaled to the OPERANDS, not relative to the result. dK is the
% difference of two K(1) values of magnitude ~2.005e-3 that yields ~2.48e-7, cancelling about
% four decimal digits, so the error floor is eps(|K(1)|) ~ 4.3e-19 in absolute terms -- i.e. a
% relative floor on dK of ~1.8e-12. A RelTol tighter than that measures cancellation noise
% rather than the identity (a RelTol of 1e-12 fails here at a measured 1.18e-12, despite the
% identity holding). 8*eps(|K|) leaves roughly an order of magnitude of headroom.
verifyEqual(testCase, dK, pred, 'AbsTol', 8*eps(abs(b.K(1))), ...
    'the dyson/bare K(1) separation must equal mu2*G0*Sigma0/(1+Sigma0)');
% Separation floor: independent of the identity above, so a collapse of the two conventions to
% the same value still fails even if `pred` were itself mis-derived. Float noise is ~1e-16, so
% the measured 1.239e-4 clears this by ~12 orders of magnitude.
verifyGreaterThan(testCase, abs(dK)/abs(b.K(1)), 1e-5);
end

% Strict PM + [nJ,nw] remains supported: slot 1 uses column 1 and slots 2:end stay legacy.
function test_strict_matrix_uses_static_column_only(testCase)
[G0, Sigma, Jnu] = fixture();
Jm = repmat(Jnu, 1, numel(G0));
Jm(:,2:end) = Jm(:,2:end) .* (1 + (1:numel(G0)-1));
leg = invz_emt_scalar(G0, Sigma, Jm, struct());
st  = invz_emt_scalar(G0, Sigma, Jm, struct('static_medium', 'strict_1z_dyson_ref'));
mom = invz_coupling_moments(Jm);
mom0 = structfun(@(x) x(1), mom, 'UniformOutput', false);
Gref = invz_static_medium_reference(G0(1), Sigma(1), 'strict_1z_dyson_ref');
K0 = invz_medium_moment_closure(Gref, mom0, 'strict_1z_dyson_ref');
verifyEqual(testCase, st.K(1), K0, 'AbsTol', 0);
verifyEqual(testCase, st.K(2:end), leg.K(2:end), 'AbsTol', 0);
verifyEqual(testCase, st.G(2:end), leg.G(2:end), 'AbsTol', 0);
end
