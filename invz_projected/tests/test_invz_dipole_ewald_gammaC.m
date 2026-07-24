function tests = test_invz_dipole_ewald_gammaC
% Gate-C (Step-4 Task 9): primitive-computable Gate-C checks C1-C5 for the
% opt-in Ewald dipolar primitive invz_dipole_ewald.m, per
% docs/invzp_ewald_prereg.md sec 5 (FROZEN) and the plan's "Task 9" section
% (docs/superpowers/plans/2026-07-24-ewald-dipolar-primitive.md).
%
%   C1/C2 -- boundary convention (per-q filter omits ONLY G=0 at exact
%            Gamma) + regular-Gamma tensor reconstruction (no fourth
%            directional/surface term).
%   C3    -- isolated reciprocal G=0 summand equals the closed-form
%            projector P0_ab(q) at M_id, and participates for nonzero q.
%   C4    -- even/odd analytic remainder about Gamma (frozen contraction
%            bounds; the odd part is REPORTED, never required to vanish).
%   C5    -- the macroscopic G=0 term couples ONLY to the uniform mode
%            v=ones(4,1)/2.
%
% Gate-C checks 6/7 (caller integration/cache/provenance) and the complete
% Cartesian q-path reconstruction in check 4 need invz_jq_modes wiring and
% remain Step 5 -- OUT OF SCOPE here.
%
% This file inlines its OWN local independent parts reconstruction
% (reconstruct_parts below) from PUBLIC geom data. test_invz_dipole_ewald.m
% has a similar-looking helper of the same name, but MATLAB local functions
% are file-scoped, so it is unreachable from here; this is an independent,
% from-scratch reimplementation (matched against the primitive at M_id by
% the very first test below, BEFORE any check dissects the parts) that also
% exposes the ISOLATED G=0 summand C3/C4/C5 need.
%
% All structural/exact-identity comparisons use invz_ewald_metrics' M_id
% (max component |A-B| <= 1e-12*T_scale, T_scale = max component of
% max(|A|,|B|)), called on ONE complete [3,3,ntau,ntau] tensor at a time
% (never a stacked multi-q array) per docs/invzp_ewald_prereg.md sec 3.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));            % invz_projected
addpath(fullfile(here, '..', '..'));      % repo root: invz_dipole_ewald, MF_dipole
addpath(fullfile(here, '..', '..', 'invz_common'));  % invz_ion, invz_const
end

% =====================================================================
% shared helpers (names deliberately free of "test" so functiontests skips
% them -- test_invz_dipole_ewald.m precedent)
% =====================================================================
function eo = mk_eopts(alpha, r_cut, g_cut, boundary)
eo = struct('alpha', alpha, 'r_cut', r_cut, 'g_cut', g_cut, 'boundary', boundary);
end

function eo = eopts_alpha_matched(alpha)
% Frozen Task-9 rays/magnitudes config (brief): alpha=alpha0,
% r_cut=6.5/alpha0, g_cut=13*alpha0 -- the same alpha-matched generous-cutoff
% rule as Gate-B (Task 8); reimplemented independently here per-file.
eo = mk_eopts(alpha, 6.5/alpha, 13*alpha, 'conducting_k0_omitted');
end

function [rays, raylabels] = frozen_rays()
% Five frozen reciprocal-space rays in REDUCED (Miller/hkl) coordinates:
% a*, b*, c*, [1 1 1], [2 1 -1] (prereg sec 5 item 3 / plan Task 9). Since
% invz_dipole_ewald's q is already expressed in units of the reciprocal
% lattice vectors, "ray a*" IS q_hkl=[1 0 0], etc. -- no separate unit
% conversion is needed before scaling by a signed magnitude s.
rays = [1 0 0; 0 1 0; 0 0 1; 1 1 1; 2 1 -1];
raylabels = {'a*', 'b*', 'c*', '[1 1 1]', '[2 1 -1]'};
end

function smags = frozen_smags()
% Frozen signed magnitudes, r.l.u. (prereg sec 5 item 3 / plan Task 9).
% Order matters: ray_mag_grid below relies on this exact
% [+1e-3, -1e-3, +1e-4, -1e-4] ordering to pair +s/-s per ray in C4.
smags = [1e-3, -1e-3, 1e-4, -1e-4];
end

function [Q, rayOf, sOf] = ray_mag_grid(rays, smags)
% Q = [nrays*nsmags, 3], row order: ray outer (r=1..nrays), smags inner
% (si=1..nsmags) -- i.e. rows (r-1)*nsmags+1 .. r*nsmags are ray r's block,
% in the exact smags order. rayOf/sOf are parallel [n,1] index/value arrays
% so callers never have to re-derive (ray,s) from a linear row index by
% hand.
nr = size(rays, 1); ns = numel(smags);
Q = zeros(nr*ns, 3); rayOf = zeros(nr*ns, 1); sOf = zeros(nr*ns, 1);
row = 0;
for r = 1:nr
    for si = 1:ns
        row = row + 1;
        Q(row, :) = smags(si)*rays(r, :);
        rayOf(row) = r;
        sOf(row) = smags(si);
    end
end
end

function [dR, dG, dS, dG0] = reconstruct_parts(geom, qrow, alpha, g_cut)
% Independent (test-only, own copy per-file) reconstruction from PUBLIC geom
% data of the real, full-reciprocal-total, self, and ISOLATED reciprocal
% G=0 contributions. Every ordered (n,m) block is formed on its own; no
% conjugate block is ever copied. dR+dG+dS reproduces the primitive's dip
% exactly (guarded by the very first test in this file); dG0 is the single
% G=[0,0,0] candidate's contribution alone (a SUBSET of dG at nonzero q;
% exactly zero at exact Gamma, where G=0 is excluded from dG by the same
% k~=0 filter the primitive itself applies).
ntau = geom.ntau; B = geom.B; Vc = geom.Vc; taucart = geom.taucart; tau = geom.tau;
Gcart = geom.Gcart; Ghkl = geom.Ghkl;
K = floor(qrow + 0.5); qbar = qrow - K; qcart = qbar*B;
k = Gcart + qcart; kk = sum(k.^2, 2); keep = (kk <= g_cut^2) & (kk > 0);
isG0 = all(Ghkl == 0, 2); g0k = isG0(keep);
selfval = 4*alpha^3/(3*sqrt(pi));
z = complex(zeros(3, 3, ntau, ntau)); dR = z; dG = z; dS = z; dG0 = z;
ksel = k(keep, :); kk2 = kk(keep); kernel = (4*pi/Vc)*exp(-kk2/(4*alpha^2))./kk2;
Gk = Gcart(keep, :);
for n = 1:ntau
    for m = 1:ntau
        d = taucart(m, :) - taucart(n, :);
        gph = exp(-1i*2*pi*(K*(tau(m, :) - tau(n, :)).'));
        X = geom.real{n, m}.x; gab = geom.real{n, m}.gab; ph = exp(-1i*(X*qcart.'));
        bR = zeros(3, 3);
        for aa = 1:3, for bb = 1:3, bR(aa, bb) = -sum(ph.*gab(:, aa, bb)); end, end %#ok<ALIGN>
        dR(:, :, n, m) = gph*bR;
        phG = exp(1i*(Gk*d.')); bG = zeros(3, 3);
        for aa = 1:3, for bb = 1:3, bG(aa, bb) = sum(kernel.*ksel(:, aa).*ksel(:, bb).*phG); end, end %#ok<ALIGN>
        dG(:, :, n, m) = gph*bG;
        if any(g0k)
            ker0 = kernel(g0k); k0 = ksel(g0k, :); ph0 = exp(1i*(Gk(g0k, :)*d.')); b0 = zeros(3, 3);
            for aa = 1:3, for bb = 1:3, b0(aa, bb) = sum(ker0.*k0(:, aa).*k0(:, bb).*ph0); end, end %#ok<ALIGN>
            dG0(:, :, n, m) = gph*b0;
        end
        if n == m, dS(:, :, n, m) = gph*(-selfval*eye(3)); end
    end
end
end

function P0 = isolated_projector(qcart, Vc, alpha)
% Closed-form isolated reciprocal G=0 summand (prereg sec 5 item 3):
%   P0_ab(q) = (4*pi/Vc) qhat_a qhat_b exp(-|q|^2/(4*alpha^2)).
% qcart must be nonzero (qhat undefined at exact Gamma; callers never
% evaluate this at q=0 -- R(0) is handled separately, without a P0 term).
qn = norm(qcart);
assert(qn > 0, 'isolated_projector: qcart must be nonzero (qhat undefined at Gamma).');
qhat = qcart/qn;
P0 = (4*pi/Vc) * (qhat.' * qhat) * exp(-(qn^2)/(4*alpha^2));
end

function P0full = projector_block(qrow, geom, eo, ntau)
% isolated_projector(...) replicated identically over every ordered
% sublattice pair (the G=0 plane wave carries no phase: G=0 dotted with ANY
% displacement is exactly 0, so exp(i*0)=1 for every pair -- verified
% structurally, not merely assumed, by test_gammaC3 below).
K = floor(qrow + 0.5); qbar = qrow - K; qcart = qbar*geom.B;
P0 = isolated_projector(qcart, geom.Vc, eo.alpha);
P0full = complex(zeros(3, 3, ntau, ntau));
for n = 1:ntau
    for m = 1:ntau
        P0full(:, :, n, m) = P0;
    end
end
end

function Vperp = helmert_basis()
% Deterministic machine-orthonormal basis for the orthogonal complement of
% v=ones(4,1)/2 (equivalently of ones(4,1)): the classical Helmert
% contrasts, written out explicitly in closed form (no null/svd sign
% ambiguity) so every column is exactly orthogonal to ones(4,1) and to
% every other column. Columns, before normalization:
%   [1 -1 0 0]', [1 1 -2 0]', [1 1 1 -3]'.
Vperp = [ 1  1  1
         -1  1  1
          0 -2  1
          0  0 -3];
Vperp = Vperp ./ vecnorm(Vperp, 2, 1);
end

% =====================================================================
% Guard (plan Task 9 / brief point 1): FIRST prove the local reconstruction
% reproduces the primitive exactly, BEFORE any check below dissects the
% parts. Every subsequent Gate-C test in this file depends on
% reconstruct_parts; this is its trust anchor.
% =====================================================================
function test_gammaC_reconstruction_matches_primitive_before_dissection(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
eo = eopts_alpha_matched(fx.alpha0);
[rays, ~] = frozen_rays();
smags = frozen_smags();
[Qrays, ~, ~] = ray_mag_grid(rays, smags);
Q = [0 0 0; Qrays];                       % exact Gamma + all 20 ray/magnitude points

[dip, ~, geom] = invz_dipole_ewald(Q, a, tau, eo);
worst_margin = -inf;
for i = 1:size(Q, 1)
    [dR, dGtot, dS, ~] = reconstruct_parts(geom, Q(i, :), eo.alpha, eo.g_cut);
    recon = dR + dGtot + dS;
    mres = M.mid(recon, dip(:, :, :, :, i));
    verifyTrue(testCase, mres.pass, sprintf(...
        ['Gate-C reconstruction guard: real+reciprocal-total+self != primitive at Q row %d ' ...
         '(M_id worst_margin=%.3e). Every later Gate-C check in this file depends on this ' ...
         'reconstruction matching first.'], i, mres.worst_margin));
    worst_margin = max(worst_margin, mres.worst_margin);
end
fprintf(['Gate-C reconstruction guard: worst M_id margin over exact Gamma + 20 frozen ray/' ...
    'magnitude points = %.3e.\n'], worst_margin);
end

% =====================================================================
% C1: boundary convention -- the per-q reciprocal filter omits ONLY G=0 at
% exact Gamma (prereg sec 5 item 2 "half"; plan Task 9 C1/C2).
% =====================================================================
function test_gammaC1_filter_omits_only_G0_at_exact_gamma(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau;
fx = invz_ewald_fixtures();
eo = eopts_alpha_matched(fx.alpha0);
[~, ~, geom] = invz_dipole_ewald([0 0 0], a, tau, eo);

isG0 = all(geom.Ghkl == 0, 2);
verifyEqual(testCase, nnz(isG0), 1, ...
    'Gate-C1: the candidate union must contain exactly one G=[0,0,0] row.');

% qcart==0 at exact Gamma, so k=Gcart+0=Gcart identically.
kk = sum(geom.Gcart.^2, 2);
keep = (kk <= eo.g_cut^2) & (kk > 0);
verifyFalse(testCase, keep(isG0), 'Gate-C1: G=0 candidate must be OMITTED at exact Gamma.');

gmag = vecnorm(geom.Gcart, 2, 2);
inside = gmag <= eo.g_cut;
verifyEqual(testCase, nnz(~keep & inside), 1, sprintf(...
    ['Gate-C1: exactly one in-sphere (|G|<=g_cut) candidate may be omitted at exact Gamma ' ...
     '(G=0); found %d.'], nnz(~keep & inside)));
verifyTrue(testCase, all(keep(inside & ~isG0)), ...
    'Gate-C1: every non-G=0 in-sphere candidate must be RETAINED at exact Gamma.');
end

% =====================================================================
% C2: regular-Gamma tensor -- dip_reg(0) reconstruction equals the returned
% primitive, and there is no directional/surface term (prereg sec 5 items
% 1-2; plan Task 9 C1/C2).
% =====================================================================
function test_gammaC2_regular_gamma_reconstruction_no_fourth_term(testCase)
% Frozen Gate-C check 1 (k=0 boundary term vs the Gate-B oracle at M_FD) is
% CITED/CONSUMED, not rerun: it is independently retained by
% test_gateB_scalar_oracle_agreement in test_invz_dipole_ewald_ref.m
% (Step-4 Task 8, docs/invzp_ewald_prereg.md sec 4 + sec 5 item 1). Task 9
% is scoped to the primitive-computable structural/M_id checks below, which
% this file does not duplicate against invz_scalar_ewald_ref.
%
% This test reconstructs
%   dip_reg(0) = dip^(r)(0) + sum_{G~=0} dip^(G)(0) + dip^(self)
% from PUBLIC geom data (reconstruct_parts above) and requires it to equal
% the returned exact-Gamma primitive at M_id. Because the reconstruction is
% built ONLY from the real-space screened sum, the reciprocal sum over
% G~=0 (G=0 is excluded by construction -- proved separately by
% test_gammaC1 above), and the closed-form self term -- none of which take
% a "sample shape" or "approach direction" argument -- an M_id pass PROVES
% the returned exact-Gamma tensor carries no fourth (directional, surface,
% shape/demag-dependent Lorentz) term: any such term would show up as an
% unaccounted-for residual.
ion = invz_ion(); a = ion.a; tau = ion.tau;
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
eo = eopts_alpha_matched(fx.alpha0);

[dip0, ~, geom] = invz_dipole_ewald([0 0 0], a, tau, eo);
[dR, dGtot, dS, dG0] = reconstruct_parts(geom, [0 0 0], eo.alpha, eo.g_cut);
verifyEqual(testCase, max(abs(dG0(:))), 0, ...
    'Gate-C2: the isolated G=0 term must be exactly absent (zero) at exact Gamma.');

recon = dR + dGtot + dS;
mres = M.mid(recon, dip0);
verifyTrue(testCase, mres.pass, sprintf(...
    ['Gate-C2: dip_reg(0) reconstruction (real+recip[G~=0]+self) != returned exact-Gamma ' ...
     'primitive at M_id (worst_margin=%.3e) -- implies an undocumented fourth term.'], ...
    mres.worst_margin));

% Structural corollary: dip_reg(0) is real (a directional/shape-dependent
% Lorentz-type addition would generically not preserve this).
Tscale = max(abs(dip0(:)));
verifyLessThanOrEqual(testCase, max(abs(imag(dip0(:)))), 1e-12*Tscale, ...
    'Gate-C2: exact-Gamma primitive must be real at M_id.');

% "No directional term": the assembly loop processes every q row
% independently (no cross-row state), so the exact-Gamma slice of a batch
% that ALSO contains the 20 frozen nonzero ray/magnitude points must be
% IDENTICAL to the standalone Gamma-only call -- the returned Gamma tensor
% cannot depend on which other q's/approach directions happen to share its
% call.
[rays, ~] = frozen_rays();
[Qrays, ~, ~] = ray_mag_grid(rays, frozen_smags());
dipBatch = invz_dipole_ewald([0 0 0; Qrays], a, tau, eo);
mdir = M.mid(dip0, dipBatch(:, :, :, :, 1));
verifyTrue(testCase, mdir.pass, sprintf(...
    ['Gate-C2: exact-Gamma primitive changed when computed alongside other q rows ' ...
     '(M_id worst_margin=%.3e) -- implies a directional/context-dependent term.'], ...
    mdir.worst_margin));

fprintf(['Gate-C2 regular-Gamma reconstruction: M_id worst_margin=%.3e; ' ...
    'context-independence M_id worst_margin=%.3e.\n'], mres.worst_margin, mdir.worst_margin);
end

% =====================================================================
% C3: isolated projector -- for every signed ray/magnitude, every ordered
% pair, every Cartesian component, the extracted isolated G=0 summand
% equals the closed-form projector at M_id, and the G=0 candidate
% participates at nonzero q (prereg sec 5 item 3; plan Task 9 C3).
% Subtracting the projector and checking only isfinite is NOT this gate.
% =====================================================================
function test_gammaC3_isolated_G0_projector_and_participation(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau, 1);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
eo = eopts_alpha_matched(fx.alpha0);
[rays, raylabels] = frozen_rays();
smags = frozen_smags();
[Q, rayOf, sOf] = ray_mag_grid(rays, smags);

[~, ~, geom] = invz_dipole_ewald(Q, a, tau, eo);
isG0 = all(geom.Ghkl == 0, 2);

worst_margin = -inf;
for i = 1:size(Q, 1)
    qrow = Q(i, :);
    lbl = sprintf('ray %s, s=%.4g', raylabels{rayOf(i)}, sOf(i));

    % ---- participation: G=0 candidate is retained by the per-q filter ----
    K = floor(qrow + 0.5); qbar = qrow - K; qcart = qbar*geom.B;
    kk = sum((geom.Gcart + qcart).^2, 2);
    keep = (kk <= eo.g_cut^2) & (kk > 0);
    verifyTrue(testCase, keep(isG0), sprintf(...
        'Gate-C3: G=0 candidate must PARTICIPATE at nonzero q (%s).', lbl));

    % ---- isolated G=0 summand, extracted from the local reconstruction ---
    [~, ~, ~, dG0] = reconstruct_parts(geom, qrow, eo.alpha, eo.g_cut);
    verifyGreaterThan(testCase, max(abs(dG0(:))), 0, sprintf(...
        'Gate-C3: isolated G=0 contribution must be nonzero at nonzero q (%s).', lbl));

    % ---- compare to the closed-form projector, replicated over every pair
    P0full = projector_block(qrow, geom, eo, ntau);
    mres = M.mid(dG0, P0full);
    verifyTrue(testCase, mres.pass, sprintf(...
        ['Gate-C3: isolated G=0 summand != (4*pi/Vc)*qhat*qhat''*exp(-|q|^2/(4a^2)) at M_id ' ...
         '(%s, worst_margin=%.3e).'], lbl, mres.worst_margin));
    worst_margin = max(worst_margin, mres.worst_margin);
end
fprintf('Gate-C3 isolated-G0 projector: worst M_id margin over all 20 ray/magnitude points = %.3e.\n', ...
    worst_margin);
end

% =====================================================================
% C4: even/odd analytic remainder about Gamma (prereg sec 5 item 4; plan
% Task 9 C4). R_odd is REPORTED, never required to vanish -- raw
% non-Bravais off-diagonal blocks carry an O(q) odd-in-q imaginary term.
% =====================================================================
function test_gammaC4_even_odd_analytic_remainder(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau, 1);
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
eo = eopts_alpha_matched(fx.alpha0);
[rays, raylabels] = frozen_rays();
smags = frozen_smags();                    % [+1e-3 -1e-3 +1e-4 -1e-4], fixed order
nrays = size(rays, 1);

% ---- R(0) = dip_reg(0): a single shared reference tensor, no P0 term ------
R0 = invz_dipole_ewald([0 0 0], a, tau, eo);

% ---- batch every signed ray/magnitude point in one call -------------------
[Q, ~, ~] = ray_mag_grid(rays, smags);
[dip, ~, geom] = invz_dipole_ewald(Q, a, tau, eo);

diff_even_3 = zeros(nrays, 1); diff_even_4 = zeros(nrays, 1);
diff_odd_3  = zeros(nrays, 1); diff_odd_4  = zeros(nrays, 1);
Tscale_even_4 = zeros(nrays, 1);

for r = 1:nrays
    base = (r-1)*4;
    Rp3 = dip(:, :, :, :, base+1) - projector_block(Q(base+1, :), geom, eo, ntau);
    Rm3 = dip(:, :, :, :, base+2) - projector_block(Q(base+2, :), geom, eo, ntau);
    Rp4 = dip(:, :, :, :, base+3) - projector_block(Q(base+3, :), geom, eo, ntau);
    Rm4 = dip(:, :, :, :, base+4) - projector_block(Q(base+4, :), geom, eo, ntau);

    Reven3 = (Rp3 + Rm3)/2; Rodd3 = (Rp3 - Rm3)/2;
    Reven4 = (Rp4 + Rm4)/2; Rodd4 = (Rp4 - Rm4)/2;

    diff_even_3(r) = max(abs(Reven3(:) - R0(:)));
    diff_even_4(r) = max(abs(Reven4(:) - R0(:)));
    diff_odd_3(r)  = max(abs(Rodd3(:)));
    diff_odd_4(r)  = max(abs(Rodd4(:)));

    me4 = M.mid(Reven4, R0);              % harvest T_scale(0;R_even(s),R(0)) at s=1e-4
    Tscale_even_4(r) = me4.T_scale;

    fprintf(['Gate-C4 (ray %s): E_even(1e-3)=%.3e E_even(1e-4)=%.3e | R_odd(1e-3) max=%.3e ' ...
        'R_odd(1e-4) max=%.3e (REPORTED, not required to vanish).\n'], raylabels{r}, ...
        diff_even_3(r), diff_even_4(r), diff_odd_3(r), diff_odd_4(r));
end

E_even_3 = max(diff_even_3); E_even_4 = max(diff_even_4);
E_odd_3  = max(diff_odd_3);  E_odd_4  = max(diff_odd_4);
Tscale_even4_agg = max(Tscale_even_4);      % max_rays T_scale(0;R_even(1e-4),R(0))
A_T4 = 1e-8*Tscale_even4_agg;

verifyLessThanOrEqual(testCase, E_even_4, 1e-6*Tscale_even4_agg, sprintf(...
    'Gate-C4: E_even(1e-4)=%.6e exceeds 1e-6*Tscale_even(1e-4)=%.6e.', ...
    E_even_4, 1e-6*Tscale_even4_agg));
verifyLessThanOrEqual(testCase, E_even_4, 0.02*E_even_3 + A_T4, sprintf(...
    'Gate-C4: E_even(1e-4)=%.6e exceeds 0.02*E_even(1e-3)+A_T(1e-4)=%.6e.', ...
    E_even_4, 0.02*E_even_3 + A_T4));
verifyLessThanOrEqual(testCase, E_odd_4, 0.20*E_odd_3 + A_T4, sprintf(...
    'Gate-C4: E_odd(1e-4)=%.6e exceeds 0.20*E_odd(1e-3)+A_T(1e-4)=%.6e.', ...
    E_odd_4, 0.20*E_odd_3 + A_T4));

fprintf(['Gate-C4 even/odd remainder summary: E_even(1e-3)=%.3e E_even(1e-4)=%.3e ' ...
    'E_odd(1e-3)=%.3e E_odd(1e-4)=%.3e (R_odd REPORTED, not required to vanish) ' ...
    'A_T(1e-4)=%.3e.\n'], E_even_3, E_even_4, E_odd_3, E_odd_4, A_T4);
end

% =====================================================================
% C5: uniform-mode support -- the macroscopic G=0 term couples ONLY to
% v=ones(4,1)/2 (prereg sec 5 item 5; plan Task 9 C5). Uses the ISOLATED
% implementation G=0 contribution Delta (reconstruct_parts' dG0 output),
% NOT dip(q)-dip(0) -- that difference also contains the finite-q analytic
% remainder C4 tests separately.
% =====================================================================
function test_gammaC5_uniform_mode_support(testCase)
ion = invz_ion(); a = ion.a; tau = ion.tau; ntau = size(tau, 1);
verifyEqual(testCase, ntau, 4, 'Gate-C5: v=ones(4,1)/2 assumes the frozen ntau=4 LiHoF4 basis.');
fx = invz_ewald_fixtures();
M  = invz_ewald_metrics();
eo = eopts_alpha_matched(fx.alpha0);
[rays, raylabels] = frozen_rays();
smags = frozen_smags();
[Q, rayOf, sOf] = ray_mag_grid(rays, smags);

v = ones(4, 1)/2;
Vperp = helmert_basis();
verifyLessThanOrEqual(testCase, max(abs(Vperp.'*Vperp - eye(3)), [], 'all'), 1e-13, ...
    'Gate-C5: Vperp must be machine-orthonormal.');
verifyLessThanOrEqual(testCase, max(abs(v.'*Vperp)), 1e-13, ...
    'Gate-C5: Vperp must be orthogonal to v (equivalently to ones(4,1)).');

[~, ~, geom] = invz_dipole_ewald(Q, a, tau, eo);
fourpiVc = 4*pi/geom.Vc;

% The prereg's compact "v'*Delta*v = 4*(4*pi/Vc)*qhat_c^2" retains the
% Gaussian screening factor implicitly: Delta_cc IS, algebraically, exactly
% the closed-form projector P0_cc(q) replicated over every pair (proved
% structurally by C3 above), so at the exact-identity (M_id) scale the
% target below MUST retain the same exp(-|q|^2/(4*alpha^2)) factor used to
% build Delta itself -- dropping it fails M_id by ~|q|^2/(4*alpha^2), which
% at these frozen magnitudes (e.g. ray [2 1 -1], s=1e-4) is ~2.7e-7, four
% orders of magnitude above the 1e-12 relative floor.
worst_vDv_margin = -inf; worst_perp_ratio = -inf;
for i = 1:size(Q, 1)
    qrow = Q(i, :);
    lbl = sprintf('ray %s, s=%.4g', raylabels{rayOf(i)}, sOf(i));

    [~, ~, ~, dG0] = reconstruct_parts(geom, qrow, eo.alpha, eo.g_cut);
    Deltacc = squeeze(dG0(3, 3, :, :));       % [ntau,ntau]; c = Cartesian index 3

    K = floor(qrow + 0.5); qbar = qrow - K; qcart = qbar*geom.B;
    qhat = qcart/norm(qcart);
    target = 4*fourpiVc*qhat(3)^2*exp(-(norm(qcart)^2)/(4*eo.alpha^2));

    vDv = v.'*Deltacc*v;
    mres = M.mid(vDv, target);
    verifyTrue(testCase, mres.pass, sprintf(...
        ['Gate-C5: v''*Delta_cc*v != 4*(4*pi/Vc)*qhat_c^2*exp(-|q|^2/(4a^2)) at M_id (%s, ' ...
         'worst_margin=%.3e).'], lbl, mres.worst_margin));
    worst_vDv_margin = max(worst_vDv_margin, mres.worst_margin);

    Mperp = Vperp.'*Deltacc*Vperp;
    perp_max = max(abs(Mperp(:)));
    verifyLessThanOrEqual(testCase, perp_max, 1e-4*fourpiVc, sprintf(...
        'Gate-C5: Vperp''*Delta_cc*Vperp exceeds 1e-4*(4*pi/Vc) (%s, max=%.3e, bound=%.3e).', ...
        lbl, perp_max, 1e-4*fourpiVc));
    worst_perp_ratio = max(worst_perp_ratio, perp_max/fourpiVc);
end
fprintf(['Gate-C5 uniform-mode support: worst v''*Delta*v M_id margin=%.3e; worst ' ...
    'Vperp''*Delta*Vperp / (4*pi/Vc) = %.3e (bound 1e-4).\n'], worst_vDv_margin, worst_perp_ratio);
end
