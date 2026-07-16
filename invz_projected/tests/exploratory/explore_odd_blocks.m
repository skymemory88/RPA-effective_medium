% explore_odd_blocks.m  (ODD preflight P0.3) -- ODD off-diagonal block symmetry.
% READ-ONLY exploratory script. Not a test (lives in tests/exploratory/).
%
% Purpose (ODD main-body plan P0.3): using MF_dipole directly (dpRng 30, geom
% reuse), measure the c-a / c-b off-diagonal dipole blocks J^{c,alpha}(q) =
% -C.gfac*dip(3,alpha,:,:) [meV]:
%   (i)   on-axis rays [q 0 0] and [0 0 q], q in {1e-1,1e-2,1e-3}: max|J^{ca}|,
%         max|J^{cb}| -> confirm C2-about-c decay to 0; rate + residual vs Jcc0.
%   (ii)  tilted ray q*[1 0 1]/sqrt(2): the NON-decaying direction-dependent
%         macroscopic limit vs the 4*pi*C.gfac/ion.Vc shape scale.
%   (iii) generic q = [0.31 0.17 0.09]: block magnitudes vs Jcc0.
%   (iv)  smallest shells of the standard 16^3 grid along a* and c*.
% Also pins info.Jcc0 / info.Jaa0 (dpRng 30) as the jcc0/jaa0 anchors.
%
% Run:  matlab -batch "run('.../invz_projected/tests/exploratory/explore_odd_blocks.m')"

here = fileparts(mfilename('fullpath'));   % .../invz_projected/tests/exploratory
addpath(fullfile(here, '..', '..'));       % invz_projected module
addpath(fullfile(here, '..', '..', '..')); % repo root: MF_dipole, exchange

ion = invz_ion();
C   = invz_const();
dpRng = 30;

fprintf('=== P0.3 explore_odd_blocks  (%s) ===  dpRng=%d\n', char(datetime('now')), dpRng);

% Prime the q-independent geometry once, reuse for every q (5-arg MF_dipole).
[~, ~, geomD] = MF_dipole([0 0 0], dpRng, ion.a, ion.tau);

% Reference couplings (info block, dpRng 30) from invz_jq_modes.
[~, info] = invz_jq_modes(ion, [0.25 0 0], struct('dpRng', dpRng, 'cache', true));
Jcc0     = info.Jcc0;         Jaa0     = info.Jaa0;
Jcc0_dip = info.Jcc0_dipole;  Jaa0_dip = info.Jaa0_dipole;
fprintf('\n-- info (dpRng %d) --\n', dpRng);
fprintf('Jcc0=%.17g meV  Jaa0=%.17g meV\n', Jcc0, Jaa0);
fprintf('Jcc0_dipole=%.17g meV (published anchor 6.821e-3)  Jaa0_dipole=%.17g meV (3.912e-3)\n', ...
    Jcc0_dip, Jaa0_dip);
fprintf('ANCHOR jcc0 = %.17g\nANCHOR jaa0 = %.17g\n', Jcc0, Jaa0);

macroscale = 4*pi*C.gfac/ion.Vc;    % direction-dependent macroscopic dipole shape scale
fprintf('4*pi*C.gfac/ion.Vc = %.6g meV = %.4g ueV  (small-q macroscopic scale)\n', ...
    macroscale, macroscale*1e3);

% ---------------------------------------------------------------------------
% (i) on-axis rays
% ---------------------------------------------------------------------------
qs = [1e-1 1e-2 1e-3];
maxca_a = zeros(1,3); maxcb_a = zeros(1,3);
maxca_c = zeros(1,3); maxcb_c = zeros(1,3);
for k = 1:3
    maxca_a(k) = blkmax([qs(k) 0 0], 3, 1, C.gfac, dpRng, ion, geomD);
    maxcb_a(k) = blkmax([qs(k) 0 0], 3, 2, C.gfac, dpRng, ion, geomD);
    maxca_c(k) = blkmax([0 0 qs(k)], 3, 1, C.gfac, dpRng, ion, geomD);
    maxcb_c(k) = blkmax([0 0 qs(k)], 3, 2, C.gfac, dpRng, ion, geomD);
end
fprintf('\n-- (i) on-axis decay  [meV] --\n');
fprintf(' q        maxJca[q00]   maxJcb[q00]   maxJca[00q]   maxJcb[00q]\n');
for k = 1:3
    fprintf('%.0e   %.6e   %.6e   %.6e   %.6e\n', qs(k), maxca_a(k), maxcb_a(k), maxca_c(k), maxcb_c(k));
end
fprintf('decay ratio maxJca[q00]: (1e-1)/(1e-2)=%.3g  (1e-2)/(1e-3)=%.3g  (linear-in-q => ~10)\n', ...
    maxca_a(1)/maxca_a(2), maxca_a(2)/maxca_a(3));
fprintf('residual maxJca[q00] @ q=1e-3 / Jcc0 = %.3e   (Task2 gate <= 1e-6*Jcc0)\n', maxca_a(3)/Jcc0);
fprintf('ANCHOR odd_onaxis_smallq.q     = [%.0e %.0e %.0e]  (along [q 0 0])\n', qs(1), qs(2), qs(3));
fprintf('ANCHOR odd_onaxis_smallq.maxca = [%.17g %.17g %.17g]\n', maxca_a(1), maxca_a(2), maxca_a(3));

% ---------------------------------------------------------------------------
% (ii) tilted ray q*[1 0 1]/sqrt(2): NON-decaying macroscopic limit
% ---------------------------------------------------------------------------
dir = [1 0 1]/sqrt(2);
maxca_t = zeros(1,3);
for k = 1:3, maxca_t(k) = blkmax(qs(k)*dir, 3, 1, C.gfac, dpRng, ion, geomD); end
fprintf('\n-- (ii) tilted ray q*[1 0 1]/sqrt(2)  [meV] --\n');
for k = 1:3
    fprintf(' q=%.0e:  maxJca=%.6e = %.4g ueV   (ratio to 4pi*gfac/Vc = %.4g)\n', ...
        qs(k), maxca_t(k), maxca_t(k)*1e3, maxca_t(k)/macroscale);
end
fprintf('=> tilted-ray limit does NOT decay: bounds small-q assertions (on-axis ONLY).\n');

% ---------------------------------------------------------------------------
% (iii) generic q
% ---------------------------------------------------------------------------
qg = [0.31 0.17 0.09];
maxca_g = blkmax(qg, 3, 1, C.gfac, dpRng, ion, geomD);
maxcb_g = blkmax(qg, 3, 2, C.gfac, dpRng, ion, geomD);
fprintf('\n-- (iii) generic q = [0.31 0.17 0.09]  [meV] --\n');
fprintf('maxJca=%.17g (%.4g ueV, %.3g*Jcc0)   maxJcb=%.6e (%.4g ueV)\n', ...
    maxca_g, maxca_g*1e3, maxca_g/Jcc0, maxcb_g, maxcb_g*1e3);
fprintf('ANCHOR odd_generic_q_max = %.17g\n', maxca_g);

% ---------------------------------------------------------------------------
% (iv) smallest shells of the standard 16^3 grid along a* and c*
% ---------------------------------------------------------------------------
q_a16 = [1/16 0 0];  q_c16 = [0 0 1/16];
fprintf('\n-- (iv) smallest 16^3 shells  [meV] --\n');
fprintf('[1/16 0 0]: maxJca=%.6e (%.3g*Jcc0)   [0 0 1/16]: maxJca=%.6e (%.3g*Jcc0)\n', ...
    blkmax(q_a16,3,1,C.gfac,dpRng,ion,geomD), blkmax(q_a16,3,1,C.gfac,dpRng,ion,geomD)/Jcc0, ...
    blkmax(q_c16,3,1,C.gfac,dpRng,ion,geomD), blkmax(q_c16,3,1,C.gfac,dpRng,ion,geomD)/Jcc0);

fprintf('\n=== explore_odd_blocks done ===\n');

% ------- local helper -------
function m = blkmax(q, mu, nu, gfac, dpRng, ion, geomD)
% max |J^{mu,nu}(q)| element of the 4x4 sublattice block, J = -gfac*dip(mu,nu,:,:).
dip = MF_dipole(q, dpRng, ion.a, ion.tau, geomD);   % [3,3,4,4], Ang^-3
blk = -gfac*squeeze(dip(mu, nu, :, :));             % 4x4, meV
m   = max(abs(blk), [], 'all');
end
