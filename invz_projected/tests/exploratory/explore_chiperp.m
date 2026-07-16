% explore_chiperp.m  (ODD preflight P0.2) -- units + transverse susceptibility chi_perp.
% READ-ONLY exploratory script (no module code touched). Not a test: lives in
% tests/exploratory/, which runtests(...,'tests') never recurses into.
%
% Purpose (ODD main-body plan P0.2):
%   * pin the transverse (a,b) single-ion susceptibility block at (1.53 K, B=0)
%     as the chi_perp anchor for Task 3 (invz_chiperp);
%   * split symmetric/antisymmetric (gyrotropic) parts and elastic on/off shares;
%   * sweep Bx = 0:6 T at 0.31 K (chi_aa, chi_bb) -- smoothness + P0 anchor;
%   * verify dimensional closure: chi0 in meV^-1, J^{ca}/Jcc0 in meV (C.gfac
%     carries mu0/4pi*(gL*muB)^2), so E1 delta-J = V*chi*V needs NO extra g-factors.
%
% Run:  matlab -batch "run('.../invz_projected/tests/exploratory/explore_chiperp.m')"
% Digits printed %.17g so the anchor literals round-trip the exact doubles.

here = fileparts(mfilename('fullpath'));   % .../invz_projected/tests/exploratory
addpath(fullfile(here, '..', '..'));       % invz_projected module
addpath(fullfile(here, '..', '..', '..')); % repo root: MF_dipole, exchange

ion = invz_ion();
C   = invz_const();

fprintf('=== P0.2 explore_chiperp  (%s) ===\n', char(datetime('now')));
fprintf('git: run `git rev-parse --short HEAD` alongside; dpRng not used here (single-ion only)\n');

% ---------------------------------------------------------------------------
% (1) (a,b) block of the full electronuclear chi0z at (1.53 K, B=0, hyp=true).
%     Build si with the DEFAULT single-ion options that invz_chiperp will use
%     (hyp=true, Jxx0=ion.Jxx0, transverse_mf='legacy_x'); at B=0, <Jx>=0 so the
%     transverse MF vanishes and the choice is immaterial.
% ---------------------------------------------------------------------------
T = 1.53;  B = [0 0 0];
si = invz_single_ion(ion, T, invz_field_vec(B), struct('hyp', true));
c0_el  = invz_chi0z(si, T, 0, struct('elastic', true));    % 3x3, meV^-1
c0_nel = invz_chi0z(si, T, 0, struct('elastic', false));

blk_el  = c0_el(1:2, 1:2, 1);      % (a,b)x(a,b) = Cartesian x,y block
blk_nel = c0_nel(1:2, 1:2, 1);

% Symmetrize exactly as invz_chiperp (Task 3) will: Hermitian part then real().
Xp      = real((blk_el  + blk_el' )/2);
Xp_nel  = real((blk_nel + blk_nel')/2);
asymM   = (blk_el - blk_el')/2;                 % anti-Hermitian (gyrotropic) part
maxasym = max(abs(asymM(:)));
maximag = max(abs(imag(blk_el(:))));            % should be ~0 at B=0
elastic_share = (Xp(1,1) - Xp_nel(1,1)) / Xp(1,1);

fprintf('\n-- chi_perp @ (1.53 K, 0 T), hyp, elastic on --\n');
fprintf('Xp(1,1)=%.17g  Xp(2,2)=%.17g  meV^-1\n', Xp(1,1), Xp(2,2));
fprintf('Xp(1,2)=%.17g  Xp(2,1)=%.17g\n', Xp(1,2), Xp(2,1));
fprintf('max|imag block|=%.6g   max|antisym|=%.6g   |Xp(1,1)-Xp(2,2)|=%.3g\n', ...
    maximag, maxasym, abs(Xp(1,1)-Xp(2,2)));
fprintf('elastic_share (chi_aa)=%.17g   (inelastic-only Xp(1,1)=%.17g)\n', ...
    elastic_share, Xp_nel(1,1));
fprintf('EXPECT chi_aa in 16-17 meV^-1 band (11 => truncation, x2 off => convention slip => STOP)\n');

fprintf('\nANCHOR chiperp_1p53K_0T = [%.17g %.17g; %.17g %.17g]\n', ...
    Xp(1,1), Xp(1,2), Xp(2,1), Xp(2,2));
fprintf('ANCHOR chiperp_asym_1p53K = %.17g\n', maxasym);
fprintf('ANCHOR chiperp_elastic_share_1p53K = %.17g\n', elastic_share);

% ---------------------------------------------------------------------------
% (2) Bx sweep at 0.31 K: chi_aa, chi_bb (elastic on, hyp, default MF).
% ---------------------------------------------------------------------------
T2 = 0.31;  Bx = 0:6;
chi_aa = zeros(1, numel(Bx));  chi_bb = zeros(1, numel(Bx));
mf_conv = false(1, numel(Bx)); mf_res = zeros(1, numel(Bx));
for i = 1:numel(Bx)
    sib  = invz_single_ion(ion, T2, invz_field_vec([Bx(i) 0 0]), struct('hyp', true));
    cb   = invz_chi0z(sib, T2, 0, struct('elastic', true));
    blkb = cb(1:2, 1:2, 1);
    symb = real((blkb + blkb')/2);
    chi_aa(i) = symb(1,1);
    chi_bb(i) = symb(2,2);
    if isfield(sib, 'mf_converged'), mf_conv(i) = sib.mf_converged; end
    if isfield(sib, 'mf_residual'),  mf_res(i)  = sib.mf_residual;  end
end
fprintf('\n-- Bx sweep @ 0.31 K (elastic on, hyp) --\n');
for i = 1:numel(Bx)
    fprintf('Bx=%d T:  chi_aa=%.17g   chi_bb=%.17g   mf_conv=%d  mf_res=%.2g\n', ...
        Bx(i), chi_aa(i), chi_bb(i), mf_conv(i), mf_res(i));
end
reldiff = max(abs(diff(chi_aa)) ./ abs(chi_aa(1:end-1)));
fprintf('max relative step |diff(chi_aa)/chi_aa| = %.4g   (Task3 gate < 0.25)\n', reldiff);
fprintf('ANCHOR chiperp_0p31K_Bx.Bx   = [%s]\n', num2str(Bx, '%d '));
fprintf('ANCHOR chiperp_0p31K_Bx.chi_aa = [');  fprintf('%.17g ', chi_aa);  fprintf(']\n');
fprintf('       chiperp_0p31K_Bx.chi_bb = [');  fprintf('%.17g ', chi_bb);  fprintf(']\n');

% ---------------------------------------------------------------------------
% (3) Dimensional closure cross-check.
% ---------------------------------------------------------------------------
fprintf('\n-- dimensional closure --\n');
fprintf('C.gfac = %.6g meV*Ang^3 (carries mu0/4pi*(gL*muB)^2)\n', C.gfac);
fprintf('C.gfac*4/ion.Vc = %.6g meV  (== J_D, matches invz_const comment 1.1654e-3)\n', C.gfac*4/ion.Vc);
fprintf('=> J^{ca} = -gfac*dip [meV], chi0 [meV^-1]; deltaJ=V*chi*V is meV. No extra g-factors.\n');

fprintf('\n=== explore_chiperp done ===\n');
