% explore_tensor_blocks.m  (A0 preflight, Task 2 of docs/superpowers/plans/2026-07-17-invz-tensor-full.md)
% READ-ONLY exploratory script. Not a test (lives in tests/exploratory/, dev-time
% only). Uses invz_odd_blocks (NOT raw MF_dipole) so the tensor anchors are
% measured through the same production entry point invz_tensor will eventually
% mirror block-for-block (Task 3's parity test test_block_parity_with_projected
% _odd_blocks). Explicit addpath below: invz_projected is NOT normally on the
% path for invz_tensor work; this script is the one dev-time exception.
%
% Purpose (plan Task 2 interface): measure, at dpRng 30 (production default),
%   (i)   on-axis rays [q 0 0] (primary) and [0 0 q] (bonus context), q in
%         {1e-1,1e-2,1e-3}: max|J^{ca}| -> LINEAR small-q decay (C2-about-c).
%   (ii)  tilted ray q*[1 0 1]/sqrt(2), same q's: the NON-decaying direction-
%         dependent macroscopic limit vs the 4*pi*C.gfac/ion.Vc shape scale.
%   (iii) generic q = [0.31 0.17 0.09]: block magnitude vs Jcc0.
%   (iv)  Gamma info (Jcc0, Jaa0) at dpRng 30.
%   (v)   smallest shells of the standard 16^3 grid along a*/c*, dpRng in
%         {20, 30, 40}: sensitivity of the ODD residual to the real-space
%         dipole-sum cutoff (NEW vs the superseded P0.3 preflight, which only
%         measured a single dpRng).
% Values printed here (format long g / %.17g) are hand-copied into the frozen
% invzt_anchors.m literals (Step 2) -- this script is never called by any suite.
%
% Run:  matlab -batch "run('.../invz_tensor/tests/exploratory/explore_tensor_blocks.m')"

format long g

here = fileparts(mfilename('fullpath'));                        % invz_tensor/tests/exploratory
addpath(fullfile(here, '..', '..', '..'));                       % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', '..', 'invz_projected'));     % invz_odd_blocks (dev-time only)
addpath(fullfile(here, '..', '..', '..', 'invz_common'));        % invz_ion, invz_const

ion = invz_ion();
C   = invz_const();
dpRng0 = 30;

fprintf('=== A0 explore_tensor_blocks  (%s) ===  git %s\n', char(datetime('now')), git_hash());
fprintf('dpRng (primary) = %d\n', dpRng0);

macroscale = 4*pi*C.gfac/ion.Vc;    % direction-dependent macroscopic dipole shape scale
fprintf('4*pi*C.gfac/ion.Vc = %.17g meV = %.6g ueV  (small-q macroscopic scale)\n', ...
    macroscale, macroscale*1e3);

% ---------------------------------------------------------------------------
% One combined invz_odd_blocks call at dpRng 30 covering (i)/(ii)/(iii)/(v):
% geometry is primed ONCE per call (dominant cost at dpRng 30), so batching all
% q rows into one call is materially faster than one call per q.
% ---------------------------------------------------------------------------
qs = [1e-1 1e-2 1e-3];
dirT = [1 0 1]/sqrt(2);
q_all = [ qs(1) 0 0; qs(2) 0 0; qs(3) 0 0; ...            % rows 1-3:  on-axis a*
          0 0 qs(1); 0 0 qs(2); 0 0 qs(3); ...            % rows 4-6:  on-axis c* (bonus)
          qs(1)*dirT; qs(2)*dirT; qs(3)*dirT; ...         % rows 7-9:  tilted [1 0 1]/sqrt(2)
          0.31 0.17 0.09; ...                              % row 10:    generic q
          1/16 0 0; ...                                    % row 11:    16^3 shell, a*
          0    0 1/16 ];                                   % row 12:    16^3 shell, c*

[Vca30, ~, ~, info30] = invz_odd_blocks(ion, q_all, struct('dpRng', dpRng0, 'cache', true));
mabs = @(k) max(abs(Vca30(:,:,k)), [], 'all');
maxca_a  = arrayfun(@(k) mabs(k), 1:3);
maxca_c  = arrayfun(@(k) mabs(k), 4:6);
maxca_t  = arrayfun(@(k) mabs(k), 7:9);
maxca_g  = mabs(10);
maxca_a16 = mabs(11);
maxca_c16 = mabs(12);

fprintf('\n-- (iv) Gamma info (dpRng %d) --\n', dpRng0);
fprintf('Jcc0=%.17g meV  Jaa0=%.17g meV\n', info30.Jcc0, info30.Jaa0);
fprintf('Jcc0_dipole=%.17g meV (published anchor 6.821e-3)  Jaa0_dipole=%.17g meV (3.912e-3)\n', ...
    info30.Jcc0_dipole, info30.Jaa0_dipole);
fprintf('ANCHOR jcc0 = %.17g\n', info30.Jcc0);
fprintf('ANCHOR jaa0 = %.17g\n', info30.Jaa0);

fprintf('\n-- (i) on-axis decay  [meV] --\n');
fprintf(' q        maxJca[q00]        maxJca[00q]\n');
for k = 1:3
    fprintf('%.0e   %.17g   %.17g\n', qs(k), maxca_a(k), maxca_c(k));
end
fprintf('decay ratio maxJca[q00]: (1e-1)/(1e-2)=%.6g  (1e-2)/(1e-3)=%.6g  (linear-in-q => ~10)\n', ...
    maxca_a(1)/maxca_a(2), maxca_a(2)/maxca_a(3));
fprintf('residual maxJca[q00] @ q=1e-3 / Jcc0 = %.6e\n', maxca_a(3)/info30.Jcc0);
fprintf('ANCHOR odd_onaxis_smallq.q     = [%.0e %.0e %.0e]  (along [q 0 0])\n', qs(1), qs(2), qs(3));
fprintf('ANCHOR odd_onaxis_smallq.maxca = [%.17g %.17g %.17g]\n', maxca_a(1), maxca_a(2), maxca_a(3));

fprintf('\n-- (ii) tilted ray q*[1 0 1]/sqrt(2)  [meV] -- NON-decaying limit\n');
for k = 1:3
    fprintf(' q=%.0e:  maxJca=%.17g = %.6g ueV   (ratio to 4pi*gfac/Vc = %.6g)\n', ...
        qs(k), maxca_t(k), maxca_t(k)*1e3, maxca_t(k)/macroscale);
end
fprintf('=> tilted-ray limit does NOT decay: small-q decay assertions stay ON-AXIS ONLY.\n');
fprintf('ANCHOR odd_tilted_limit.q       = %.0e  (along [1 0 1]/sqrt(2))\n', qs(3));
fprintf('ANCHOR odd_tilted_limit.maxca   = %.17g\n', maxca_t(3));
fprintf('ANCHOR odd_tilted_limit.macroscale = %.17g\n', macroscale);
fprintf('ANCHOR odd_tilted_limit.ratio   = %.17g\n', maxca_t(3)/macroscale);

fprintf('\n-- (iii) generic q = [0.31 0.17 0.09]  [meV] --\n');
fprintf('maxJca=%.17g (%.6g ueV, %.6g*Jcc0)\n', maxca_g, maxca_g*1e3, maxca_g/info30.Jcc0);
fprintf('ANCHOR odd_generic_q_max = %.17g\n', maxca_g);

fprintf('\n-- (v) smallest 16^3 shells, dpRng sensitivity  [meV] --\n');
fprintf('dpRng=%d (from the combined call above): maxJca[1/16,0,0]=%.17g (%.6g*Jcc0)   maxJca[0,0,1/16]=%.17g (%.6g*Jcc0)\n', ...
    dpRng0, maxca_a16, maxca_a16/info30.Jcc0, maxca_c16, maxca_c16/info30.Jcc0);

dpSens = [20 40];
sensRows = struct('dpRng', dpRng0, 'maxca_a16', maxca_a16, 'maxca_c16', maxca_c16, ...
    'Jcc0', info30.Jcc0, 'Jaa0', info30.Jaa0);
for dp = dpSens
    q_shell = [1/16 0 0; 0 0 1/16];
    [Vcas, ~, ~, infos] = invz_odd_blocks(ion, q_shell, struct('dpRng', dp, 'cache', true));
    m_a16 = max(abs(Vcas(:,:,1)), [], 'all');
    m_c16 = max(abs(Vcas(:,:,2)), [], 'all');
    fprintf('dpRng=%d: maxJca[1/16,0,0]=%.17g (%.6g*Jcc0)   maxJca[0,0,1/16]=%.17g (%.6g*Jcc0)   Jcc0=%.17g  Jaa0=%.17g\n', ...
        dp, m_a16, m_a16/infos.Jcc0, m_c16, m_c16/infos.Jcc0, infos.Jcc0, infos.Jaa0);
    sensRows(end+1) = struct('dpRng', dp, 'maxca_a16', m_a16, 'maxca_c16', m_c16, ...
        'Jcc0', infos.Jcc0, 'Jaa0', infos.Jaa0); %#ok<AGROW>
end
[~, ord] = sort([sensRows.dpRng]);
sensRows = sensRows(ord);
fprintf('\nANCHOR dpRng_sensitivity.dpRng     = [%s]\n', num2str([sensRows.dpRng]));
fprintf('ANCHOR dpRng_sensitivity.maxca_a16 = [%s]\n', num2str([sensRows.maxca_a16], '%.17g '));
fprintf('ANCHOR dpRng_sensitivity.maxca_c16 = [%s]\n', num2str([sensRows.maxca_c16], '%.17g '));
fprintf('ANCHOR dpRng_sensitivity.Jcc0      = [%s]\n', num2str([sensRows.Jcc0], '%.17g '));
fprintf('ANCHOR dpRng_sensitivity.Jaa0      = [%s]\n', num2str([sensRows.Jaa0], '%.17g '));
rel_a16 = (sensRows(1).maxca_a16 - sensRows(3).maxca_a16) / sensRows(3).maxca_a16;
rel_c16 = (sensRows(1).maxca_c16 - sensRows(3).maxca_c16) / sensRows(3).maxca_c16;
fprintf('relative spread (dpRng 20 vs 40): a16 = %.6g   c16 = %.6g\n', rel_a16, rel_c16);

fprintf('\n=== explore_tensor_blocks done ===\n');

% ------- local helper -------
function h = git_hash()
[status, out] = system('git rev-parse HEAD');
if status == 0, h = strtrim(out); else, h = 'unknown'; end
end
