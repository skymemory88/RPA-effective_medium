function diag4()
% Anatomy of the outer Sigma<->EMT map at omega_n = 0 in the ordered region.
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag4.log'); if exist(LOG,'file'), delete(LOG); end

ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:);  J0 = info.Jcc0;
Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end
T = 0.31;  Ecut = 40;
[wn, wts, beta] = invz_matsubara(T, Ecut);

for B = [1.0 3.0 3.6 4.0 4.4]
    Bv = [B 0 0];
    si = invz_single_ion(ion, T, Bv, struct('hyp',true,'Jxx0',Jaa0,'transverse_mf','legacy_x'));
    c0 = invz_chi0z(si, T, 1i*wn, struct('elastic', true));
    G0 = -real(squeeze(c0(3,3,:)));
    tl = invz_twolevel(ion, T, Bv, struct('Jxx0',Jaa0,'transverse_mf','legacy_x'));
    g  = real(invz_g(tl, 1i*wn));
    chi0 = -G0(1);
    say(LOG, sprintf('===== B = %.2f T ==================================================', B));
    say(LOG, sprintf(' chi0cc(0) = %.6g meV^-1 ; Delta = %.5g meV ; M2 = %.5g ; n01 = %.6g ; g0 = %.5g', ...
        chi0, tl.Delta, tl.M2, tl.n01, tl.g0));
    say(LOG, sprintf(' MF/RPA:  1 - J0*chi0 = %+.5g   ;  1 - max(Jq)*chi0 = %+.5g', 1-J0*chi0, 1-max(Jf)*chi0));
    say(LOG, sprintf(' Sigma0 that makes the WORST grid mode marginal:  Sigma0* = max(Jq)*chi0 - 1 = %+.6g', max(Jf)*chi0-1));

    % --- scan Sigma0 (holding the rest of Sigma at its converged-ish shape = 0) --------------
    say(LOG,'   Sigma0      #(D+J*G0<0)      A(0)          K(0)        lam1        lam2      F(Sigma)(0)');
    for s0 = [-0.5 -0.2 0 0.05 0.1 0.2 0.5 1.0 2.0 max(Jf)*chi0-1-0.01 max(Jf)*chi0-1+0.01]
        Sig = zeros(size(wn));  Sig(1) = s0;
        den = 1 + s0 + Jf*G0(1);
        nneg = nnz(den < 0);
        med = invz_emt_scalar(G0, Sig, Jf, struct());
        A0 = mean(Jf./den);
        lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
        sg = invz_sigma(tl, lam, med.K, g, beta);
        say(LOG, sprintf('%10.5g %10d %14.6g %13.6g %11.4g %11.4g %13.6g', ...
            s0, nneg, A0, med.K(1), lam(1), lam(2), sg.Sigma(1)));
    end

    % --- the actual PM iteration trace -------------------------------------------------------
    Sig = zeros(size(wn));  hist = nan(1,60);
    say(LOG, '  --- PM outer iteration (mix 0.7) : iter, dS, Sigma0, K(0), #neg-denominator ---');
    for it = 1:60
        med = invz_emt_scalar(G0, Sig, Jf, struct());
        lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
        sg  = invz_sigma(tl, lam, med.K, g, beta);
        dS  = max(abs(sg.Sigma - Sig));
        nneg = nnz(1 + Sig(1) + Jf*G0(1) < 0);
        hist(it) = Sig(1);
        if it <= 30 || mod(it,10)==0
            say(LOG, sprintf('  %3d  dS=%11.4g  Sigma0=%12.6g  K0=%12.5g  nneg=%6d', it, dS, Sig(1), med.K(1), nneg));
        end
        Sig = Sig + 0.7*(sg.Sigma - Sig);
        if dS < 1e-8, say(LOG,sprintf('  CONVERGED at iter %d', it)); break; end
    end
    say(LOG,'');
end
say(LOG,'done');
end
function say(LOG,s), fid=fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid); end
