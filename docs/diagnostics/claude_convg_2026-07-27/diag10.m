function diag10()
% Verify the Stieltjes identities that section 10 rests on:
%   (1) A(w) = (1/G0)*[1 + x*S(x)],  x = -(1+Sigma)/G0,  S(x) = mean(1./(J-x))
%   (2) ordered static: mean_q Gq = S(y), K0 = 1/S(y)+y,  y = K0 - 1/Gstat
%   (3) x crosses max(J) exactly at the measured resummed window edge
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root); addpath(fullfile(root,'invz_projected')); addpath(fullfile(root,'invz_common'));
warning('off','MATLAB:nearlySingularMatrix');
SD = '/private/tmp/claude-501/-Users-yikaiyang-Library-CloudStorage-OneDrive-Nexus365-Programming-scripts-Matlab-Simulation-invZ-expansion/fe404c1e-1dd6-4ec7-868e-75af95bf66a2/scratchpad';
LOG = fullfile(SD,'diag10.log'); if exist(LOG,'file'), delete(LOG); end
ion = invz_ion();
[Jnu, info] = invz_bz_couplings(ion, struct('grid',[16 16 16],'dpRng',30,'cache',true));
Jf = Jnu(:); Jaa0 = ion.Jxx0; if isfield(info,'Jaa0'), Jaa0 = info.Jaa0; end
T = 0.31; [wn,~,~] = invz_matsubara(T,40);

% ---- (1) dynamic sector identity, at a real converged state -------------------------
p = invz_solve_point(ion, T, [4.4 0 0], Jf, struct('J0eff',info.Jcc0,'Jxx0',Jaa0));
si = invz_single_ion(ion,T,[4.4 0 0],struct('hyp',true,'Jxx0',Jaa0,'transverse_mf','legacy_x'));
c0 = invz_chi0z(si,T,1i*wn,struct('elastic',true));  G0 = -real(squeeze(c0(3,3,:)));
Sig = p.Sigma(:);  D = 1 + Sig;
med = invz_emt_scalar(G0, Sig, Jf, struct());
A_direct = zeros(size(D)); A_stielt = zeros(size(D));
for n = 1:numel(D)
    A_direct(n) = mean(Jf./(D(n) + Jf*G0(n)));
    x = -D(n)/G0(n);
    S = mean(1./(Jf - x));
    A_stielt(n) = (1/G0(n))*(1 + x*S);
end
say(LOG,'--- (1) A = (1/G0)[1 + x S(x)] vs direct mean, over all Matsubara slots (B=4.4 T, T=0.31 K) ---');
say(LOG,sprintf('  max abs diff = %.3e   max rel diff = %.3e', ...
    max(abs(A_direct-A_stielt)), max(abs(A_direct-A_stielt)./max(abs(A_direct),realmin))));
% and K rebuilt from S only
K_stielt = A_stielt.*D./(1 - A_stielt.*G0);
say(LOG,sprintf('  K from S vs med.K : max rel diff = %.3e', max(abs(K_stielt-med.K)./max(abs(med.K),realmin))));

% ---- (2) ordered static sector identity ---------------------------------------------
pj = invz_solve_point_ordered(ion,T,[4.0 0 0],Jf,struct('J0eff',info.Jcc0,'Jxx0',Jaa0,'ordered_mode','jensen'));
Gstat = pj.G(1);  K0 = pj.K(1);
y  = K0 - 1/Gstat;
Sy = mean(1./(Jf - y));
Gq = Gstat./(1 + (Jf-K0).*Gstat);
say(LOG,'');
say(LOG,'--- (2) ordered static: mean_q Gq  ==  S(y),  y = K0 - 1/Gstat  (B=4.0 T jensen state) ---');
say(LOG,sprintf('  mean(Gq)    = %.12g', mean(Gq)));
say(LOG,sprintf('  S(y)        = %.12g   (rel diff %.3e)', Sy, abs(Sy-mean(Gq))/abs(mean(Gq))));
say(LOG,sprintf('  Gstat       = %.12g   closure |S(y)-Gstat|/|Gstat| = %.3e', Gstat, abs(Sy-Gstat)/abs(Gstat)));
say(LOG,sprintf('  1/S(y)+y    = %.12g   vs K0 = %.12g  (rel diff %.3e)', 1/Sy+y, K0, abs(1/Sy+y-K0)/max(abs(K0),realmin)));

% ---- (3) where does x = (1+Sigma0)/chi0 cross max(J)? --------------------------------
say(LOG,'');
say(LOG,sprintf('--- (3) x = (1+Sigma(0))/chi0cc(0) vs max_q J = %.6g ueV ---', 1e3*max(Jf)));
say(LOG,'    B      Sigma0       chi0cc(0)      x (ueV)    x-maxJ (ueV)  outside?');
for B = [3.0 3.4 3.6 3.8 4.0 4.2 4.4]
    q = invz_solve_point(ion, T, [B 0 0], Jf, struct('J0eff',info.Jcc0,'Jxx0',Jaa0));
    s2 = invz_single_ion(ion,T,[B 0 0],struct('hyp',true,'Jxx0',Jaa0,'transverse_mf','legacy_x'));
    cc = invz_chi0z(s2,T,0,struct('elastic',true));  chi0 = real(cc(3,3,1));
    x = (1+q.Sigma0)/chi0;
    say(LOG,sprintf('  %4.2f  %10.5g  %12.6g  %10.5g  %12.4g   %s   (PMconv=%d)', ...
        B, q.Sigma0, chi0, 1e3*x, 1e3*(x-max(Jf)), tf(x>max(Jf)), q.converged));
end
say(LOG,'done');
end
function s = tf(b), if b, s='YES'; else, s='no '; end, end
function say(LOG,s), fid=fopen(LOG,'a'); fprintf(fid,'%s\n',s); fclose(fid); end
