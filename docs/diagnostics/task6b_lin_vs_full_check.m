here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'invz_common'));
addpath(fullfile(here, '..', '..', 'invz_projected'));
T = 0.31; C = invz_const(); beta = 1/(C.kB*T);
D0 = 0.2; M0op = 3.0;
Jnu = linspace(-4e-3, 4.0e-3, 24).';
Jz2 = [0, M0op; M0op, 0]; H = [D0/2 0; 0 -D0/2];
[V,E] = eig((H+H')/2,'vector'); [E,ix]=sort(real(E)); V=V(:,ix);
p = exp(-beta*(E-E(1))); p=p/sum(p);
Mz = V'*Jz2*V;
tl = struct('m', real(Mz(1,1)), 'M2', abs(Mz(1,2))^2, 'n01', p(1)-p(2), 'Delta', E(2)-E(1), 'g0', 0);
tl.g0 = 2*tl.n01/tl.Delta;
[wn,wts,~] = invz_matsubara(T,40);
g = real(invz_g(tl, 1i*wn));
G0 = -(tl.M2*g);
Sigma = zeros(size(wn)); ok=false;
for outer=1:200
  med = invz_emt_scalar(G0, Sigma, Jnu, struct());
  if ~med.converged, break; end
  lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
  sg = invz_sigma(tl, lam, med.K, g, beta);
  dS = max(abs(sg.Sigma-Sigma)); Sigma = Sigma + 0.7*(sg.Sigma-Sigma);
  if dS<1e-8, ok=true; break; end
end
fprintf('converged: %d\n', ok);
med = invz_emt_scalar(G0, Sigma, Jnu, struct());
lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
sg = invz_sigma(tl, lam, med.K, g, beta);
K = med.K;
G_full = G0 ./ (1 + Sigma + K.*G0);
G_lin  = G0 - G0.*(Sigma + K.*G0);
dU37 = 0.5*( sg.alpha*tl.n01*tl.Delta/(1+sg.alpha) - tl.M2*lam(1) + real(sum(wts.*K.*(G_full-G0)))/beta );
dU21_full = 0.5*real(sum(wts.*K.*G_full))/beta - (tl.Delta/(2*beta*tl.M2))*real(sum(wts.*(G_full-G0)));
dU21_lin  = 0.5*real(sum(wts.*K.*G_lin))/beta  - (tl.Delta/(2*beta*tl.M2))*real(sum(wts.*(G_lin-G0)));
fprintf('dU37 (eq37, full G in K(G-G0) term) = %.6e\n', dU37);
fprintf('dU21 with FULL resummed G           = %.6e\n', dU21_full);
fprintf('dU21 with LINEARIZED (Born) G       = %.6e\n', dU21_lin);
fprintf('rel diff dU37 vs dU21_lin  = %.4g%%\n', 100*abs(dU37-dU21_lin)/max(abs(dU37),abs(dU21_lin)));
fprintf('rel diff dU37 vs dU21_full = %.4g%%\n', 100*abs(dU37-dU21_full)/max(abs(dU37),abs(dU21_full)));
maxG = max(abs(G_full-G0)); maxGlin = max(abs(G_lin-G0));
fprintf('max|G_full-G0| = %.4e, max|G_lin-G0| = %.4e (ratio %.4g)\n', maxG, maxGlin, maxG/maxGlin);
fprintf('alpha = %.6e, lam1 = %.6e, lam2 = %.6e\n', sg.alpha, lam(1), lam(2));
