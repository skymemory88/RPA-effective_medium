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
Sigma = zeros(size(wn));
for outer=1:200
  med = invz_emt_scalar(G0, Sigma, Jnu, struct());
  if ~med.converged, break; end
  lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
  sg = invz_sigma(tl, lam, med.K, g, beta);
  dS = max(abs(sg.Sigma-Sigma)); Sigma = Sigma + 0.7*(sg.Sigma-Sigma);
  if dS<1e-8, break; end
end
med = invz_emt_scalar(G0, Sigma, Jnu, struct());
lam = invz_lambdas(med.K, g, wts, beta, [1 2]);
sg = invz_sigma(tl, lam, med.K, g, beta);
K = med.K;
G = G0 ./ (1 + Sigma + K.*G0);
alpha = sg.alpha;
fprintf('alpha = %.6e, Sigma(end) [large wn] = %.6e (should -> alpha)\n', alpha, Sigma(end));
fprintf('Sigma(1) [wn=0] = %.6e\n', Sigma(1));

n1 = (1-tl.n01)/2;
dU37 = 0.5*( sg.alpha*tl.n01*tl.Delta/(1+sg.alpha) - tl.M2*lam(1) + real(sum(wts.*K.*(G-G0)))/beta );

% naive treatment (no alpha-tail correction): drop odd part of (Delta+wn)(G-G0) wholesale
dU21_naive = 0.5*real(sum(wts.*K.*G))/beta - (tl.Delta/(2*beta*tl.M2))*real(sum(wts.*(G-G0)));

% alpha-tail-corrected treatment: R = (G-G0) + alpha*G0 (faster-decaying remainder)
R = (G - G0) + alpha*G0;
fprintf('max|G-G0| = %.4e, max|R| = %.4e (ratio %.4g -- R should be MUCH smaller if alpha-tail dominates)\n', ...
    max(abs(G-G0)), max(abs(R)), max(abs(R))/max(abs(G-G0)));
fprintf('G-G0(end) [large wn] = %.4e, -alpha*G0(end) = %.4e, R(end) = %.4e\n', ...
    G(end)-G0(end), -alpha*G0(end), R(end));

dU21_corrected = 0.5*real(sum(wts.*K.*G))/beta - alpha*n1*tl.Delta - (tl.Delta/(2*beta*tl.M2))*real(sum(wts.*R));

fprintf('\ndU37 (eq37)                 = %.6e\n', dU37);
fprintf('dU21_naive (no tail fix)    = %.6e  (rel diff %.4g%%)\n', dU21_naive, 100*abs(dU37-dU21_naive)/max(abs(dU37),abs(dU21_naive)));
fprintf('dU21_corrected (tail fix)   = %.6e  (rel diff %.4g%%)\n', dU21_corrected, 100*abs(dU37-dU21_corrected)/max(abs(dU37),abs(dU21_corrected)));
