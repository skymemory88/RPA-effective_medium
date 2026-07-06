function [wn, wts, beta] = invz_matsubara(T, Ecut)
%INVZ_MATSUBARA Bosonic Matsubara grid up to Ecut (meV), n>=0 with doubling weights.
C = invz_const();  beta = 1/(C.kB*T);
nmax = max(8, ceil(Ecut*beta/(2*pi)));
wn  = 2*pi*(0:nmax).'/beta;
wts = [1; 2*ones(nmax,1)];
end
