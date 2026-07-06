function sig = invz_sigma(tl, lam, K, g, beta)
%INVZ_SIGMA Jensen self-energy Sigma = alpha + gamma.*g  (HTML eqs 22, 23, 27).
pref = tl.M2 / tl.n01^2;
sig.alpha = pref * (lam(2) - 0.5*(tl.g0 + beta*(1 - tl.n01^2))*lam(1));
sig.gamma = pref * (lam(1) - (1 - tl.n01^2)*K(:));
sig.Sigma = sig.alpha + sig.gamma .* g(:);
end
