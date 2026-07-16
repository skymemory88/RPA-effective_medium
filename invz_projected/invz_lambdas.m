function lam = invz_lambdas(K, g, wts, beta, plist)
%INVZ_LAMBDAS lambda_p = (1/beta) sum_n K(iwn) g(iwn)^p   (HTML eq 21 / J 2.19).
K = K(:); g = g(:); wts = wts(:);
lam = zeros(numel(plist),1);
for ip = 1:numel(plist)
    lam(ip) = real(sum(wts .* K .* g.^plist(ip))) / beta;
end
end
