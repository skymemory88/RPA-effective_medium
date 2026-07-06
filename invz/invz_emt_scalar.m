function med = invz_emt_scalar(G0, Sigma, Jnu_flat, opts)
%INVZ_EMT_SCALAR Scalar effective-medium fixed point at fixed self-energy.
% G  = G0./(1 + Sigma + K.*G0)            (Dyson, R eq 9)
% Gq = G ./ (1 + (J_nu - K).*G)           (R eq 7)
% K  = (1/N) sum_{q,nu} J_nu.*Gq ./ G     (R eq 8)
if nargin < 4, opts = struct(); end
tol  = 1e-10; if isfield(opts,'tol'),      tol  = opts.tol;      end
mit  = 500;   if isfield(opts,'max_iter'), mit  = opts.max_iter; end
mix  = 0.5;   if isfield(opts,'mix'),      mix  = opts.mix;      end
G0 = G0(:);  Sigma = Sigma(:);  Jf = Jnu_flat(:);
K = zeros(size(G0)); if isfield(opts,'K0') && ~isempty(opts.K0), K = opts.K0(:); end
nJ = numel(Jf);
converged = false;
for it = 1:mit
    G  = G0 ./ (1 + Sigma + K.*G0);
    Gq = G.' ./ (1 + (Jf - K.').*G.');           % [nJ, nw]
    Kn = (Jf.'*Gq).'/nJ ./ G;
    res = max(abs(Kn - K) ./ max(abs(Kn), 1e-14));
    K = K + mix*(Kn - K);
    if res < tol, converged = true; break; end
end
G  = G0 ./ (1 + Sigma + K.*G0);
Gq = G.' ./ (1 + (Jf - K.').*G.');
med.G = G;  med.K = K;
med.closure = max(abs(mean(Gq,1).' - G) ./ max(abs(G), 1e-14));
med.converged = converged;  med.iters = it;
end
