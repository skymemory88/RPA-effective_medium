function Sc = invz_sigma_crit(J0, Jnu_flat)
%INVZ_SIGMA_CRIT Closed-form critical self-energy, Rønnow 2007 eq (10):
%   Sigma_c(0) = (1/N) sum_{q,nu} J_nu(q) / (J(0) - J_nu(q)),
% valid at the zero-frequency critical point of the degenerate (elastic) doublet.
x = Jnu_flat(:);
keep = (J0 - x) > 1e-12;
if ~all(keep)
    warning('invz:sigmaCritExcluded', 'Excluded %d uniform-mode entries.', sum(~keep));
end
Sc = mean(x(keep) ./ (J0 - x(keep)));
end
