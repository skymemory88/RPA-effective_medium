function tf = invz_is_gamma_equiv(q, tau)
%INVZ_IS_GAMMA_EQUIV True where q (r.l.u.) is a Gamma-equivalent point of the lattice.
% Shared by invz_jq_modes (Lorentz cavity placement) and invz_jq_path (Gamma-limit guard).
tf = abs(real(sum(exp(2i*pi*(tau*q.'))))/size(tau,1) - 1) < 1e-9;
end
