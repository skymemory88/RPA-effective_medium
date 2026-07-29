function out = invz_reciprocal_static_closure(go, Jf, K0)
%INVZ_RECIPROCAL_STATIC_CLOSURE Static EMT q-average in the reciprocal (defactored) coordinate.
%
% ONE definition of the ordered static closure's q-average, written so that the
% removable singularity at the local Gstat pole is evaluated instead of divided
% into Inf/Inf. Both the audit path (INVZ_ORDERED_RESIDUAL block B under
% eso.audit_coordinate = 'defactored') and the production iteration
% (INVZ_EMT_STATIC_ORDERED under eso.closure_coordinate = 'defactored') call this,
% so the two can never drift apart.
%
% This is an EXACT ALGEBRAIC REASSOCIATION of
%       Gq   = Gstat ./ (1 + (J(q) - K0).*Gstat)
%       Gbar = mean(Gq),   Jloc = mean(J.*Gq)/mean(Gq)
% obtained by dividing numerator and denominator by Gstat:
%       Gq   = 1 ./ (1/Gstat + J(q) - K0),      z := 1/Gstat.
% NOT a regulariser: no broadening, no added tolerance, no clipped denominator,
% no sign change, no pole floor. Away from the pole the two forms agree to
% rounding; at the pole the factored form evaluates Inf/Inf = NaN while this one
% evaluates the finite limit Gq -> 1/(J(q) - K0).
%
% z is built from Gstat's own parts rather than as 1/Gstat, so no infinity is
% formed at any stage. With
%       Gstat = G0inel0/d0 + xi*G0el0,   d0 = gstat_local_denom = 1 + Sigma0 + K0*G0inel0,
% one has Gstat = (G0inel0 + xi*G0el0*d0)/d0, hence
%       z = 1/Gstat = d0 / (G0inel0 + xi*G0el0*d0),
% which is exactly 0 at the pole d0 = 0 and loses no digits near it.
%
% Inputs
%   go  the second output of INVZ_GSTAT_ORDERED (needs .gstat_local_denom, .xi,
%       .G0bare and the caller's G0inel0/G0el0 split recovered from it -- see below)
%   Jf  flat coupling vector J(q) (column or row; used as Jf(:))
%   K0  current static medium constant
%
% The G0inel0/G0el0 split is NOT recoverable from `go` alone (go.G0bare is their
% sum), so callers pass it explicitly through go.G0inel0/go.G0el0; INVZ_GSTAT_ORDERED
% records both.
%
% Output struct:
%   .z        1/Gstat in the stable representation (0 at the local pole, Inf when Gstat = 0)
%   .Gbar     mean_q G(q,0)
%   .Jloc     mean_q(J*G)/mean_q(G) -- the closure's K0 update map
%   .Jscale   max|J| (the scale the raw residual is measured against)
% G = -chi (meV^-1), ferromagnetic positive J.
Jf = Jf(:);
d0 = go.gstat_local_denom;
Hz = go.G0inel0 + go.xi*go.G0el0*d0;
z  = d0/Hz;
Jscale = max(abs(Jf));
if isinf(z)
    % Gstat = 0: every Gq = 0, so the weighted mean degenerates to the plain
    % coupling mean (the limit of Jloc as Gstat -> 0), and Gbar is exactly 0.
    Jloc = mean(Jf);
    Gbar = 0;
elseif isfinite(z)
    % `scale` only keeps the intermediate weights in range; it cancels exactly
    % between Gbar and Jloc and changes no value.
    scale = max(abs(z), Jscale);
    weights = scale./(z + Jf - K0);
    meanWeights = mean(weights);
    Gbar = meanWeights/scale;
    Jloc = mean(Jf.*weights)/meanWeights;
else
    Gbar = NaN;
    Jloc = NaN;
end
out = struct('z', z, 'Gbar', Gbar, 'Jloc', Jloc, 'Jscale', Jscale);
end
