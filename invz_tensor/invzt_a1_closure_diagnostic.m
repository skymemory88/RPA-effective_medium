function out = invzt_a1_closure_diagnostic(cdom, crest, Sigma, lat, opts)
%INVZT_A1_CLOSURE_DIAGNOSTIC Side-effect-free split/whole A1 medium comparison.
%   OUT = INVZT_A1_CLOSURE_DIAGNOSTIC(CDOM, CREST, SIGMA, LAT, OPTS) evaluates
%   one frozen A1 medium without updating SIGMA. CDOM and CREST must be the
%   byte-identical transition split supplied to both candidates.
%
%   opts.closure  'split' (default): CDOM/(1+Sigma) + CREST
%                 'whole'          : (CDOM+CREST)/(1+Sigma)
%   opts.odd      true (default): false applies INVZT_ODD_MASK to a copy of LAT
%   opts.rank_tol 1e-12: static active-subspace clipping threshold
%
%   This is a diagnostic primitive, not a second solver and not a production
%   representation selector. It returns the full frozen chi_tilde plus the
%   BZ-local medium, strict-Gamma mass, rank/clipped mass, and uniform response.
if nargin < 5, opts = struct(); end
closure = char(getf(opts, 'closure', 'split'));
odd = ~isfield(opts, 'odd') || isempty(opts.odd) || ~isequal(opts.odd, false);
rank_tol = getf(opts, 'rank_tol', 1e-12);

if ~ismember(closure, {'split', 'whole'})
    error('invzt:closureDiagnosticMode', ...
        'opts.closure must be ''split'' or ''whole''.');
end
if ~(isnumeric(cdom) && isnumeric(crest) && isequal(size(cdom), size(crest)) ...
        && size(cdom,1) == 3 && size(cdom,2) == 3)
    error('invzt:closureDiagnosticShape', ...
        'cdom and crest must be equal-sized numeric [3,3,nwn] arrays.');
end
nwn = size(cdom, 3);
if ~(isnumeric(Sigma) && isreal(Sigma) && numel(Sigma) == nwn ...
        && all(isfinite(Sigma(:))))
    error('invzt:closureDiagnosticSigma', ...
        'Sigma must contain %d finite real entries.', nwn);
end
Sigma = Sigma(:);
if any((1 + Sigma) <= 0)
    error('invzt:closureDiagnosticDomain', ...
        'The frozen state crosses the required 1+Sigma > 0 domain.');
end
if ~(isstruct(lat) && isfield(lat, 'Jt') && isfield(lat, 'JtGamma') ...
        && isfield(lat, 'w') && isfield(lat, 'info') ...
        && isfield(lat.info, 'Jaa0') && isfield(lat.info, 'Jcc0'))
    error('invzt:closureDiagnosticLattice', ...
        'lat must be a full INVZT_JQ_TENSOR result with Jt/JtGamma/w/info.');
end

denom = reshape(1 + Sigma, 1, 1, nwn);
switch closure
    case 'split'
        ctil = cdom ./ denom + crest;
    case 'whole'
        ctil = (cdom + crest) ./ denom;
end

lat_eff = lat;
if ~odd, lat_eff.Jt = invzt_odd_mask(lat.Jt); end
[Gcc, diag4] = invzt_gcc_lattice(ctil, lat_eff);
Gloc = -Gcc(:);
G0til = -real(squeeze(ctil(3,3,:)));
K = 1 ./ Gloc - 1 ./ G0til;

ctil0 = real((ctil(:,:,1) + ctil(:,:,1)')/2);
[crit, clipped_mass, active_rank] = ...
    invzt_crit_static(ctil0, lat.JtGamma, rank_tol);
C12 = kron(eye(4), ctil0);
[Uc, Dc] = eig((C12 + C12')/2);
dc = real(diag(Dc));
dc(dc < rank_tol) = 0;
Sc = Uc * diag(sqrt(max(dc, 0))) * Uc';
Mc = eye(12) - Sc*lat.JtGamma*Sc;
gamma_mass_eigenvalues = sort(real(eig((Mc + Mc')/2)));
if abs(gamma_mass_eigenvalues(1) - crit) > 100*eps(max(1, abs(crit)))
    error('invzt:closureDiagnosticCritMismatch', ...
        'Reconstructed Gamma mass spectrum does not reproduce INVZT_CRIT_STATIC.');
end

Jd = kron(ones(4)/4, diag([lat.info.Jaa0, lat.info.Jaa0, lat.info.Jcc0]));
u = kron(ones(4,1)/2, eye(3));
Xu = invzt_chi_rpa(ctil0, Jd);
Xg = invzt_chi_rpa(ctil0, lat.JtGamma);
chi_uniform = u' * Xu * u;
chi_gamma_projected = u' * Xg * u;
full0 = cdom(:,:,1) + crest(:,:,1);

out = struct();
out.schema = 'invzt_a1_closure_diagnostic/v1';
out.closure = closure;
out.odd = odd;
out.rank_tol = rank_tol;
out.chi_tilde = ctil;
out.chi_tilde0 = ctil0;
out.local_static_eigenvalues = sort(real(eig(ctil0)));
out.G = Gloc;
out.Gloc0 = Gloc(1);
out.K = K;
out.K0 = K(1);
out.diag4 = diag4;
out.diag4_spread = max(max(diag4, [], 1) - min(diag4, [], 1));
out.crit = crit;
out.gamma_mass_eigenvalues = gamma_mass_eigenvalues;
out.crit_clipped_mass = clipped_mass;
out.crit_active_rank = active_rank;
out.chi_uniform0 = real((chi_uniform + chi_uniform')/2);
out.chi_gamma_projected0 = real((chi_gamma_projected + chi_gamma_projected')/2);
out.uniform_static_eigenvalues = sort(real(eig(out.chi_uniform0)));
out.boundary_handoff_ratio = real(ctil0(3,3)) / real(full0(3,3));
out.min_one_plus_sigma = min(1 + Sigma);
end

function v = getf(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
