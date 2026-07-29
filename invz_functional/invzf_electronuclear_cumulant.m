function out = invzf_electronuclear_cumulant(ion, T, B, hz, labels, opts)
%INVZF_ELECTRONUCLEAR_CUMULANT Exact connected Jz cumulant of the full 136-state local ion.
%
%   OUT = INVZF_ELECTRONUCLEAR_CUMULANT(ION, T, B, HZ, LABELS, OPTS) returns the exact
%   connected rank-2/3/4 Matsubara cumulant of Jz for the source-biased single-ion
%   problem, i.e. the local generator of invzp_convg_fix.md WP1 evaluated on all 136
%   electronuclear states rather than on a two-level projection.
%
%   This is a WIRING function, not new physics: invzf_multilevel_cumulant is already an
%   exact finite-level cumulant engine for any (E, O, beta), and invz_single_ion already
%   diagonalises the full electronuclear Hamiltonian and exports the eigenvalues E and
%   the Jz matrix Mz in that same eigenbasis. Before this function existed the two were
%   never connected, so C3 and C4 for the 136-state problem were available only as
%   finite-difference h-derivatives of C2 (invzf_electronuclear_local), never directly.
%
%   THE SOURCE. HZ (meV) is imposed through invz_single_ion's hz_fixed mode, which adds a
%   FIXED -HZ*Jz term. OPTS.transverse_mf defaults to 'none' here, unlike production:
%   the object WP1 needs is the LOCAL generator, so the transverse mean field must not be
%   re-solved as HZ varies. Solving it would make every h-derivative a mean-field-dressed
%   derivative and would silently break the static-derivative identity that grades this
%   oracle. Pass 'legacy_x' deliberately if a dressed generator is what you want.
%
%   LABELS are bosonic Matsubara INDICES (not frequencies) and must sum to zero; their
%   count selects the rank. OPTS is forwarded to invzf_multilevel_cumulant (local_rank,
%   degeneracy_tolerance, reality_tolerance, ...), plus:
%     .hyp          (true)     include the I = 7/2 nuclear space (dim 136)
%     .transverse_mf('none')   see above
%
%   OUT: .value, .cumulant (the engine's full output struct), .m = <Jz>, .beta,
%        .E, .nlevels, .discarded_boltzmann_weight, .si (the single-ion solution).
%
%   Truncation caveat inherited from the engine: discarded_boltzmann_weight reports only
%   omitted thermal weight and is NOT a bound on virtual intermediate-state truncation.
%   A local_rank ladder is mandatory before any truncated result is used.
if nargin < 6 || isempty(opts), opts = struct(); end
validateattributes(hz, {'numeric'}, {'real','scalar','finite'});

hyp = getf_local(opts, 'hyp', true);
tmf = getf_local(opts, 'transverse_mf', 'none');

sopts = struct('hyp', hyp, 'transverse_mf', tmf, 'hz_fixed', hz);
si = invz_single_ion(ion, T, B, sopts);

kB = 8.617333262e-2;                       % meV/K, matching invz_const's convention
beta = 1/(kB*T);

eopts = rmfield_if(opts, {'hyp','transverse_mf','allow_bad_status'});
cum = invzf_multilevel_cumulant(si.E, si.Mz, beta, labels, eopts);

% MEASURED DEFECT (2026-07-29, WP1 gate). invzf_multilevel_cumulant's default
% dense_beta_span_limit = 2000 selects the dense block-exponential backend at
% beta*energy_span up to 2000, but that backend returns status = 'nonfinite' (a NaN
% cumulant) for RANK-2 cumulants at NONZERO Matsubara labels once beta*span exceeds
% roughly 1e3. Reproduction: LiHoF4 ion, B = [4 0 0] T, hz = 0.02 meV, labels [3 -3];
% ok at T >= 0.6 K (beta*span 813) and nonfinite at T = 0.3 K (beta*span 1627); the
% 17-state electronic-only problem fails identically at the same span, so it is a span
% effect, not a dimension effect. Rank 2 at zero frequency and rank 3 at nonzero
% frequency are unaffected. Setting dense_beta_span_limit below beta*span selects the
% exponential-action backend, which returns a finite and limit-independent value.
% A NaN must never leave this function silently: a caller that does not inspect
% .status would propagate it into a functional coefficient. Pass
% opts.allow_bad_status = true to receive the raw result instead of an error.
if ~strcmp(cum.status, 'ok') && ~getf_local(opts, 'allow_bad_status', false)
    error('invzf:electronuclearCumulantStatus', ...
        ['invzf_multilevel_cumulant returned status ''%s'' (backend ''%s'', ' ...
         'beta*energy_span = %.4g, dense_beta_span_limit = %.4g) for rank %d labels %s. ' ...
         'The cumulant is not usable. Selecting the exponential-action backend by ' ...
         'setting opts.dense_beta_span_limit below beta*energy_span is the known ' ...
         'workaround; opts.allow_bad_status = true returns the raw result unchecked.'], ...
        cum.status, char(string(cum.evaluation_backend)), cum.beta_energy_span, ...
        cum.dense_beta_span_limit, numel(labels), mat2str(labels(:).'));
end

% invzf_multilevel_cumulant returns the connected cumulant in .connected (schema v2);
% .value does not exist. Re-exported as .value here so every caller reads one name.
out = struct('value', cum.connected, 'cumulant', cum, 'm', si.Jexp(3), ...
    'beta', beta, 'E', si.E, 'nlevels', numel(si.E), 'hz', hz, ...
    'labels', labels(:).', 'transverse_mf', tmf, 'si', si);
if isfield(cum, 'discarded_boltzmann_weight')
    out.discarded_boltzmann_weight = cum.discarded_boltzmann_weight;
else
    out.discarded_boltzmann_weight = NaN;
end
end

function v = getf_local(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = rmfield_if(s, names)
for k = 1:numel(names)
    if isfield(s, names{k}), s = rmfield(s, names{k}); end
end
end
