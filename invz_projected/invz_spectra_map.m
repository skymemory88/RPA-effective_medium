function S = invz_spectra_map(ion, T, fields, w, opts)
%INVZ_SPECTRA_MAP Uniform-mode chi''_cc over a transverse-field sweep.
% The 1/z leg uses the Jensen-consistent ordered/PM dispatch. The RPA leg is
% phased independently from its own bare-PM mass and bare ordered mean field;
% it never reuses the 1/z-selected state.
% S.phase / S.phase_rpa use 1=ordered, 2=paramagnet, 0=masked.
% S.Bc (alias S.Bc_1z) and S.Bc_rpa are separate sweep-midpoint boundaries;
% S.rpa_mass_pm is 1-Jcc0*chi0cc0 on the independently computed PM state.
%
% opts:
%   hyp       true
%   grid      [16 16 16]
%   dpRng     30
%   cache     false
%   eta       5e-3 meV
%   parallel  false
%   peak_wmin 0 meV
%   verbose   true
%   solve_opts struct passed to the point solvers
%   Jnu/info optional precomputed coupling pair
% Coupling backend/grid fields are forwarded to invz_bz_couplings.
if nargin < 5, opts = struct(); end

hyp      = getf(opts, 'hyp', true);
eta      = getf(opts, 'eta', 5e-3);
parallel = getf(opts, 'parallel', false);
wmin     = getf(opts, 'peak_wmin', 0);
verbose  = getf(opts, 'verbose', true);
sxtra    = getf(opts, 'solve_opts', struct());
invz_check_solve_opts(sxtra);

if ~(isnumeric(fields) && isreal(fields) && all(isfinite(fields(:))) && all(fields(:) >= 0))
    error('invz:fields', 'fields must contain finite nonnegative transverse fields.');
end
fields = fields(:).';
w = w(:);
nB = numel(fields);
nw = numel(w);

hasJnu = isfield(opts, 'Jnu');
hasInfo = isfield(opts, 'info');
if xor(hasJnu, hasInfo)
    error('invz:couplings', 'opts.Jnu and opts.info must be supplied together.');
end
if hasJnu
    Jnu = opts.Jnu(:);
    info = opts.info;
else
    [Jnu, info] = invz_bz_couplings(ion, opts);
end
Jcc0 = info.Jcc0;
Jaa0 = ion.Jxx0;
if isfield(info, 'Jaa0'), Jaa0 = info.Jaa0; end
Jshape = 0;
if isfield(info, 'Jshape_cc'), Jshape = info.Jshape_cc; end

sopts = sxtra;
sopts.hyp = hyp;
sopts.J0eff = Jcc0;
sopts.Jxx0 = Jaa0;

chiz = nan(nw, nB);
chirpa = nan(nw, nB);
Sigma0 = nan(1, nB);
phase = zeros(1, nB);
phaseRpa = zeros(1, nB);
rpaMass = nan(1, nB);
moment = nan(1, nB);
Dord = nan(1, nB);

nWorkers = 0;
if parallel && ~isempty(ver('parallel')), nWorkers = Inf; end
parfor (k = 1:nB, nWorkers)
    [chiz(:,k), chirpa(:,k), Sigma0(k), phase(k), phaseRpa(k), ...
        rpaMass(k), moment(k), Dord(k)] = ...
        one_field(ion, T, fields(k), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts);
    if verbose
        labels = {'masked', 'ordered', 'paramagnet'};
        fprintf('  B = %5.2f T : 1/z %-11s RPA %-11s Sigma0 = %s\n', ...
            fields(k), labels{phase(k)+1}, labels{phaseRpa(k)+1}, num2str(Sigma0(k)));
    end
end

S = struct('fields', fields, 'w', w, 'T', T, 'info', info, ...
    'chiz', chiz, 'chirpa', chirpa, 'Sigma0', Sigma0, ...
    'phase', phase, 'phase_rpa', phaseRpa, 'rpa_mass_pm', rpaMass, ...
    'moment', moment, 'D_ord', Dord);
S.hmf_integral_mode = getf(sxtra,'hmf_integral_mode','full_profile');
S.visual_only = ismember(S.hmf_integral_mode, ...
    {'endpoint_trapezoid_visual','filtered_profile_visual'});
S.Bc = invz_boundary_field(fields, phase == 1, phase == 2);
S.Bc_1z = S.Bc;
S.Bc_rpa = invz_boundary_field(fields, phaseRpa == 1, phaseRpa == 2);
S.Epeak = invz_peak_energy(chiz, w, wmin);
S.Epeak_rpa = invz_peak_energy(chirpa, w, wmin);
end

function [chiz, chirpa, Sigma0, phase, phaseRpa, rpaMass, moment, Dord] = ...
    one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp, sopts)
nw = numel(w);
chiz = nan(nw,1);
chirpa = nan(nw,1);
Sigma0 = NaN;
rpaMass = NaN;
moment = NaN;
Dord = NaN;

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);
if phase ~= 0 && ~isempty(pt)
    copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, ...
        'Jshape', Jshape, 'hyp', hyp, 'si', pt.si);
    try
        out = invz_chi_realaxis(ion, T, B, pt, w, copts);
        chiz = imag(out.chi_cc_q(1,:)).';
        Sigma0 = pt.Sigma0;
        if phase == 1
            moment = pt.m0;
            Dord = pt.D_uni;
        end
    catch err
        if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        phase = 0;
        Sigma0 = NaN;
        moment = NaN;
        Dord = NaN;
    end
end

% Independent bare-RPA leg. J0eff selects/constructs the MF state; the exact same
% scalar is Jsel in the denominator. This remains available when the 1/z column is masked.
ropts = struct('Jsel', Jcc0, 'J0eff', Jcc0, 'bare_rpa', true, ...
    'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp, 'npass', 1);
try
    out0 = invz_chi_realaxis(ion, T, B, [], w, ropts);
    chirpa = imag(out0.chi_cc_q(1,:)).';
    phaseRpa = out0.phase_rpa;
    rpaMass = out0.rpa_mass_pm;
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    phaseRpa = 0;
end
end
