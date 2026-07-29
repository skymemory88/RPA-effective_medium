function contract = invzf_wp0_contract(opts)
%INVZF_WP0_CONTRACT WP0 "freeze the microscopic contract" (invzp_convg_fix.md Sec.9, L389-401).
% Freezes and RECORDS the microscopic anchors -- NO solving happens here:
%   (a) ion/Hamiltonian anchors (invz_ion, verbatim) + which hyperfine/exchange/dipolar/
%       demagnetization terms are active under the production configuration;
%   (b) unit convention (meV/K/T) + the Boltzmann/mu_B constants actually used, read from
%       invz_common/invz_const.m (never hardcoded here);
%   (c) the Matsubara convention actually implemented, read VERBATIM from invz_matsubara.m;
%   (d) the q-grid/Fourier convention and self-site handling, read VERBATIM from the
%       invz_bz_couplings.m / invz_jq_modes.m header comment blocks;
%   (e) a two-backend (bruteforce/ewald) coupling-set fingerprint on the SAME grid: the
%       exact-byte digest (invz_common/invz_exact_numeric_digest) of the flattened Jnu, plus
%       Jcc0/Jaa0, for EACH backend.
%
% opts.grid  (default [16 16 16] -- invz_bz_couplings' own production default)
% opts.dpRng (default 30         -- invz_jq_modes' own production default)
% opts.ewald (default alpha=0.3 Ang^-1, r_cut=16 Ang, g_cut=3 Ang^-1, boundary=
%   'conducting_k0_omitted'). This is NOT a repo-frozen value: invz_jq_modes.m and
%   invz_dipole_ewald.m headers both state this function "does not synthesize frozen
%   defaults" for {alpha,r_cut,g_cut}, and a repo-wide search found no concrete numeric
%   default anywhere else. The values above satisfy invz_dipole_ewald.m's documented
%   accuracy guard (alpha*r_cut >= 4.5, here 4.8; g_cut/(2*alpha) >= 4.5, here 5.0) with a
%   modest margin, chosen BEFORE inspecting any backend-comparison result (see
%   contract.ewald_params_note).
if nargin < 1 || isempty(opts), opts = struct(); end
grid  = lopt(opts, 'grid',  [16 16 16]);
dpRng = lopt(opts, 'dpRng', 30);
eopts = lopt(opts, 'ewald', struct('alpha',0.3,'r_cut',16,'g_cut',3,'boundary','conducting_k0_omitted'));

root = fileparts(fileparts(mfilename('fullpath')));
ion = invz_ion();
C   = invz_const();

contract = struct();
contract.schema  = 'invzf_wp0_contract/v1';
contract.created = datestr(now);

% ---- (a) ion / Hamiltonian anchors -------------------------------------------------------
contract.ion = ion;   % every ion field, verbatim (scalars + lattice/basis matrices)
contract.active_terms = struct( ...
    'hyperfine', local_term(ion.A ~= 0, sprintf(['ion.A = %.6g meV (nonzero); wired via ' ...
        'opts.hyp into the 136-state J=8 (2J+1=17) x I=3.5 (2I+1=8) electronuclear basis ' ...
        '(invz_common/invz_single_ion.m)'], ion.A)), ...
    'exchange', local_term(ion.J12 ~= 0, sprintf(['ion.J12 = %.6g meV nn exchange (nonzero), ' ...
        'sign(J12) = %d enters Jcc = -gfac*dip_cc + sign(J12)*|J12|*ex_cc ' ...
        '(invz_projected/invz_jq_modes.m header)'], ion.J12, sign(ion.J12))), ...
    'dipolar', local_term(true, sprintf(['C.gfac = %.6g meV*Ang^3 always enters Jcc/Jaa ' ...
        '(invz_projected/invz_jq_modes.m); production backend = bruteforce'], C.gfac)), ...
    'demagnetization', local_term(ion.demag ~= 0, sprintf(['ion.demag = %.6g (0 = off/' ...
        'intrinsic, the R2007 benchmark, and the production default); nonzero folds a shape ' ...
        'term into info.Jshape_cc/info.Jaa0 only (invz_projected/invz_jq_modes.m L115-125), ' ...
        'never into info.Jcc0/Jnu (demag-invariant criticality)'], ion.demag)), ...
    'odd', local_term(ion.odd ~= 0, sprintf(['ion.odd = %.6g (0 = off/published scalar-cc ' ...
        'model, the production default); drivers opt in via opts.odd, brute-force-only'], ...
        ion.odd)) );

% ---- (b) unit convention + constants -----------------------------------------------------
contract.units = struct('energy','meV', 'temperature','K', 'field','T', ...
    'beta_formula', 'beta = 1/(kB*T)  (invz_common/invz_matsubara.m)', ...
    'constants', C, 'constants_source', 'invz_common/invz_const.m (read, not hardcoded)');

% ---- (c) Matsubara convention, verbatim --------------------------------------------------
contract.matsubara = struct('source_file', 'invz_common/invz_matsubara.m', ...
    'source_verbatim', fileread(fullfile(root, 'invz_common', 'invz_matsubara.m')));

% ---- (d) q-grid / Fourier convention, verbatim headers -----------------------------------
contract.qgrid = struct( ...
    'grid_used', grid, 'dpRng_used', dpRng, ...
    'production_default_grid', [16 16 16], 'production_default_dpRng', 30, ...
    'gamma_policy_default_route', ['legacy/absent-grid-policy route (no gridConvention/' ...
        'gridOffset/gammaPolicy option present): exact rows with all(abs(q))<=1e-12 are ' ...
        'dropped unconditionally -- invz_bz_couplings.m: "qc = qc(any(abs(qc) > 1e-12, 2), :)"'], ...
    'grid_convention_default_route', ['qVec_generator(ion.a,''mode'',''grid'',''grid'',grid,' ...
        '''range'',[-0.5 0.5]) -- the pre-Task-4 legacy call (invz_bz_couplings.m L46-60), ' ...
        'bit-identical to the absent-grid-policy path'], ...
    'bz_couplings_header_verbatim', local_header(fullfile(root,'invz_projected','invz_bz_couplings.m')), ...
    'jq_modes_header_verbatim',     local_header(fullfile(root,'invz_projected','invz_jq_modes.m')), ...
    'gamma_formula_bruteforce', ['Jgamma_cc = Jpath_base_cc + lorz*ones(4);  Jpath_base_cc = ' ...
        '-gfac*dip_sphere_cc(0) + Jex_cc(0)   (invz_jq_modes.m header)'], ...
    'gamma_formula_ewald', ['Jgamma_cc = Jpath_base_cc + lorz*ones(4);  Jpath_base_cc = ' ...
        '-gfac*dip_reg_cc(0) + Jex_cc(0) - lorz*ones(4)  ("Ewald adds 0 at Gamma", ' ...
        'invz_jq_modes.m header)'] );

% ---- (e) coupling-set fingerprint, both backends, SAME grid ------------------------------
common = struct('grid', grid, 'dpRng', dpRng, 'cache', false);  % cache=false: exercise the
                                                                 % constructor, not a cache hit
bfOpts = common;  bfOpts.dipole = 'bruteforce';
ewOpts = common;  ewOpts.dipole = 'ewald';  ewOpts.ewald = eopts;

contract.ewald_params_note = struct('opts_used', eopts, 'is_repo_frozen_default', false, ...
    'rationale', sprintf(['no concrete {alpha,r_cut,g_cut} default is frozen anywhere in ' ...
        'this repo (invz_jq_modes.m/invz_dipole_ewald.m headers say so explicitly, and a ' ...
        'repo-wide grep for constructed opts.ewald structs found none); chosen to satisfy ' ...
        'the documented accuracy guard alpha*r_cut>=4.5 (used %.4g) and g_cut/(2*alpha)>=4.5 ' ...
        '(used %.4g), BEFORE running any backend comparison'], eopts.alpha*eopts.r_cut, ...
        eopts.g_cut/(2*eopts.alpha)));

contract.coupling_opts = struct('bruteforce', bfOpts, 'ewald', ewOpts);
contract.backends = struct( ...
    'bruteforce', local_build_backend(ion, bfOpts), ...
    'ewald',      local_build_backend(ion, ewOpts));
end

% ==============================================================================================
function s = local_term(active, evidence)
s = struct('active', logical(active), 'evidence', evidence);
end

function v = lopt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name)), v = opts.(name); else, v = default; end
end

function txt = local_header(fname)
% Verbatim leading contiguous comment block (the "header") right after the function line.
lines = cellstr(splitlines(string(fileread(fname))));
out = {};  started = false;
for i = 1:numel(lines)
    l = lines{i};
    if ~started
        if startsWith(strtrim(l), 'function'), started = true; end
        continue;
    end
    lt = strtrim(l);
    if isempty(lt) || startsWith(lt, '%')
        out{end+1} = l; %#ok<AGROW>
    else
        break;
    end
end
txt = strjoin(out, newline);
end

function rec = local_build_backend(ion, bzOpts)
% Calls invz_bz_couplings exactly ONCE for this backend. Any error is recorded honestly
% (rec.ok = false) -- never faked, never silently downgraded to a default.
rec = struct('opts', bzOpts, 'ok', false, 'error_id', '', 'error_msg', '');
try
    [Jnu, info, Jaa0, detail] = invz_bz_couplings(ion, bzOpts);
    Jflat = Jnu(:);
    rec.ok         = true;
    rec.digest     = invz_exact_numeric_digest(Jflat);
    rec.Jcc0       = info.Jcc0;
    rec.Jaa0       = Jaa0;
    rec.Jnu_flat   = Jflat;
    rec.info_summary = local_info_summary(info);
    rec.detail_summary = struct('qvec', detail.qvec, 'weights', detail.weights, ...
        'Jnu_unflat', detail.Jnu_unflat, 'Juni', detail.Juni, 'flattening', detail.flattening);
catch ME
    rec.error_id  = ME.identifier;
    rec.error_msg = ME.message;
end
end

function s = local_info_summary(info)
keep = {'Jcc0_dipole','Jaa0_dipole','Jcc0','Jaa0','Jshape_cc','dpRng','dipole', ...
        'Jpath_base_cc','Jgamma_cc'};
s = struct();
for i = 1:numel(keep)
    if isfield(info, keep{i}), s.(keep{i}) = info.(keep{i}); end
end
end
