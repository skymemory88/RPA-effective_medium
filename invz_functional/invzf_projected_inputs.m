function pilot = invzf_projected_inputs(ion, T, B, Jnu, info, Jaa0, opts)
%INVZF_PROJECTED_INPUTS Map production scalar inputs into the strict pilot.
%
%   PILOT = INVZF_PROJECTED_INPUTS(ION,T,B,JNU,INFO,JAA0,OPTS) uses the
%   production transverse-doublet constructor and coupling spectrum, but
%   returns only inputs for the disconnected invz_functional prototype.
%
%   The fixed transverse doublet is subsequently diagonalized exactly under
%   the projected longitudinal source.  Thus Bz enters as
%       H = ion.gL*muB*Bz
%   while the doublet itself is constructed at [Bx By 0].  Hyperfine,
%   ODD/retarded couplings, and the production ordered tanh/xi replacement
%   are not part of this scalar pilot.
%
%   opts.transverse_mf (default 'legacy_x') is forwarded to invz_twolevel.
%   opts.mode_weights may supply a nonuniform normalized BZ measure.

if nargin < 7 || isempty(opts), opts = struct(); end
validateattributes(T, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(Jnu, {'numeric'}, {'real','vector','finite','nonempty'});
validateattributes(Jaa0, {'numeric'}, {'real','scalar','finite'});
if ~(isstruct(info) && isscalar(info) && isfield(info,'Jcc0'))
    error('invzf:projectedInfo', 'info.Jcc0 is required from invz_bz_couplings.');
end
validateattributes(info.Jcc0, {'numeric'}, {'real','scalar','finite'});

B = invz_field_vec(B);
tmf = get_opt(opts,'transverse_mf','legacy_x');
tlo = struct('Jxx0',Jaa0,'transverse_mf',tmf);
tl = invz_twolevel(ion,T,[B(1:2) 0],tlo);
C = invz_const();

model = struct('Delta',tl.Delta,'M',sqrt(tl.M2), ...
    'beta',1/(C.kB*T));
lattice = struct('J0',info.Jcc0,'Jmodes',Jnu(:));
if isfield(opts,'mode_weights') && ~isempty(opts.mode_weights)
    validateattributes(opts.mode_weights, {'numeric'}, ...
        {'real','vector','finite','nonnegative','numel',numel(Jnu)});
    lattice.mode_weights = opts.mode_weights(:);
end

pilot = struct();
pilot.model = model;
pilot.lattice = lattice;
pilot.H = ion.gL*C.muB*B(3);
pilot.local_twolevel = tl;
pilot.T = T;
pilot.B = B;
pilot.Jaa0 = Jaa0;
pilot.info = info;
pilot.provenance = struct( ...
    'schema','invzf_projected_inputs/v1', ...
    'local_model','fixed_transverse_electronic_doublet', ...
    'longitudinal_source','ion.gL*muB*Bz', ...
    'hyperfine',false, ...
    'production_dispatch',false);
end

function v = get_opt(opts,name,default)
if isfield(opts,name) && ~isempty(opts.(name)), v = opts.(name); else, v = default; end
end
