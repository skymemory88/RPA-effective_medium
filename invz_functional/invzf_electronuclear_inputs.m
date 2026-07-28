function pilot = invzf_electronuclear_inputs(ion, T, B, Jnu, info, opts)
%INVZF_ELECTRONUCLEAR_INPUTS Build the full-local strict pilot inputs.
%
%   PILOT = INVZF_ELECTRONUCLEAR_INPUTS(ION,T,B,JNU,INFO,OPTS) separates
%   the physical longitudinal field into
%       H = ion.gL*muB*Bz
%   and constructs model.local_function(h,wn) with h the total local source.
%   The local Hamiltonian therefore receives [Bx By 0] plus -h*Jz, avoiding
%   double counting of Bz.
%
%   opts.local is forwarded to invzf_electronuclear_local.  Its current
%   scalar contract requires transverse_mf='none'.  Production dispatch is
%   not installed.

if nargin < 6 || isempty(opts), opts = struct(); end
validateattributes(T, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(Jnu, {'numeric'}, {'real','vector','finite','nonempty'});
if ~(isstruct(info) && isscalar(info) && isfield(info,'Jcc0'))
    error('invzf:electronuclearInfo', ...
        'info.Jcc0 is required from invz_bz_couplings.');
end
validateattributes(info.Jcc0, {'numeric'}, {'real','scalar','finite'});
B = invz_field_vec(B);
C = invz_const();
lo = get_opt(opts,'local',struct());
if isfield(lo,'transverse_mf') && ~strcmp(lo.transverse_mf,'none')
    error('invzf:transverseMFUnclosed', ...
        'The scalar electronuclear pilot requires opts.local.transverse_mf=''none''.');
end
lo.transverse_mf = 'none';
Btrans = [B(1:2) 0];

pilot = struct();
pilot.model = struct('local_function', ...
    @(h,wn) invzf_electronuclear_local(ion,T,Btrans,h,wn,lo));
pilot.lattice = struct('J0',info.Jcc0,'Jmodes',Jnu(:));
if isfield(opts,'mode_weights') && ~isempty(opts.mode_weights)
    validateattributes(opts.mode_weights, {'numeric'}, ...
        {'real','vector','finite','nonnegative','numel',numel(Jnu)});
    pilot.lattice.mode_weights = opts.mode_weights(:);
end
pilot.H = ion.gL*C.muB*B(3);
pilot.T = T;
pilot.B = B;
pilot.Btrans = Btrans;
pilot.info = info;
pilot.provenance = struct('schema','invzf_electronuclear_inputs/v1', ...
    'local_model','full_source_biased_electronuclear', ...
    'longitudinal_source','ion.gL*muB*Bz', ...
    'transverse_mf','none','production_dispatch',false);
end

function v = get_opt(opts,name,default)
if isfield(opts,name) && ~isempty(opts.(name)), v = opts.(name); else, v = default; end
end
