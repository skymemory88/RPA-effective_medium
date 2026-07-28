function ctx = invz_ordered_node_context(ion, T, Bx, Jnu_flat, opts)
%INVZ_ORDERED_NODE_CONTEXT Freeze h-independent ordered-node inputs.
%
%   Diagnostic-only companion to INVZ_ORDERED_NODE_NEWTON. The context
%   mirrors INVZ_HMF_ORDERED's fixed-h node construction for the unchanged
%   resummed equations, but keeps full states available to continuation.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isstruct(opts) || ~isscalar(opts)
    error('invzp:OrderedContext:InvalidInput','opts must be a scalar struct.');
end
if ~isfield(opts,'J0eff')
    error('invzp:OrderedContext:InvalidInput','opts.J0eff is required.');
end
J0eff = opts.J0eff;
if ~isnumeric(J0eff) || ~isreal(J0eff) || ~isscalar(J0eff) || ~isfinite(J0eff)
    error('invzp:OrderedContext:InvalidInput', ...
        'opts.J0eff must be a finite real scalar.');
end
if ~isnumeric(Jnu_flat) || ~isreal(Jnu_flat) || ~isvector(Jnu_flat) || ...
        isempty(Jnu_flat) || any(~isfinite(Jnu_flat),'all')
    error('invzp:OrderedContext:InvalidInput', ...
        'Jnu_flat must be a nonempty finite real vector.');
end
Jnu_flat = Jnu_flat(:);

Jxx0 = getf(opts,'Jxx0',ion.Jxx0);
hyp = getf(opts,'hyp',true);
tmf = getf(opts,'transverse_mf','legacy_x');
Ecut = getf(opts,'Ecut',40);
eopts = getf(opts,'emt',struct());
eso = getf(opts,'emt_static',struct());
[sm,eopts,eso] = invz_check_static_medium(opts,eopts,eso);
if ~strcmp(sm.scheme,'resummed')
    error('invzp:OrderedContext:Scheme', ...
        'Continuation context currently supports static_medium=''resummed'' only.');
end
eso.warn = false;
[wn,wts,beta] = invz_matsubara(T,Ecut);
Jmom = getf(opts,'Jmom',[]);
if isempty(Jmom), Jmom = invz_coupling_moments(Jnu_flat); end

sibase = struct('hyp',hyp,'Jxx0',Jxx0,'transverse_mf',tmf);
for name = {'mf_maxit','mf_mix'}
    if isfield(opts,name{1}), sibase.(name{1}) = opts.(name{1}); end
end

ctx = struct( ...
    'schema','invzp_ordered_node_context/v1', ...
    'ion',ion,'T',T,'Bx',Bx,'Jnu_flat',Jnu_flat, ...
    'J0eff',J0eff,'Jxx0',Jxx0,'transverse_mf',tmf, ...
    'wn',wn,'wts',wts,'beta',beta,'Ecut',Ecut, ...
    'eopts',eopts,'eso',eso,'Jmom',Jmom,'sibase',sibase, ...
    'Jscale',max(abs(Jnu_flat)));
if ~(isfinite(ctx.Jscale) && ctx.Jscale > 0)
    error('invzp:OrderedContext:InvalidInput', ...
        'The coupling vector must have a positive finite scale.');
end
end

function value = getf(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end
