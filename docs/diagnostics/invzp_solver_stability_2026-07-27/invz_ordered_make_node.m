function [node, meta] = invz_ordered_make_node(ctx, h)
%INVZ_ORDERED_MAKE_NODE Reproduce one fixed-h ordered node outside HMF.
%
%   The formulas mirror INVZ_HMF_ORDERED's nested eval_node constructor.
%   No nonlinear solve is performed.

if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx,'schema') || ...
        ~strcmp(ctx.schema,'invzp_ordered_node_context/v1')
    error('invzp:OrderedNode:InvalidContext', ...
        'ctx must come from invz_ordered_node_context.');
end
if ~isnumeric(h) || ~isreal(h) || ~isscalar(h) || ~isfinite(h)
    error('invzp:OrderedNode:InvalidInput','h must be a finite real scalar.');
end

sio = ctx.sibase;
sio.hz_fixed = h;
si = invz_single_ion(ctx.ion,ctx.T,ctx.Bx,sio);
if abs(si.hz-h) > 1e-12
    error('invz:hzFixed','hz_fixed not held: si.hz = %.6g vs %.6g',si.hz,h);
end
tl = invz_twolevel_ordered(ctx.ion,ctx.T,ctx.Bx,h, ...
    struct('Jxx0',ctx.Jxx0,'transverse_mf',ctx.transverse_mf, ...
    'domain_policy','return'));
meta = struct('status','ok','h',h,'m',si.Jexp(3), ...
    'Delta',tl.Delta,'G0bare0',NaN);
if ~tl.valid
    node = [];
    meta.status = 'degenerate_doublet';
    return
end

c0 = invz_chi0z(si,ctx.T,1i*ctx.wn,struct('elastic',true));
G0 = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si,ctx.T,1i*ctx.wn(1),struct('elastic',false));
G0inel0 = -real(c0i(3,3,1));
X = real(c0(:,:,1));
switch ctx.transverse_mf
    case 'none'
        feedback = 0;
    case 'legacy_x'
        feedback = X(3,1)*(ctx.Jxx0/(1-ctx.Jxx0*X(1,1)))*X(1,3);
    case 'vector_ab'
        transverse = [1 2];
        feedback = X(3,transverse)*(ctx.Jxx0* ...
            ((eye(2)-ctx.Jxx0*X(transverse,transverse))\X(transverse,3)));
    otherwise
        error('invz:transverseMF', ...
            'unknown transverse_mf ''%s''',ctx.transverse_mf);
end
G0bare0 = -(X(3,3)+feedback);
G0el0 = G0bare0-G0inel0;
g = real(invz_g(tl,1i*ctx.wn));

node = struct('tl',tl,'G0',G0,'g',g,'wts',ctx.wts, ...
    'wn',ctx.wn,'beta',ctx.beta,'J0eff',ctx.J0eff, ...
    'G0inel0',G0inel0,'G0el0',G0el0,'G0bare0',G0bare0, ...
    'eso',ctx.eso,'eopts',ctx.eopts,'Jnu_flat',ctx.Jnu_flat, ...
    'Jmom',ctx.Jmom);
meta.G0bare0 = G0bare0;
end
