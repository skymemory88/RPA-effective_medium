function loc = invzf_electronuclear_local(ion, T, B, h, wn, opts)
%INVZF_ELECTRONUCLEAR_LOCAL Source-biased full local oracle for the strict pilot.
%
%   LOC = INVZF_ELECTRONUCLEAR_LOCAL(ION,T,B,H,WN,OPTS) diagonalizes the
%   full electronic or electronuclear single-ion Hamiltonian with the total
%   projected longitudinal source -H*Jz and returns a positive connected
%   C2(iwn) from the same Hamiltonian.  B must be purely transverse; an
%   external Bz belongs in the common functional as Hext=gL*muB*Bz, while H
%   here is the total local source varied by that functional.
%
%   To keep one exact local generator, this first oracle requires
%   transverse_mf='none'.  Reusing the production transverse mean field
%   would require adding its transverse moment and conjugate field as
%   independent functional variables; a partial feedback correction is not
%   admitted here.
%
%   opts.hyp (true) selects the full 136-state electronuclear manifold.
%   opts.h_step overrides the source finite-difference scale.  Seven source
%   nodes supply nested five-point stencils and sixth-order Richardson
%   extrapolation for dC2/dh and d2C2/dh2.  opts.beta_rel_step (5e-4)
%   controls an analogous spectrum-only beta derivative.

if nargin < 6 || isempty(opts), opts = struct(); end
validateattributes(T, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(h, {'numeric'}, {'real','scalar','finite'});
validateattributes(wn, {'numeric'}, {'real','vector','finite','integer'});
B = invz_field_vec(B);
if abs(B(3)) > 64*eps(max(1,max(abs(B))))
    error('invzf:electronuclearBz', ...
        ['B must be purely transverse. Put gL*muB*Bz in the common-functional ' ...
         'external source and pass the total local h separately.']);
end
tmf = get_opt(opts,'transverse_mf','none');
if ~strcmp(tmf,'none')
    error('invzf:transverseMFUnclosed', ...
        ['The scalar electronuclear oracle requires transverse_mf=''none''. ' ...
         'Other modes need additional transverse functional variables.']);
end
hyp = get_opt(opts,'hyp',true);
valid_hyp = isscalar(hyp) && (islogical(hyp) || ...
    (isnumeric(hyp) && isreal(hyp) && isfinite(hyp) && any(hyp == [0 1])));
if ~valid_hyp
    error('invzf:electronuclearOpts', ...
        'opts.hyp must be scalar logical or numeric 0/1.');
end
hyp = logical(hyp);

C = invz_const();
beta = 1/(C.kB*T);
if hyp
    if ~isfield(ion,'A')
        error('invzf:electronuclearIon', ...
            'ion.A is required when opts.hyp=true.');
    end
    validateattributes(ion.A, {'numeric'}, {'real','scalar','finite'});
    hyperfine_scale = abs(ion.A);
else
    hyperfine_scale = 0;
end
scale = max([abs(h),C.kB*T,hyperfine_scale,1e-3]);
dh = get_opt(opts,'h_step',eps('double')^(1/6)*scale);
validateattributes(dh, {'numeric'}, {'real','scalar','finite','positive'});
br = get_opt(opts,'beta_rel_step',5e-4);
validateattributes(br, {'numeric'}, {'real','scalar','finite','positive','<',0.2});

offset = [-2 -1 -0.5 0 0.5 1 2].';
states = cell(numel(offset),1);
for k = 1:numel(offset)
    states{k} = local_state(ion,T,B,h+dh*offset(k),wn,hyp,beta);
end
ic = find(offset == 0,1);
s0 = states{ic};

f = cellfun(@(s) s.f0,states);
m = cellfun(@(s) s.m,states);
C2 = cell2mat(cellfun(@(s) s.C2.',states,'UniformOutput',false));
[dfdh,dfdh_err] = richardson_first(f,dh);
[d2fdh2,d2fdh2_err] = richardson_second(f,dh);
[dmdh,dmdh_err] = richardson_first(m,dh);
[dC2dh,dC2dh_err] = richardson_first(C2,dh);
[d2C2dh2,d2C2dh2_err] = richardson_second(C2,dh);

db = br*beta;
cb = zeros(7,numel(wn));
for k = 1:numel(offset)
    bk = beta+db*offset(k);
    if bk <= 0
        error('invzf:electronuclearBetaStep', 'beta stencil crossed zero.');
    end
    cb(k,:) = c2_spectrum(s0.E,s0.O,bk,wn).';
end
[dC2dbeta,dC2dbeta_err] = richardson_first(cb,db);

chi = c2_spectrum(s0.E,s0.O,beta,0);
finite_values = [s0.f0;s0.u0;s0.m;chi;s0.C2(:);dC2dh(:); ...
    d2C2dh2(:);dC2dbeta(:);s0.tail_A];
status = 'ok';
if any(~isfinite(finite_values)), status = 'nonfinite'; end

loc = struct('status',status,'T',T,'B',B,'h',h,'beta',beta, ...
    'wn',wn(:),'omega',2*pi*wn(:)/beta,'f0',s0.f0,'u0',s0.u0, ...
    'm',s0.m,'chi',real(chi),'C2',real(s0.C2), ...
    'dC2dh',real(dC2dh(:)),'d2C2dh2',real(d2C2dh2(:)), ...
    'dC2dbeta',real(dC2dbeta(:)),'tail_A',s0.tail_A, ...
    'single_ion',s0.si,'source_step',dh,'beta_step',db, ...
    'derivative_error',struct( ...
        'dfdh',dfdh_err,'d2fdh2',d2fdh2_err,'dmdh',dmdh_err, ...
        'dC2dh',dC2dh_err(:),'d2C2dh2',d2C2dh2_err(:), ...
        'dC2dbeta',dC2dbeta_err(:)), ...
    'identity_error',struct( ...
        'm_plus_dfdh',s0.m+dfdh, ...
        'chi_plus_d2fdh2',chi+d2fdh2, ...
        'chi_minus_dmdh',chi-dmdh), ...
    'provenance',struct('schema','invzf_electronuclear_local/v1', ...
        'hyp',hyp,'transverse_mf','none','production_dispatch',false));
end

function s = local_state(ion,T,B,h,wn,hyp,beta)
si = invz_single_ion(ion,T,B, ...
    struct('hyp',hyp,'hz_fixed',h,'transverse_mf','none'));
logZrel = logsumexp(-beta*si.E);
f0 = si.E0-logZrel/beta;
u0 = sum(si.P.*(si.E0+si.E));
O = si.Mz;
C2 = c2_spectrum(si.E,O,beta,wn);
E = si.E(:);
p = si.P(:);
d = E.'-E;
dp = p-p.';
tail_S = real(sum(dp.*d.*abs(O).^2,'all'));
tail_A = max(0,tail_S)*beta^2/(4*pi^2);
s = struct('si',si,'E',E,'O',O,'f0',f0,'u0',u0, ...
    'm',si.Jexp(3),'C2',C2,'tail_A',tail_A);
end

function C2 = c2_spectrum(E,O,beta,wn)
E = E(:);
n = numel(E);
p = exp(-beta*(E-min(E)));
p = p/sum(p);
m = real(diag(O)).'*p;
d = E.'-E;
pa = repmat(p,1,n);
pb = pa.';
dp = pa-pb;
weight = abs(O).^2;
x = beta*d;
q0 = zeros(size(d));
small = abs(x) <= 1e-4;
q0(small) = beta*pa(small).*exprel_local(-x(small));
q0(~small) = dp(~small)./d(~small);

wn = wn(:);
C2 = zeros(size(wn));
for k = 1:numel(wn)
    if wn(k) == 0
        C2(k) = real(sum(q0.*weight,'all')-beta*m^2);
    else
        omega = 2*pi*wn(k)/beta;
        num = dp;
        num(small) = q0(small).*d(small);
        C2(k) = real(sum(num.*weight./(d-1i*omega),'all'));
    end
end
end

function y = exprel_local(x)
y = ones(size(x));
small = abs(x) < 1e-5;
xs = x(small);
y(small) = 1+xs.*(1/2+xs.*(1/6+xs.*(1/24+xs.*(1/120+xs/720))));
y(~small) = expm1(x(~small))./x(~small);
end

function z = logsumexp(x)
xmax = max(x);
z = xmax+log(sum(exp(x-xmax)));
end

function [d,err] = richardson_first(y,h)
dfull = (y(1,:)-8*y(2,:)+8*y(6,:)-y(7,:))/(12*h);
dhalf = (y(2,:)-8*y(3,:)+8*y(5,:)-y(6,:))/(6*h);
d = (16*dhalf-dfull)/15;
err = abs(dhalf-dfull)/15;
end

function [d2,err] = richardson_second(y,h)
d2full = (-y(7,:)+16*y(6,:)-30*y(4,:)+16*y(2,:)-y(1,:))/(12*h^2);
d2half = 4*(-y(6,:)+16*y(5,:)-30*y(4,:)+16*y(3,:)-y(2,:))/(12*h^2);
d2 = (16*d2half-d2full)/15;
err = abs(d2half-d2full)/15;
end

function v = get_opt(opts,name,default)
if isfield(opts,name) && ~isempty(opts.(name)), v = opts.(name); else, v = default; end
end
