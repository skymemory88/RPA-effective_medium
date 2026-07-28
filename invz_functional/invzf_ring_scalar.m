function ring = invzf_ring_scalar(loc, Jmodes, mode_weights)
%INVZF_RING_SCALAR Strict scalar ring functional and analytic derivatives.
%
%   RING = INVZF_RING_SCALAR(LOC, JMODES, MODE_WEIGHTS) evaluates, per site,
%
%     f_ring = (1/(2*beta)) sum_n w_n mean_q[
%                  log(1 - J_q*C_n) + J_q*C_n ].
%
%   LOC is returned by INVZF_TWOLEVEL_LOCAL.  MODE_WEIGHTS defaults to a
%   uniform normalized measure.  A nonnegative Matsubara list uses w0=1 and
%   wn=2; a signed list uses unit weights and must not contain duplicates.
%
%   The function fails closed with status='pole' when any 1-J_q*C_n <= 0.
%   dBetaf_dBeta is the explicit derivative of beta*f_ring at fixed source
%   and integer Matsubara labels.  For a complete 0:N or -N:N grid,
%   f_tail_bound rigorously bounds the omitted scalar ring free-energy tail.

if nargin < 3 || isempty(mode_weights)
    mode_weights = ones(numel(Jmodes), 1)/numel(Jmodes);
end
validateattributes(Jmodes, {'numeric'}, {'real','vector','finite','nonempty'});
validateattributes(mode_weights, {'numeric'}, ...
    {'real','vector','finite','nonnegative','numel',numel(Jmodes)});

J = Jmodes(:);
qw = mode_weights(:);
qsum = sum(qw);
if ~(isfinite(qsum) && qsum > 0)
    error('invzf:ringWeights', 'mode_weights must have a finite positive sum.');
end
qw = qw/qsum;

n = loc.wn(:);
if ~isfield(loc,'status') || ~strcmp(loc.status,'ok') ...
        || any(~isfinite([loc.C2(:);loc.dC2dh(:);loc.d2C2dh2(:);loc.dC2dbeta(:)]))
    fw = nan(size(n));
    ring = failed_ring('local_nonfinite', NaN, fw, qw);
    return
end
if numel(unique(n)) ~= numel(n)
    error('invzf:ringFrequencies', 'Matsubara indices must not contain duplicates.');
end
if all(n >= 0)
    fw = ones(size(n));
    fw(n > 0) = 2;
elseif any(n > 0) && any(n < 0)
    if ~all(ismember(-n, n))
        error('invzf:ringFrequencies', ...
            'A signed Matsubara list must contain every +/-n partner.');
    end
    check_even_frequency(loc,n);
    fw = ones(size(n));
else
    error('invzf:ringFrequencies', ...
        'Use a nonnegative list or a signed list containing +/-n partners.');
end

nf = numel(n);
term = nan(nf,1);
qret = nan(nf,1);
pret = nan(nf,1);
min_den = inf;
for k = 1:nf
    C = loc.C2(k);
    den = 1 - J*C;
    min_den = min(min_den, min(den));
    if any(~isfinite(den)) || any(den <= 0)
        ring = failed_ring('pole', min_den, fw, qw);
        return
    end
    x = J*C;
    term(k) = sum(qw.*log1m_plus_x(x));
    qret(k) = sum(qw.*(J.^2*C./den));
    pret(k) = sum(qw.*(J.^2./den.^2));
end

beta = loc.beta;
f = sum(fw.*term)/(2*beta);
dfdh = -sum(fw.*qret.*loc.dC2dh)/(2*beta);
d2fdh2 = -sum(fw.*(pret.*loc.dC2dh.^2 ...
                   + qret.*loc.d2C2dh2))/(2*beta);
dBetaf_dBeta = -sum(fw.*qret.*loc.dC2dbeta)/2;

mu2 = sum(qw.*J.^2);
f_J2 = -mu2*sum(fw.*loc.C2.^2)/(4*beta);
dfdh_J2 = -mu2*sum(fw.*loc.C2.*loc.dC2dh)/(2*beta);
d2fdh2_J2 = -mu2*sum(fw.*(loc.dC2dh.^2 ...
                          + loc.C2.*loc.d2C2dh2))/(2*beta);

finite_out = [term;qret;pret;f;dfdh;d2fdh2;dBetaf_dBeta; ...
    mu2;f_J2;dfdh_J2;d2fdh2_J2];
if any(~isfinite(finite_out))
    ring = failed_ring('nonfinite', min_den, fw, qw);
    return
end

[f_tail_bound,tail_rho,tail_cutoff] = scalar_tail_bound(loc,J,qw,n);

ring = struct('status', 'ok', 'f', f, 'dfdh', dfdh, ...
    'd2fdh2', d2fdh2, 'fixed_reference_chi', -d2fdh2, ...
    'dBetaf_dBeta', dBetaf_dBeta, ...
    'term', term, 'q_return', qret, 'p_return', pret, ...
    'freq_weights', fw, 'mode_weights', qw, 'min_den', min_den, ...
    'mu2', mu2, 'f_J2', f_J2, 'dfdh_J2', dfdh_J2, ...
    'd2fdh2_J2', d2fdh2_J2, 'f_tail_bound', f_tail_bound, ...
    'tail_rho', tail_rho, 'tail_cutoff', tail_cutoff);
end

function y = log1m_plus_x(x)
% Stable log(1-x)+x.  The direct expression cancels its leading -x^2/2.
y = log1p(-x)+x;
small = abs(x) <= 1e-3;
xs = x(small);
y(small) = -xs.^2.*(1/2+xs.*(1/3+xs.*(1/4+xs.*(1/5 ...
    +xs.*(1/6+xs.*(1/7+xs/8))))));
end

function check_even_frequency(loc,n)
fields = {'C2','dC2dh','d2C2dh2','dC2dbeta'};
tol = 128*eps;
for k = find(n > 0).'
    partner = find(n == -n(k),1);
    for f = 1:numel(fields)
        a = loc.(fields{f})(k);
        b = loc.(fields{f})(partner);
        if abs(a-b) > tol*max([1,abs(a),abs(b)])
            error('invzf:ringFrequencySymmetry', ...
                'Signed local data violate %s(-n)=%s(n) at |n|=%d.', ...
                fields{f},fields{f},n(k));
        end
    end
end
end

function [bound,rho,N] = scalar_tail_bound(loc,J,qw,n)
% Rigorous positive-frequency remainder for the scalar ring value.  For n>0,
% C_n=A/(n^2+a^2)<=A/n^2 and
% |log(1-x)+x|<=x^2/[2(1-|x|)].
bound = NaN; rho = NaN; N = max(abs(n));
if N < 1, return; end
if all(n >= 0)
    complete = isequal(n,(0:N).');
else
    complete = isequal(sort(n),(-N:N).');
end
if ~complete, return; end
delta = loc.Delta/2;
b = delta^2/loc.E^2;
A = loc.M^2*b*4*loc.E*tanh(loc.beta*loc.E)*loc.beta^2/(4*pi^2);
rho = max(abs(J))*A/(N+1)^2;
if ~(isfinite(rho) && rho < 1), return; end
mu2 = sum(qw.*J.^2);
bound = mu2*A^2/(6*loc.beta*(1-rho)*N^3);
end

function ring = failed_ring(status, min_den, fw, qw)
ring = struct('status', status, 'f', NaN, 'dfdh', NaN, ...
    'd2fdh2', NaN, 'fixed_reference_chi', NaN, ...
    'dBetaf_dBeta', NaN, ...
    'term', [], 'q_return', [], 'p_return', [], ...
    'freq_weights', fw, 'mode_weights', qw, 'min_den', min_den, ...
    'mu2', NaN, 'f_J2', NaN, 'dfdh_J2', NaN, 'd2fdh2_J2', NaN, ...
    'f_tail_bound', NaN, 'tail_rho', NaN, 'tail_cutoff', NaN);
end
