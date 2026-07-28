function out = invzf_centered_cluster_exact(Delta, M, h, J, beta, wn, sites)
%INVZF_CENTERED_CLUSTER_EXACT Exact cluster with fixed local centring.
%
%   OUT = INVZF_CENTERED_CLUSTER_EXACT(DELTA,M,H,J,BETA,WN,SITES)
%   diagonalizes
%
%     Hc = sum_i [DELTA*|1><1|_i-H*X_i]
%          -sum_{i<j} J_ij*(X_i-m0)*(X_j-m0),
%
%   where m0 is the exact uncoupled one-site moment at source H and is held
%   fixed for the whole coupling scan.  The implementation maps this to
%   site sources h_i=H-m0*sum_j J_ij plus the exact constant
%   -m0^2*sum_{i<j}J_ij.

if nargin < 7, sites = []; end
validateattributes(J, {'numeric'}, {'real','finite','2d','square'});
Jc = (J+J.')/2;
Jc(1:size(Jc,1)+1:end) = 0;
loc0 = invzf_twolevel_local(Delta,M,h,beta,0);
m0 = loc0.m;
hsite = h-m0*sum(Jc,2);
cluster = invzf_cluster_exact(Delta,M,hsite,J,beta,wn,sites);
constant = -0.5*m0^2*sum(cluster.J,'all');

out = struct('F',cluster.F+constant,'U',cluster.U+constant, ...
    'deltaF',cluster.F+constant-size(J,1)*loc0.f0, ...
    'm0',m0,'h',h,'h_site',hsite,'constant',constant, ...
    'local_reference',loc0,'cluster',cluster, ...
    'provenance',struct('schema','invzf_centered_cluster_exact/v1', ...
        'reference_moment_fixed',true,'physical_continuation',false));
end
