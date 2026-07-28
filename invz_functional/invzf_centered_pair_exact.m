function out = invzf_centered_pair_exact(Delta, M, h, j, beta, wn)
%INVZF_CENTERED_PAIR_EXACT Exact fixed-reference C3-C3 pair oracle.
%
%   OUT = INVZF_CENTERED_PAIR_EXACT(DELTA,M,H,J,BETA,WN) diagonalizes
%
%     Hc(j) = sum_i [DELTA*|1><1|_i - H*X_i]
%             - J*(X_1-m0)*(X_2-m0),
%
%   where X=M*sigma_x and m0 is the exact one-site moment at J=0.  The
%   reference m0 is held fixed as J is varied.  Consequently the linear
%   coupling coefficient vanishes, and [J^3]F contains only the bare
%   three-line C3-C3 vacuum topology.
%
%   The centred Hamiltonian is evaluated through the existing uncentred
%   pair oracle using h_pair=H-J*m0 and the additive constant -J*m0^2.
%   This is an exact-cluster fixture, not a physical source-continuation
%   prescription.

validateattributes(j, {'numeric'}, {'real','scalar','finite'});
loc0 = invzf_twolevel_local(Delta,M,h,beta,0);
m0 = loc0.m;
hpair = h-j*m0;
pair = invzf_two_site_exact(Delta,M,hpair,j,beta,wn);
F = pair.F-j*m0^2;

out = struct('F',F,'deltaF',F-2*loc0.f0,'m0',m0, ...
    'h',h,'j',j,'h_pair',hpair,'constant',-j*m0^2, ...
    'local_reference',loc0,'pair',pair, ...
    'provenance',struct('schema','invzf_centered_pair_exact/v1', ...
        'reference_moment_fixed',true,'physical_continuation',false));
end
