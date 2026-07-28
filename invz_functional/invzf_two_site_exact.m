function out = invzf_two_site_exact(Delta, M, h, j, beta, wn)
%INVZF_TWO_SITE_EXACT Dense exact oracle for the scalar two-level pair.
%
%   H2 = Hloc(h) x I + I x Hloc(h) - j X x X,
%   Hloc = Delta*|1><1| - h*X,  X=M*sigma_x.
%
%   OUT contains pair free energy, internal energy, total moment, total static
%   susceptibility, and the local connected C11(iwn).  No finite differences
%   are used.

validateattributes(Delta, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(M, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(h, {'numeric'}, {'real','scalar','finite'});
validateattributes(j, {'numeric'}, {'real','scalar','finite'});
validateattributes(beta, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(wn, {'numeric'}, {'real','vector','finite','integer'});

cl = invzf_cluster_exact(Delta,M,h,[0 j;j 0],beta,wn,1);
out = struct('F',cl.F,'U',cl.U,'m_pair',cl.m_total, ...
    'chi_pair',cl.chi_uniform,'C11',cl.C_local(:,1),'wn',cl.wn, ...
    'E',cl.E,'p',cl.p);
end
