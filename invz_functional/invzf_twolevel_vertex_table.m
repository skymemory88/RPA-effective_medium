function tab = invzf_twolevel_vertex_table(Delta, M, h, beta, cutoff)
%INVZF_TWOLEVEL_VERTEX_TABLE Freeze dynamic local vertices on a signed grid.
%
%   TAB = INVZF_TWOLEVEL_VERTEX_TABLE(...,CUTOFF) precomputes C2(n),
%   gamma3(n,m,-n-m), and gamma4(n,-n,m,-m) for n,m=-CUTOFF:CUTOFF.
%   gamma3 entries whose third frequency lies outside the stored grid are
%   still exact; valid_convolution identifies entries retained by the
%   unpadded finite-grid Phi33 contraction.

validateattributes(cutoff, {'numeric'}, ...
    {'real','scalar','finite','integer','nonnegative'});
labels = (-cutoff:cutoff).';
nf = numel(labels);
loc = invzf_twolevel_local(Delta,M,h,beta,labels);
gamma3 = complex(zeros(nf));
gamma4 = complex(zeros(nf));
valid = false(nf);
for i = 1:nf
    n = labels(i);
    for j = 1:nf
        m = labels(j);
        v3 = invzf_twolevel_1pi_vertex(Delta,M,h,beta,[n,m,-n-m]);
        v4 = invzf_twolevel_1pi_vertex(Delta,M,h,beta,[n,-n,m,-m]);
        if ~strcmp(v3.status,'ok') || ~strcmp(v4.status,'ok')
            error('invzf:twolevelVertexTable', ...
                'Local 1PI construction failed at n=%d, m=%d.',n,m);
        end
        gamma3(i,j) = v3.gamma;
        gamma4(i,j) = v4.gamma;
        valid(i,j) = abs(n+m) <= cutoff;
    end
end
if max(abs(imag(gamma3)),[],'all') <= 1024*eps(max(1,max(abs(gamma3),[],'all')))
    gamma3 = real(gamma3);
end
if max(abs(imag(gamma4)),[],'all') <= 1024*eps(max(1,max(abs(gamma4),[],'all')))
    gamma4 = real(gamma4);
end
tab = struct('status','ok','schema','invzf_twolevel_vertex_table/v1', ...
    'Delta',Delta,'M',M,'h',h,'beta',beta,'cutoff',cutoff, ...
    'labels',labels,'C2',loc.C2,'gamma3',gamma3,'gamma4',gamma4, ...
    'valid_convolution',valid,'local',loc, ...
    'convolution_policy','retain n,m,n+m only when all lie on signed grid');
end
