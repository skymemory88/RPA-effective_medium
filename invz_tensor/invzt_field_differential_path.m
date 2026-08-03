function out = invzt_field_differential_path(Cbare, Ctilde, Jgamma, qpath, opts)
%INVZT_FIELD_DIFFERENTIAL_PATH Tensor modified-field differential on one path.
%   OUT = INVZT_FIELD_DIFFERENTIAL_PATH(CBARE,CTILDE,JGAMMA,QPATH,OPTS)
%   solves, on CTILDE's active subspace,
%       CTILDE * du0/dh = CBARE * duMF/dh.
%   QPATH=duMF/dh declares the constrained single-ion path. With fixed
%   applied transverse field, the transverse components must additionally
%   satisfy du0_t/dh = [JGAMMA*dm/dh]_t. OUT reports that compatibility
%   residual; it does not assume it vanishes.
if nargin<5,opts=struct();end
rank_tol=getf(opts,'rank_tol',1e-12);
transverse=getf(opts,'transverse_indices',[1 2]);
for x={Cbare,Ctilde,Jgamma}
    if ~(isnumeric(x{1})&&isequal(size(x{1}),[3 3])&&all(isfinite(x{1}(:))))
        error('invzt:fieldDifferentialMatrix','Response/coupling inputs must be finite 3x3 matrices.');
    end
end
if ~(isnumeric(qpath)&&isreal(qpath)&&numel(qpath)==3&&all(isfinite(qpath(:))))
    error('invzt:fieldDifferentialPath','qpath must contain three finite real entries.');
end
if ~(isnumeric(transverse)&&isvector(transverse) ...
        && all(ismember(transverse,[1 2])) ...
        && numel(unique(transverse))==numel(transverse))
    error('invzt:fieldDifferentialTransverse', ...
        'transverse_indices must be unique indices from [1 2].');
end
Cbare=real((Cbare+Cbare')/2);
Ctilde=real((Ctilde+Ctilde')/2);
Jgamma=real((Jgamma+Jgamma')/2);
qpath=qpath(:); transverse=transverse(:).';
[V,D]=eig(Ctilde); d=real(diag(D)); scale=max(1,max(abs(d)));
active=d>rank_tol*scale;
if ~any(active)
    error('invzt:fieldDifferentialRank','Ctilde has no positive active subspace.');
end
Cpinv=V(:,active)*diag(1./d(active))*V(:,active)';
projector=V(:,active)*V(:,active)';
A=Cpinv*Cbare;
dm=Cbare*qpath;
du0=A*qpath;
equation_residual=Ctilde*du0-dm;
range_residual=(eye(3)-projector)*Cbare;
interaction_slope=Jgamma*dm;
transverse_residual=du0(transverse)-interaction_slope(transverse);

% Locally consistent fixed-applied-transverse path. The modified transverse
% fields are not assumed to obey the bare MF constraint: solve for q_t such
% that the required du0_t equals Jgamma_t*dm while q_3=1.
R=A-Jgamma*Cbare;
Rt=R(transverse,transverse);
fixed=struct('defined',false,'rcond',NaN,'qpath',nan(3,1), ...
    'moment_slope',nan(3,1),'required_u0_slope',nan(3,1), ...
    'transverse_residual',nan(numel(transverse),1), ...
    'transverse_residual_inf',NaN,'longitudinal_integrand',NaN);
if ~isempty(transverse)
    fixed.rcond=rcond(Rt);
    if fixed.rcond>=rank_tol
        qfixed=zeros(3,1); qfixed(3)=1;
        qfixed(transverse)=-(Rt\R(transverse,3));
        dmfixed=Cbare*qfixed; du0fixed=A*qfixed;
        jdmfixed=Jgamma*dmfixed;
        tres=du0fixed(transverse)-jdmfixed(transverse);
        fixed=struct('defined',true,'rcond',rcond(Rt),'qpath',qfixed, ...
            'moment_slope',dmfixed,'required_u0_slope',du0fixed, ...
            'transverse_residual',tres, ...
            'transverse_residual_inf',max(abs(tres)), ...
            'longitudinal_integrand',du0fixed(3));
    end
end

out=struct('schema','invzt_field_differential_path/v1', ...
    'Cbare',Cbare,'Ctilde',Ctilde,'Jgamma',Jgamma,'qpath',qpath, ...
    'active_rank',nnz(active),'active_eigenvalues',d(active), ...
    'range_residual_fro',norm(range_residual,'fro'), ...
    'A_pseudoinverse',A,'moment_slope',dm,'required_u0_slope',du0, ...
    'equation_residual_inf',max(abs(equation_residual)), ...
    'interaction_u0_slope',interaction_slope, ...
    'transverse_indices',transverse, ...
    'transverse_compatibility_residual',transverse_residual, ...
    'transverse_compatibility_inf',max(abs(transverse_residual)), ...
    'longitudinal_integrand',du0(3), ...
    'fixed_applied_transverse_path',fixed);
end

function v=getf(s,name,default)
if isfield(s,name)&&~isempty(s.(name)),v=s.(name);else,v=default;end
end
