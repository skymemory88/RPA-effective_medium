function result = invzp_outer_zero_census()
%INVZP_OUTER_ZERO_CENSUS Map the admissible outer domain at Sigma=0.
% Diagnostic only: a defined/undefined zero-state map is not a coupled-root
% existence verdict.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

F = load(fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat'));
D = load(fullfile(repo,'diag_rev3_check.mat'));
ion = invz_ion();
T = D.T;
[wn,wts,beta] = invz_matsubara(T,40);
opts = struct('emt_static',struct('Jsup',F.info.Jcc0,'warn',false));

nrow = sum(arrayfun(@(q)numel(q.nodes),D.out));
Bx_col = nan(nrow,1);
node_col = nan(nrow,1);
h_col = nan(nrow,1);
status_col = strings(nrow,1);
defined_col = false(nrow,1);
nroot_col = nan(nrow,1);
resid_col = nan(nrow,1);
Duni_col = nan(nrow,1);
mesh_margin_col = nan(nrow,1);
row = 0;
for ib = 1:numel(D.out)
    Bx = D.out(ib).Bx;
    for inode = 1:numel(D.out(ib).nodes)
        row = row+1;
        h = D.out(ib).nodes{inode}.h;
        si = invz_single_ion(ion,T,[Bx 0 0], ...
            struct('hyp',true,'Jxx0',F.info.Jaa0,'hz_fixed',h));
        tl = invz_twolevel_ordered(ion,T,Bx,h,struct('Jxx0',F.info.Jaa0));
        ctx = make_context(si,tl,T,wn,wts,beta,F.J,F.info);
        oo = invz_ordered_outer_map(zeros(size(wn)),ctx,opts);
        Bx_col(row) = Bx;
        node_col(row) = inode-1;
        h_col(row) = h;
        status_col(row) = string(oo.status);
        defined_col(row) = oo.defined;
        nroot_col(row) = oo.static.n_admissible_roots;
        if oo.defined
            resid_col(row) = oo.residual_norm;
            Duni_col(row) = oo.static.D_uni;
            mesh_margin_col(row) = ...
                oo.static.root_table.min_mesh_medium_signed(oo.static.selected_index);
        end
    end
end

tab = table(Bx_col,node_col,h_col,status_col,defined_col,nroot_col, ...
    resid_col,Duni_col,mesh_margin_col, ...
    'VariableNames',{'Bx','node','h','status','defined','n_admissible_roots', ...
                     'residual_norm','D_uni','mesh_margin'});
[Bx_summary,~,group] = unique(tab.Bx);
n_nodes = splitapply(@numel,tab.defined,group);
n_defined = splitapply(@nnz,tab.defined,group);
summary = table(Bx_summary,n_nodes,n_defined, ...
    'VariableNames',{'Bx','n_nodes','n_defined'});
result = struct('table',tab,'summary',summary, ...
    'note',['Sigma=0 domain census only; undefined does not prove no coupled root, ' ...
            'and defined does not prove a coupled fixed point.'], ...
    'coupling_provenance',F.provenance.coupling_opts);
save(fullfile(here,'wp2_outer_zero_census.mat'),'result','-v7');
disp(summary);
end

function ctx = make_context(si,tl,T,wn,wts,beta,J,info)
c0 = invz_chi0z(si,T,1i*wn,struct('elastic',true));
G0 = -real(squeeze(c0(3,3,:)));
c0i = invz_chi0z(si,T,1i*wn(1),struct('elastic',false));
G0i = -real(c0i(3,3,1));
X = real(c0(:,:,1));
feedback = X(3,1)*(info.Jaa0/(1-info.Jaa0*X(1,1)))*X(1,3);
G0e = -(X(3,3)+feedback)-G0i;
g = real(invz_g(tl,1i*wn));
ctx = struct('G0',G0,'g',g,'tl',tl,'wts',wts,'beta',beta, ...
    'Jnu_flat',J,'J0eff',info.Jcc0,'G0inel0',G0i,'G0el0',G0e);
end
