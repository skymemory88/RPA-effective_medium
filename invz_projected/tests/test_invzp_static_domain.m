function results = test_invzp_static_domain()
%TEST_INVZP_STATIC_DOMAIN Focused acceptance gates for bounded ordered statics.
here = fileparts(mfilename('fullpath'));
projected = fileparts(here);
repo = fileparts(projected);
addpath(projected, fullfile(repo,'invz_common'), repo);

names = strings(0,1); metrics = zeros(0,1); tolerances = zeros(0,1);
notes = strings(0,1);

% Generic enumeration: multiple roots, touching root, and a discontinuity.
grid = linspace(1e-6, 1-1e-6, 401);
[rm, im] = invz_bounded_roots(@(x)(x-.2).*(x-.7), grid, ...
    struct('resid_tol',1e-12));
add_check("enumerate_two_roots", max(abs(rm-[.2;.7])), 1e-8, ...
    "both isolated roots found");
[rt, it] = invz_bounded_roots(@(x)(x-.4).^2, grid, ...
    struct('resid_tol',1e-12));
add_check("touching_root", abs(rt-.4), 1e-8, ...
    "even-multiplicity root found without a sign change");
grid_disc = linspace(1e-4, .9998, 400);
pole = sqrt(2)/3;
[rd, idisc] = invz_bounded_roots(@(x)1./(x-pole), grid_disc, ...
    struct('resid_tol',1e-12));
disc_ok = isempty(rd) && size(idisc.discontinuity_brackets,1) == 1;
add_check("reject_discontinuity", double(~disc_ok), 0, ...
    "pole sign change is not returned as a root");

% A sign bracket must be polished to x_tol, not stopped at the first point
% below resid_tol; downstream equivalent closures can have different local
% conditioning at an only-marginally polished point.
root_exact = sqrt(3)/5;
[rp, ip] = invz_bounded_roots(@(x)1e3*(x-root_exact), ...
    linspace(0,1,100),struct('resid_tol',1e-6,'x_tol',1e-14));
polish_err = max(abs([rp-root_exact;ip.root_residual/1e3]));
add_check("sign_root_polished_to_xtol", polish_err, 1e-12, ...
    "sign root is refined beyond its acceptance residual");

% m -> 0 identity against the ordinary paramagnetic scalar medium.
tl0 = struct('m',0,'M2',1,'n01',1,'g0',1);
Jtoy = [-1;-.25;.5]; Sigma0 = .5; G0i = -.2;
to = struct('Jsup',1,'warn',false,'resid_tol',1e-10,'scan_points',1025);
[Kpm, Gpm, opm] = invz_emt_static_ordered( ...
    tl0, [0;0], Sigma0, Jtoy, -1e6, 10, 1, G0i, 0, to);
med = invz_emt_scalar(G0i, Sigma0, Jtoy, struct('debug',true));
pm_err = max([abs(Kpm-med.K), abs(Gpm-med.G), ...
              abs(opm.Gtil0-G0i/(1+Sigma0))]);
add_check("m0_paramagnetic_identity", pm_err, 1e-10, ...
    "bounded ordered static map equals invz_emt_scalar");
[Kbad, Gbad, obad] = invz_emt_static_ordered( ...
    tl0, [0;-1], Sigma0, Jtoy, 0, 10, 1, G0i, 0, to);
bad_xi_ok = isnan(Kbad) && isnan(Gbad) && obad.n_roots == 1 && ...
    obad.n_admissible_roots == 0 && ...
    contains(obad.root_table.reject_reason(1),"nonphysical_xi");
add_check("reject_nonphysical_xi", double(~bad_xi_ok), 0, ...
    "finite mathematical root is recorded but not exported");
scalar_bitwise = legacy_scalar_bitwise_check();
add_check("gstat_scalar_bitwise", double(~scalar_bitwise), 0, ...
    "array-capable production closure preserves legacy scalar arithmetic");

% Frozen production lattice and healthy 4 T legacy fixture.
fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
F = load(fixture_path);
[~,~,beta] = invz_matsubara(F.provenance.T, F.provenance.solve_opts.Ecut);
C = invz_const();
ion = invz_ion();
Greg = F.info.Jpath_base_cc;
z2 = linspace(0,1,10001);
directional_top = nan(size(z2));
for k = 1:numel(z2)
    Jm = Greg + C.gfac*(4*pi/ion.Vc)*(1/3-z2(k))*ones(4);
    Jm = (Jm+Jm')/2;
    directional_top(k) = max(real(eig(Jm)));
end
[directional_max, imax] = max(directional_top);
sup_err = max([max(F.J)-F.info.Jcc0, abs(directional_max-F.info.Jcc0)]);
add_check("production_Jsup", sup_err, 5e-13, ...
    sprintf("directional maximum at kz^2=%.6g",z2(imax)));

[G0i4, G0e4] = static_split(F.pt.si, F.provenance.T, F.info.Jaa0);
so = struct('Jsup',F.info.Jcc0,'warn',false,'scan_points',4097, ...
    'resid_tol',1e-10);
[K4, G4, o4] = invz_emt_static_ordered(F.pt.tl, F.pt.lambda(1:2), ...
    F.pt.Sigma0, F.J, -1e9, beta, F.info.Jcc0, G0i4, G0e4, so);
healthy_err = max(abs([K4-F.pt.K(1), G4-F.pt.G(1)]));
add_check("healthy_4T_reproduction", healthy_err, 1e-9, ...
    "one admissible root, compared with pre-change fixture");
if ~(o4.converged && o4.n_admissible_roots == 1)
    error('invz:testStaticDomain', 'healthy 4 T fixture was not uniquely admissible.');
end

% Seed independence: the public compatibility seed cannot alter any result.
[K4b, G4b, o4b] = invz_emt_static_ordered(F.pt.tl, F.pt.lambda(1:2), ...
    F.pt.Sigma0, F.J, 1e9, beta, F.info.Jcc0, G0i4, G0e4, so);
seed_err = max(abs([K4-K4b, G4-G4b, o4.D_uni-o4b.D_uni]));
add_check("seed_independence", seed_err, 0, ...
    "opposite extreme K0 seeds give identical exported root");

% Frozen node 9 and 3 T nodes 21--27: unchanged pseudo-roots are rejected,
% and the complete bounded scan reports whether an alternative exists.
D = load(fullfile(repo,'diag_rev3_check.mat'));
node9 = D.out(1).nodes{10}; % stored index is zero-based in the diagnostic
[~,~,on9] = invz_emt_static_ordered(node9.tl, node9.lam(1:2), ...
    node9.Sigma0, F.J, node9.K0, beta, F.info.Jcc0, ...
    node9.G0inel0, node9.G0el0, so);
node9_ok = strcmp(on9.status,'no_admissible_static_root') && on9.n_roots == 0;
add_check("stored_1T_node9_rejected", double(~node9_ok), 0, ...
    "no alternative admissible frozen-outer static root found");

threeT_status = strings(7,1);
threeT_nroots = nan(7,1);
for inode = 21:27
    nd = D.out(2).nodes{inode+1};
    [~,~,oo] = invz_emt_static_ordered(nd.tl, nd.lam(1:2), nd.Sigma0, ...
        F.J, nd.K0, beta, F.info.Jcc0, nd.G0inel0, nd.G0el0, so);
    threeT_status(inode-20) = string(oo.status);
    threeT_nroots(inode-20) = oo.n_roots;
end
threeT_ok = all(threeT_status == "no_admissible_static_root") && ...
            all(threeT_nroots == 0);
add_check("stored_3T_nodes21_27_rejected", double(~threeT_ok), 0, ...
    "all seven frozen outer states have no admissible static root");

% Every accepted production/toy root has positive mesh and uniform margins.
accepted_tables = [opm.root_table; o4.root_table];
margin_min = min([accepted_tables.supremum_mass; accepted_tables.D_uni; ...
    accepted_tables.min_mesh_x_signed; accepted_tables.min_mesh_medium_signed]);
add_check("accepted_positive_margins", -margin_min, 0, ...
    sprintf("minimum signed margin %.9g",margin_min));

passed = metrics <= tolerances;
summary = table(names, passed, metrics, tolerances, notes);
if ~all(passed)
    disp(summary(~passed,:));
    error('invz:testStaticDomain', '%d bounded-static acceptance checks failed.', nnz(~passed));
end

results = struct();
results.summary = summary;
results.lattice = struct('J0eff',F.info.Jcc0,'mesh_max',max(F.J), ...
    'directional_max',directional_max,'directional_arg_kz2',z2(imax), ...
    'coupling_opts',F.provenance.coupling_opts);
results.healthy_4T = struct('K0',K4,'Gstat',G4,'out',o4);
results.pm_limit = struct('K0',Kpm,'Gstat',Gpm,'out',opm);
results.node9 = struct('status',on9.status,'n_roots',on9.n_roots, ...
    'n_admissible_roots',on9.n_admissible_roots);
results.threeT = table((21:27).',threeT_status,threeT_nroots, ...
    'VariableNames',{'node','status','n_roots'});
results.root_enumerator = struct('multiple',im,'touching',it,'discontinuity',idisc);
results.provenance = struct('test',mfilename,'legacy_fixture',fixture_path, ...
    'baseline_commit',F.provenance.commit);
diag_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'wp1_static_gate.mat');
save(diag_path,'results','-v7');
fprintf('test_invzp_static_domain: %d/%d checks passed; saved %s\n', ...
    nnz(passed),numel(passed),diag_path);

    function add_check(name, metric, tolerance, note)
        names(end+1,1) = name;
        metrics(end+1,1) = metric;
        tolerances(end+1,1) = tolerance;
        notes(end+1,1) = note;
    end
end

function [G0i, G0e] = static_split(si, T, Jxx0)
c0i = invz_chi0z(si, T, 0, struct('elastic',false));
G0i = -real(c0i(3,3,1));
c0 = invz_chi0z(si, T, 0, struct('elastic',true));
X = real(c0(:,:,1));
feedback = X(3,1)*(Jxx0/(1-Jxx0*X(1,1)))*X(1,3);
G0bare = -(X(3,3)+feedback);
G0e = G0bare-G0i;
end

function ok = legacy_scalar_bitwise_check()
stream = RandStream('mt19937ar','Seed',7);
ok = true;
for k = 1:100
    tl = struct('m',randn(stream),'M2',.1+rand(stream), ...
        'n01',.2+.8*rand(stream),'g0',randn(stream));
    lam = randn(stream,2,1)*.01;
    K = randn(stream)*.01; Sigma = randn(stream)*.02;
    beta = .1+100*rand(stream); G0i = -.1-10*rand(stream); G0e = randn(stream);
    xi = (1+tanh(tl.m^2*tl.n01^2*beta*K-tl.M2*beta*lam(1))) / ...
        (1+(4*tl.n01^2*K*tl.g0+2*lam(2)+tl.g0*lam(1))*tl.M2/tl.n01^2);
    local_denom = 1+Sigma+K*G0i;
    G = G0i/local_denom+xi*G0e;
    Gtil = G/(1-K*G);
    r = (G0i+G0e)/Gtil;
    [Gn,o] = invz_gstat_ordered(tl,lam,K,Sigma,beta,G0i,G0e);
    ok = ok && isequaln([Gn,o.xi,o.Gtil0,o.r,o.gstat_local_denom], ...
                        [G,xi,Gtil,r,local_denom]);
end
end
