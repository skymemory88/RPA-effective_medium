function result = invzp_m2_asymptotic_check()
%INVZP_M2_ASYMPTOTIC_CHECK Condition the removable M2 factor in ordered Sigma.
% The test freezes a healthy production state and varies only M2. This is a
% numerical conditioning check, not a claim that M2 varies independently in
% the physical model.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(here)));
addpath(repo,fullfile(repo,'invz_common'),fullfile(repo,'invz_projected'));

fixture_path = fullfile(repo,'docs','diagnostics','invzp_static_wp1', ...
    'legacy_4T_fixture.mat');
census_path = fullfile(here,'wp2_low_field_m2_census.mat');
F = load(fixture_path);
C = load(census_path);
pt = F.pt;
[wn,~,beta] = invz_matsubara(F.provenance.T,F.provenance.solve_opts.Ecut);
g = real(invz_g(pt.tl,1i*wn));

M2 = unique([pt.tl.M2; C.result.summary.min_M2; ...
    logspace(-2,-300,60).'; 1e-305; 1e-307; 1e-308; 1e-309; ...
    1e-310; realmin; 0],'stable');
n = numel(M2);
current_finite = false(n,1);
max_abs_delta = nan(n,1);
max_rel_delta = nan(n,1);
current_norm = nan(n,1);
stable_norm = nan(n,1);

for k = 1:n
    tl = pt.tl;
    tl.M2 = M2(k);
    sc = invz_sigma_ordered(tl,pt.lambda,pt.K,g,beta);
    ss = stable_sigma_ordered(tl,pt.lambda,pt.K,g,beta);
    current_finite(k) = all(isfinite(sc.Sigma));
    current_norm(k) = max(abs(sc.Sigma));
    stable_norm(k) = max(abs(ss.Sigma));
    if current_finite(k)
        max_abs_delta(k) = max(abs(sc.Sigma-ss.Sigma));
        max_rel_delta(k) = max_abs_delta(k)/max(max(abs(ss.Sigma)),realmin);
    end
end

tab = table(M2,current_finite,current_norm,stable_norm,max_abs_delta, ...
    max_rel_delta);
positive = M2 > 0;
first_nonfinite = find(positive & ~current_finite,1);
finite_rows = positive & current_finite;
result = struct('table',tab,'healthy_fixture_M2',pt.tl.M2, ...
    'low_field_min_M2',min(C.result.summary.min_M2), ...
    'max_relative_delta_while_finite',max(max_rel_delta(finite_rows)), ...
    'largest_nonfinite_M2',max(M2(positive & ~current_finite),[],'omitnan'), ...
    'ratio_overflow_M2',2*pt.tl.m^2/realmax, ...
    'first_nonfinite_index',first_nonfinite, ...
    'stable_M2_zero_norm',stable_norm(M2==0), ...
    'provenance',struct('fixture',fixture_path,'low_field_census',census_path, ...
    'commit',git_commit(repo)), ...
    'note',['With finite m,n01,Delta,lambda,K: alpha and gamma are O(M2), ' ...
    'Q=(2m^2/M2)gamma0 has the finite form ' ...
    '2m^2[lambda1-(1-n01^2)K0]/n01^2, and Sigma tends to ' ...
    '-alpha_m-Q*g. The production expression becomes nonfinite only near ' ...
    'floating-point underflow or at M2=0; compare that scale with the ' ...
    'measured low-field minimum before considering an arithmetic rewrite.']);
save(fullfile(here,'wp2_m2_asymptotic_check.mat'),'result','-v7');

fprintf('healthy M2 %.6g, low-field minimum %.6g\n', ...
    result.healthy_fixture_M2,result.low_field_min_M2);
fprintf('max relative current/stable delta while finite %.3g\n', ...
    result.max_relative_delta_while_finite);
fprintf(['direct ratio overflows below M2 %.3g in this frozen state; ' ...
    'largest sampled nonfinite M2 %.3g; stable ||Sigma(0)|| %.6g\n'], ...
    result.ratio_overflow_M2,result.largest_nonfinite_M2, ...
    result.stable_M2_zero_norm);
end

function sig = stable_sigma_ordered(tl,lam,K,g,beta)
% Algebraically equivalent to invz_sigma_ordered for M2>0, with the only
% M2 cancellation performed before floating-point evaluation.
s0 = invz_sigma(tl,lam,K,g,beta);
K = K(:);
g = g(:);
K0 = K(1);
g0 = tl.g0;
alpha_m = (tl.m^2/tl.n01^2) * (lam(2)-g0*lam(1)+ ...
    (4/g0)*lam(3)-(1-tl.n01^2)* ...
    (1+0.5*beta*tl.Delta*tl.n01)*K0*g0);
Q = (2*tl.m^2/tl.n01^2) * ...
    (lam(1)-(1-tl.n01^2)*K0);
sig = struct('alpha',s0.alpha,'alpha_m',alpha_m, ...
    'gamma',s0.gamma,'Q',Q, ...
    'Sigma',(s0.alpha-alpha_m)+(s0.gamma-Q).*g);
end

function commit = git_commit(repo)
[status,text] = system(sprintf('git -C "%s" rev-parse HEAD',repo));
if status == 0, commit = strtrim(text); else, commit = 'unknown'; end
end
