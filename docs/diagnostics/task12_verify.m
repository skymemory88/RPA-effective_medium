% Task 12 controller verification.
%  (1) Is the DOMAIN record branch really unreachable? hmf_ordered:71 says opts.ref_margin is
%      resolved per call, and Task 10's test forced ref_denom_small with ref_margin = 1e9.
%  (2) Independently reproduce the strict-vs-resummed hstar shift (the headline physics finding).
%  (3) Confirm the two path integrals differ substantively (the spec's binding caution).
root = '/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/Programming scripts/Matlab/Simulation/invZ expansion';
addpath(root);
addpath(fullfile(root, 'invz_projected'));
addpath(fullfile(root, 'invz_common'));

ion = invz_ion();
T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
J0eff = 6.42e-3;
base = struct('J0eff', J0eff, 'Jxx0', ion.Jxx0, 'hyp', true);

fprintf('===== (2) strict vs resummed =====\n');
for scheme = {'resummed', 'strict_1z_dyson_ref'}
    o = base;  o.static_medium = scheme{1};
    [hs, p] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
    fprintf('VER %-20s status=%-12s hstar=%.17g\n', scheme{1}, p.status, hs);
    fprintf('VER %-20s slope0=%.17g  crit_star=%.17g\n', scheme{1}, p.slope0, p.crit_star);
    fprintf('VER %-20s int_Sigma0=%.17g  int_r_minus_1=%.17g  ratio=%.6g\n', scheme{1}, ...
            p.int_Sigma0, p.int_r_minus_1, p.int_r_minus_1/p.int_Sigma0);
    fprintf('VER %-20s nNodes=%d n_extend=%d redensified=%d\n', scheme{1}, ...
            numel(p.hgrid), p.n_extend, p.redensified);
    if ~strcmp(scheme{1}, 'resummed')
        fprintf('VER %-20s omit_max range=[%.6g %.6g]  omit_cubic range=[%.6g %.6g]\n', ...
                scheme{1}, min(p.omit_max), max(p.omit_max), ...
                min(p.omit_cubic), max(p.omit_cubic));
        us = unique(p.medium_status);
        fprintf('VER %-20s medium_status set = %s\n', scheme{1}, strjoin(us, ','));
    end
end

fprintf('\n===== (1) is the DOMAIN branch reachable via ref_margin? =====\n');
o = base;  o.static_medium = 'strict_1z_dyson_ref';  o.ref_margin = 1e9;
try
    [hs2, p2] = invz_hmf_ordered(ion, T, Bx, Jnu, o);
    fprintf('VER refmargin1e9 status=%s  hstar=%.17g\n', p2.status, hs2);
    us2 = unique(p2.medium_status);
    fprintf('VER refmargin1e9 medium_status set = %s\n', strjoin(us2, ','));
    fprintf('VER refmargin1e9 ref_denom(1)=%.17g  ref_margin(1)=%.17g\n', ...
            p2.ref_denom(1), p2.ref_margin(1));
    nbad = sum(~strcmp(p2.medium_status, 'ok') & ~strcmp(p2.medium_status, 'not_applicable'));
    fprintf('VER refmargin1e9 nodes with a NON-ok medium_status = %d / %d\n', ...
            nbad, numel(p2.medium_status));
    if nbad > 0
        fprintf('VER *** DOMAIN BRANCH IS REACHABLE -- the gap is closable ***\n');
    else
        fprintf('VER *** domain branch NOT reached even at ref_margin=1e9 ***\n');
    end
    fprintf('VER refmargin1e9 node_term_reason set = %s\n', ...
            strjoin(unique(p2.node_term_reason), ','));
catch err
    fprintf('VER refmargin1e9 THREW %s : %s\n', err.identifier, err.message);
end
fprintf('\nVER DONE\n');
