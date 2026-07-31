function result = test_invzp_hmf_derivative_identity()
%TEST_INVZP_HMF_DERIVATIVE_IDENTITY Guard the exact equation-(45) derivative.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(fullfile(repo,'invz_common'));

tl = struct('m',1.2,'M2',2.3,'n01',0.97,'g0',0.8);
lambda = [0.003;-0.002];
K0 = 0.0013;
Sigma0 = -0.1;
beta = 20;
G0inel0 = -120;
G0el0 = -40;
J0eff = 0.006421661809416939;
[Gstat,o] = invz_gstat_ordered(tl,lambda,K0,Sigma0,beta, ...
    G0inel0,G0el0,struct('stable_form',true));

direct = o.r+J0eff*o.G0bare;
exact = o.r*(1+J0eff*o.Gtil0);
wrong = o.r*(1+J0eff*Gstat);
identity_error = abs(direct-exact);
wrong_gap = abs(wrong-exact);
assert(identity_error < 1e-12, ...
    'F'' direct and Gtilde forms must agree algebraically.');
assert(wrong_gap > 1e-3, ...
    'The fixture must distinguish Gstat from Gtilde in F''.');

result = struct('direct',direct,'exact_Gtilde',exact, ...
    'incorrect_Gstat',wrong,'identity_error',identity_error, ...
    'incorrect_substitution_gap',wrong_gap);
fprintf(['test_invzp_hmf_derivative_identity: direct %.12g, ' ...
    'Gtilde %.12g, Gstat-substitution %.12g\n'],direct,exact,wrong);
end
