function [si, phase, massPm, out] = invz_bare_mf_state( ...
    ion, T, Bx, J0eff, Jxx0, hyp)
%INVZ_BARE_MF_STATE Select the stable bare-MF PM or ordered single-ion state.
%   [si,phase,massPm,out] = invz_bare_mf_state(ion,T,Bx,J0eff,Jxx0,hyp)
%   first evaluates the paramagnetic state and its uniform longitudinal mass
%
%       massPm = 1 - J0eff*chi0cc(0).
%
%   massPm >= 0 selects the PM state (phase=2). On the ordered side, the
%   nonzero root of hz-J0eff*<Jz>(hz)=0 is found from fixed-hz single-ion
%   evaluations and a sign bracket (phase=1). This avoids the critical
%   slowing of direct ordered Picard iteration without changing the MF
%   equation or accepting a loose last iterate.

if nargin < 6, hyp = true; end
if ~(isnumeric(Bx) && isreal(Bx) && isscalar(Bx) && isfinite(Bx))
    error('invz:bareMFField', 'Bx must be a finite real scalar.');
end
if ~(isnumeric(J0eff) && isreal(J0eff) && isscalar(J0eff) && isfinite(J0eff))
    error('invz:bareMFCoupling', 'J0eff must be a finite real scalar.');
end
if ~(isnumeric(Jxx0) && isreal(Jxx0) && isscalar(Jxx0) && isfinite(Jxx0))
    error('invz:bareMFCoupling', 'Jxx0 must be a finite real scalar.');
end
if ~(isscalar(hyp) && (islogical(hyp) || ...
        (isnumeric(hyp) && isreal(hyp) && isfinite(hyp) && any(hyp == [0 1]))))
    error('invz:bareMFHyperfine', 'hyp must be a finite logical scalar.');
end
hyp = logical(hyp);

pmOpts = struct('hyp', hyp, 'Jxx0', Jxx0);
siPm = invz_single_ion(ion, T, [Bx 0 0], pmOpts);
if ~siPm.mf_converged
    error('invz:bareMFNotConverged', ...
        'Bare-MF transverse PM mean field did not converge at B = %.12g T.', Bx);
end
cstat = invz_chi0z(siPm, T, 0, struct('elastic', true));
chi0cc0 = real(cstat(3,3,1));
massPm = 1 - J0eff*chi0cc0;
if ~isfinite(massPm)
    error('invz:bareMFMass', 'Bare-MF PM mass is non-finite at B = %.12g T.', Bx);
end

out = struct('phase', 2, 'mass_pm', massPm, 'chi0cc0_pm', chi0cc0, ...
    'method', "pm_mass", 'root', NaN, 'root_residual', NaN, ...
    'root_bracket', [NaN NaN], 'bracket_halvings', 0);
if massPm >= 0
    si = siPm;
    phase = 2;
    return;
end

% The PM mass says that hz=0 is unstable. Find the nonzero symmetry-broken
% root without iterating toward the nearby trivial root.
hhi = 1.05*abs(J0eff)*ion.J;
[fhi, ~] = bare_mf_residual(hhi);
if ~(isfinite(fhi) && fhi > 0)
    error('invz:bareMFBracket', ...
        'Could not establish the upper bare-MF bracket at B = %.12g T.', Bx);
end

b = hhi;
fb = fhi;
bracketed = false;
for k = 1:80
    a = 0.5*b;
    [fa, ~] = bare_mf_residual(a);
    if isfinite(fa) && fa < 0
        bracketed = true;
        break;
    end
    b = a;
    fb = fa;
end
if ~bracketed || ~(isfinite(fb) && fb >= 0)
    error('invz:bareMFBracket', ...
        ['PM mass predicts bare order at B = %.12g T, but the nonzero ' ...
         'hz=J0eff*<Jz> root could not be bracketed.'], Bx);
end

rootOpts = optimset('Display', 'off', 'TolX', 1e-12);
hz = fzero(@(h) bare_mf_residual(h), [a b], rootOpts);
[froot, si] = bare_mf_residual(hz);
if ~(isfinite(froot) && abs(froot) < 1e-10)
    error('invz:bareMFNotConverged', ...
        'Bare-MF ordered root residual is %.6g meV at B = %.12g T.', froot, Bx);
end
si.bare_mf_residual = froot;
si.rpa_mf_residual = froot; % compatibility with the former RPA-local solver
phase = 1;
out.phase = phase;
out.method = "bracketed_nonzero_root";
out.root = hz;
out.root_residual = froot;
out.root_bracket = [a b];
out.bracket_halvings = k;

    function [f, state] = bare_mf_residual(hzTrial)
        state = invz_single_ion(ion, T, [Bx 0 0], ...
            struct('hyp', hyp, 'Jxx0', Jxx0, 'hz_fixed', hzTrial));
        if ~state.mf_converged
            error('invz:bareMFNotConverged', ...
                ['Bare-MF fixed-hz transverse mean field failed at ' ...
                 'B = %.12g T, hz = %.6g meV.'], Bx, hzTrial);
        end
        f = hzTrial - J0eff*state.Jexp(3);
    end
end
