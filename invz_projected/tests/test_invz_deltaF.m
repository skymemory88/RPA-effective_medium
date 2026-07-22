function tests = test_invz_deltaF
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
addpath(fullfile(here, '..', '..'));
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function test_bare_limit_deltaF_zero(testCase)
% Route A with force_bare: dh = h0 - hmf = 0 identically, so BOTH dF and the
% non-claiming endpoint diagnostic are exactly zero (P0-3: the round-1 draft asserted
% tail > 0, which contradicts dh == 0). Field names per the dF_partial contract
% (stage-2 task 6b, step 1): out.dF_partial, out.endpoint_dh.
ion = invz_ion();  T = 0.31;
Jnu = linspace(-2e-3, 6.0e-3, 24).';
[dF, out] = invz_deltaF_ordered(ion, T, [2.85 0 0], Jnu, ...
    struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true, 'force_bare', true));
verifyEqual(testCase, dF, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, out.dF_partial, 0, 'AbsTol', 1e-12);
verifyEqual(testCase, out.endpoint_dh, 0, 'AbsTol', 1e-15);
end

function test_two_routes_agree(testCase)
% Framework SS9.4: route A (field dependence of Sigma) vs route B (frequency structure,
% J 2.22 temperature-integrated) must agree -- the stringent global check.
%   dF(T) = + T * int_T^inf dU(T')/T'^2 dT'      (sign derived in plan SS7, P0-3)
% Tolerance 10% relative with BOTH truncation tails required small first; tighten grids
% before tolerances. FAILURE BEYOND TOLERANCE IS A BLOCKED ESCALATION (SS9).
assumeTrue(testCase, ~isempty(getenv('INVZ_SLOW')), 'INVZ_SLOW=1 to run (T-grid sweep).');
ion = invz_ion();  T = 0.31;  Bx = [2.85 0 0];
Jnu = linspace(-2e-3, 6.0e-3, 24).';
o = struct('J0eff', 6.4e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
[dFA, outA] = invz_deltaF_ordered(ion, T, Bx, Jnu, o);
verifyTrue(testCase, isfinite(dFA));
verifyLessThan(testCase, outA.tail_est, 0.05*abs(dFA));      % route-A tail must be small
% P2-G: tail_est is only honest if route A is CONVERGED in the saturation cutoff --
% doubling hmax_fac must not move dF beyond 5% relative.
[dFA2, ~] = invz_deltaF_ordered(ion, T, Bx, Jnu, setfield(o, 'hmax_fac', 8)); %#ok<SFLD>
verifyEqual(testCase, dFA, dFA2, 'RelTol', 0.05);
% Route B on a geometric T grid; dU per J 2.22 with G reconstructed via J 2.30, on the
% SOLVER'S OWN Matsubara grid rebuilt per temperature (round-2 P0-B, made concrete in
% round 3: invz_matsubara exists -- verify the solver calls it, then reuse).
Tg = T * (1.35.^(0:17));  dU = nan(size(Tg));        % 18 points: last 3 are the Tmax check
for k = 1:numel(Tg)
    pt = invz_solve_point(ion, Tg(k), Bx, Jnu, o);
    if ~pt.converged, continue; end
    [wnk, wtsk, betak] = invz_matsubara(Tg(k), getf(o, 'Ecut', 40));
    c0k = invz_chi0z(pt.si, Tg(k), 1i*wnk, struct('elastic', true));
    G0k = -real(squeeze(c0k(3,3,:)));
    Gk  = G0k ./ (1 + pt.Sigma + pt.K .* G0k);       % J 2.30
    tlk = pt.tl;
    dU(k) = 0.5*( pt.alpha*tlk.n01*tlk.Delta/(1 + pt.alpha) - tlk.M2*pt.lambda(1) ...
                  + real(sum(wtsk .* pt.K .* (Gk - G0k)))/betak );
end
fin = isfinite(dU);
verifyGreaterThanOrEqual(testCase, nnz(fin), 15);
Tf = Tg(fin);  Uf = dU(fin);
n15 = nnz(Tf <= T*1.35^14);                          % the shorter grid's reach
dFB15 = + T * trapz(Tf(1:n15), Uf(1:n15) ./ Tf(1:n15).^2);   % CORRECTED SIGN (P0-3)
dFB   = + T * trapz(Tf, Uf ./ Tf.^2);
% round-3 P2-7: the tail is an EMPIRICAL estimate, not a proven bound -- so require
% Tmax convergence explicitly (extending Tmax by 1.35^3 moves dFB < 5%), and report
% tail_est_B alongside.
verifyEqual(testCase, dFB15, dFB, 'RelTol', 0.05);   % Tmax-converged
tail_est_B = + T * abs(Uf(end)) / Tf(end);           % empirical remainder estimate
verifyLessThan(testCase, tail_est_B, 0.05*max(abs(dFA), abs(dFB)));
verifyEqual(testCase, dFA, dFB, 'RelTol', 0.10, ...
    sprintf('two-route dF mismatch (A=%.4g, B=%.4g, tailA~%.2g, tailB~%.2g) -- BLOCKED escalation', ...
            dFA, dFB, outA.tail_est, tail_est_B));
end
