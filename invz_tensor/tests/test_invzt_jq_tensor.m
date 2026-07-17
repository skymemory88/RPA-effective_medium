function tests = test_invzt_jq_tensor
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root: MF_dipole, exchange, qVec_generator
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
end

function test_qgrid_conventions(testCase)
g1 = invzt_qgrid(8, 'halfopen');
verifyEqual(testCase, size(g1.qvec, 1), 8^3 - 1);          % Gamma dropped
verifyEqual(testCase, sum(g1.w), 1, 'AbsTol', 1e-14);
verifyLessThan(testCase, max(abs(g1.qvec(:))), 1 - 1/16);  % half-open: no +0.5 face duplicate range
g2 = evalc_call(@() invzt_qgrid(8, 'legacy_inclusive'));
[qv, ~, ~] = evalc_call(@() qVec_generator(invz_ion().a, 'mode', 'grid', ...
    'grid', [8 8 8], 'range', [-0.5 0.5], 'verbose', false));
qv = qv(any(abs(qv) > 1e-12, 2), :);
verifyEqual(testCase, g2.qvec, qv, 'AbsTol', 1e-15);       % EXACT legacy reproduction
verifyFalse(testCase, strcmp(g1.hash, g2.hash));
end

function test_shape_hermiticity_conjsym(testCase)
ion = invz_ion();
q = [0.25 0 0; 0 0 0.25; 0.31 0.17 0.09; -0.25 0 0; 0 0 -0.25; -0.31 -0.17 -0.09];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 15, 'cache', false));
verifyEqual(testCase, size(lat.Jt), [12 12 6]);
for iq = 1:6
    verifyLessThan(testCase, norm(lat.Jt(:,:,iq) - lat.Jt(:,:,iq)', 'fro'), 1e-13);
end
for iq = 1:3
    verifyLessThan(testCase, norm(lat.Jt(:,:,iq+3) - conj(lat.Jt(:,:,iq)), 'fro'), 1e-12);
end
end

function test_gamma_lorentz_exact_placement(testCase)
% full - dipole = exchange + Lorentz. Subtract the INDEPENDENTLY built exchange
% (root function) and require the remainder == lorz on Cartesian diagonals of
% every pair, 0 elsewhere. (The old ">= lorz" form fails: J12 < 0 is AFM.)
ion = invz_ion();  C = invz_const();
lorz = 4*pi/(3*ion.Vc)*C.gfac;
latF = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false));
latD = invzt_jq_tensor(ion, [0 0 0], struct('dpRng', 15, 'cache', false, 'parts', 'dipole'));
[ex, ~] = exchange([0 0 0], abs(ion.J12), ion.a, ion.tau);
R = latF.Jt - latD.Jt;                                   % exchange + Lorentz
for s = 1:4, for sp = 1:4
    blk = R(3*(s-1)+(1:3), 3*(sp-1)+(1:3)) - sign(ion.J12)*ex(:,:,s,sp);
    verifyEqual(testCase, blk, lorz*eye(3), 'AbsTol', 1e-14);
end, end
end

function test_onaxis_smallq_odd_decay(testCase)
ion = invz_ion();  A = invzt_anchors();
q = [1e-1 0 0; 1e-2 0 0; 1e-3 0 0];
lat = invzt_jq_tensor(ion, q, struct('dpRng', 30, 'cache', false));
m = arrayfun(@(iq) max(abs(lat.Jt(3:3:12, 1:3:12, iq)), [], 'all'), 1:3);
verifyEqual(testCase, m(:), A.odd_onaxis_smallq.maxca(:), 'RelTol', 1e-6);
verifyEqual(testCase, m(2)/m(1), 0.1, 'RelTol', 0.3);    % ~linear decay
end

function test_cache_roundtrip_selfverifying(testCase)
ion = invz_ion();
q = [0.25 0 0; 0.31 0.17 0.09];
l1 = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', true));
l2 = invzt_jq_tensor(ion, q, struct('dpRng', 10, 'cache', true));
verifyEqual(testCase, l2.Jt, l1.Jt);
ion2 = ion;  ion2.J12 = ion.J12 * 1.05;
l3 = invzt_jq_tensor(ion2, q, struct('dpRng', 10, 'cache', true));
verifyFalse(testCase, isequal(l3.Jt, l1.Jt));
cdir = fullfile(fileparts(mfilename('fullpath')), '..', 'cache');
verifyTrue(testCase, ~isempty(dir(fullfile(cdir, 'jqt*.mat'))));
end

% ------- local helper (not a test: different signature than test_*(testCase)) -------
function varargout = evalc_call(fn)
%EVALC_CALL Run fn() with stdout captured/discarded; return its outputs.
%   Wraps MATLAB's [T, out1, ...] = evalc('expr') form so printing functions
%   (qVec_generator, etc.) stay quiet on the console while their return
%   values still propagate to the caller. fn must be a zero-argument handle.
n = max(nargout, 0);
outs = cell(1, n);
[~, outs{:}] = evalc('fn()');
varargout = outs;
end
