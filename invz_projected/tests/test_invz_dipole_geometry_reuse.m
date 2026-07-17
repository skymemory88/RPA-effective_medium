function tests = test_invz_dipole_geometry_reuse
% Pins the q-independent geometry precompute in MF_dipole/exchange: the 5-arg
% call that reuses a prebuilt `geom` must be BIT-IDENTICAL to the 4-arg call
% that rebuilds the geometry from scratch, for every q. This guards the
% efficiency refactor (Codex finding #4) against any accidental numeric drift.
tests = functiontests(localfunctions);
end

function setupOnce(~)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..'));        % invz module
addpath(fullfile(here,'..','..'));   % repo root: MF_dipole, exchange, invz_ion
addpath(fullfile(here, '..', '..', 'invz_common'));
end

function qs = qset()
qs = [0 0 0
      0.5 0 0
      0.13 0.27 0.41
      1 0 0
      0.5 0.5 0.5
      0.25 0 0.1
     -0.3 0.2 -0.15];
end

function test_dipole_geom_reuse_bit_identical(testCase)
ion = invz_ion();
N = 8;                                          % small cutoff, fast
qs = qset();
[~, ~, geom] = MF_dipole(qs(1,:), N, ion.a, ion.tau);   % build once
for iq = 1:size(qs,1)
    q = qs(iq,:);
    [dip4, nN4] = MF_dipole(q, N, ion.a, ion.tau);        % 4-arg: rebuild
    [dip5, nN5] = MF_dipole(q, N, ion.a, ion.tau, geom);  % 5-arg: reuse
    verifyEqual(testCase, dip5, dip4, 'AbsTol', 0);       % bit-identical
    verifyEqual(testCase, nN5,  nN4,  'AbsTol', 0);
end
end

function test_exchange_geom_reuse_bit_identical(testCase)
ion = invz_ion();
qs = qset();
[~, geom] = exchange(qs(1,:), abs(ion.J12), ion.a, ion.tau);   % build once
for iq = 1:size(qs,1)
    q = qs(iq,:);
    d4 = exchange(q, abs(ion.J12), ion.a, ion.tau);            % 4-arg: rebuild
    d5 = exchange(q, abs(ion.J12), ion.a, ion.tau, geom);      % 5-arg: reuse
    verifyEqual(testCase, d5, d4, 'AbsTol', 0);                % bit-identical
end
end

function test_geom_is_q_independent(testCase)
% The returned geom must not depend on the q it was first built at.
ion = invz_ion();
N = 6;
[~, ~, gA] = MF_dipole([0 0 0],       N, ion.a, ion.tau);
[~, ~, gB] = MF_dipole([0.37 0.1 0.8], N, ion.a, ion.tau);
verifyEqual(testCase, gA.b,  gB.b,  'AbsTol', 0);
verifyEqual(testCase, gA.nN, gB.nN, 'AbsTol', 0);
verifyEqual(testCase, gA.r,  gB.r);
verifyEqual(testCase, gA.Tf, gB.Tf);

[~, xA] = exchange([0 0 0],        abs(ion.J12), ion.a, ion.tau);
[~, xB] = exchange([0.37 0.1 0.8], abs(ion.J12), ion.a, ion.tau);
verifyEqual(testCase, xA.b,    xB.b,    'AbsTol', 0);
verifyEqual(testCase, xA.r,    xB.r);
verifyEqual(testCase, xA.mask, xB.mask);
end
