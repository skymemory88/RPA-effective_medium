function tests = test_invzt_a3d_budget
%TEST_INVZT_A3D_BUDGET  7C Step-2: budget guards + dom_basis plumbing on the compact
% dominant vertex path of invzt_sigma_tensor. Error-path only (instant):
%   (1) every over-budget dominant call fails with invzt:orderedA3Budget BEFORE any
%       allocation (max_vertex_states / max_vertex_work / max_vertex_bytes);
%   (2) the guards are SCOPED to the dominant/compact path -- the legacy full-dress
%       path IGNORES the budget knobs;
%   (3) opts.dom_basis validates its E/p/Mz surface (invzt:domBasis).
% Runs with invz_projected OFF the path (CORE isolation).
tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));                          % invz_tensor
addpath(fullfile(here, '..', '..'));                    % repo root
addpath(fullfile(here, '..', '..', 'invz_common'));     % shared single-ion engine
ion = invz_ion();  T = 2.0;  B = [0.5 0 0];
si = invzt_threestate(ion, T, B, struct());
g6 = invzt_qgrid(6, 'halfopen');
lat = invzt_jq_tensor(ion, g6, struct('dpRng', 10, 'cache', true));
lat_eff = lat;  lat_eff.Jt = invzt_odd_mask(lat.Jt);    % odd = false (PM baseline)
[wn, ~, beta] = invz_matsubara(T, 40);
tc.TestData = struct('si', si, 'T', T, 'lat_eff', lat_eff, 'wn', wn, 'beta', beta);
end

function test_over_budget_work_errors(tc)
% 6*N^4*nwn*nl >> 1e3 -> refuse before building G4cc.
d = tc.TestData;
verifyError(tc, @() invzt_sigma_tensor(d.si, d.T, d.lat_eff, d.wn, d.beta, ...
    struct('dress', 'dominant', 'max_vertex_work', 1e3, 'max_outer', 3)), ...
    'invzt:orderedA3Budget');
end

function test_over_budget_states_errors(tc)
% dominant subspace has >= 2 states; max_vertex_states = 1 refuses it.
d = tc.TestData;
verifyError(tc, @() invzt_sigma_tensor(d.si, d.T, d.lat_eff, d.wn, d.beta, ...
    struct('dress', 'dominant', 'max_vertex_states', 1, 'max_outer', 3)), ...
    'invzt:orderedA3Budget');
end

function test_over_budget_bytes_errors(tc)
% compact [nwn,nl] + tile temporaries >> 1 byte -> refuse before allocation.
d = tc.TestData;
verifyError(tc, @() invzt_sigma_tensor(d.si, d.T, d.lat_eff, d.wn, d.beta, ...
    struct('dress', 'dominant', 'max_vertex_bytes', 1, 'max_outer', 3)), ...
    'invzt:orderedA3Budget');
end

function test_guards_scoped_to_dominant_only(tc)
% The legacy full-dress path must IGNORE the budget knobs: a tiny max_vertex_work that
% would refuse the dominant path leaves 'full' running (no invzt:orderedA3Budget).
% Small local Matsubara grid (Ecut = 6) so the full 81-combo full-dress build is fast.
d = tc.TestData;
[wn, ~, beta] = invz_matsubara(d.T, 6);
out = invzt_sigma_tensor(d.si, d.T, d.lat_eff, wn, beta, ...
    struct('dress', 'full', 'max_vertex_work', 1e3, 'max_outer', 1));
verifyTrue(tc, isstruct(out) && isequal(size(out.Vmat), [3 3 numel(wn)]));
end

function test_default_dominant_within_budget(tc)
% Sanity: the toy dominant path (Nd small) fits ALL default guards and runs.
d = tc.TestData;
out = invzt_sigma_tensor(d.si, d.T, d.lat_eff, d.wn, d.beta, ...
    struct('dress', 'dominant', 'max_outer', 1));
verifyTrue(tc, isstruct(out) && isequal(size(out.Vmat), [3 3 numel(d.wn)]));
end

function test_dom_basis_missing_field_errors(tc)
% opts.dom_basis must carry E, p AND Mz.
d = tc.TestData;
db = struct('E', [0; 0.1], 'p', [0.5; 0.5]);            % missing Mz
verifyError(tc, @() invzt_sigma_tensor(d.si, d.T, d.lat_eff, d.wn, d.beta, ...
    struct('dress', 'dominant', 'dom_basis', db, 'max_outer', 3)), 'invzt:domBasis');
end

function test_dom_basis_wrong_size_Mz_errors(tc)
% dom_basis.Mz must be numel(E) x numel(E).
d = tc.TestData;
db = struct('E', [0; 0.1], 'p', [0.5; 0.5], 'Mz', eye(3));   % 3x3 for a 2-state basis
verifyError(tc, @() invzt_sigma_tensor(d.si, d.T, d.lat_eff, d.wn, d.beta, ...
    struct('dress', 'dominant', 'dom_basis', db, 'max_outer', 3)), 'invzt:domBasis');
end
