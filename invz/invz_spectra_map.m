function S = invz_spectra_map(ion, T, fields, w, opts)
%INVZ_SPECTRA_MAP chi''_cc(omega, B) at the uniform (q=0) mode over a transverse-field sweep.
%   S = invz_spectra_map(ion, T, fields, w) sweeps the transverse field Bx over `fields`
%   (Tesla) at fixed temperature T (K) and returns the dissipative susceptibility on the
%   (omega, B) grid -- covering BOTH phases. At each field it first tries the ferromagnetic
%   (ordered) 1/z solve; if the point has a spontaneous moment it uses that, otherwise it
%   falls back to the paramagnetic solve. The result is a single soft-mode map continuous
%   across the transition (imagesc(S.fields, S.w, S.chiz)).
%
%   Each column is one field. Returned maps:
%     S.chiz   [nw x nB]  1/z-renormalized chi''_cc  (FM below Bc, PM above)
%     S.chirpa [nw x nB]  bare-RPA chi''_cc          (Sigma = 0, matching phase)
%   Per-field diagnostics:
%     S.phase  [1 x nB]   0 = no solution (masked), 1 = ferromagnet, 2 = paramagnet
%     S.Sigma0 [1 x nB]   static self-energy at that field
%   Columns with no solution at all (e.g. the degenerate doublet at Bx -> 0, where the
%   two-level treatment breaks down in BOTH phases) are left NaN and masked out.
%
%   opts fields (all optional):
%     .hyp      (true)        include the nuclear hyperfine manifold in chi0
%     .grid     ([16 16 16])  q-grid for the effective-medium lattice sum
%     .dpRng    (30)          real-space dipole cutoff for invz_jq_modes
%     .eta      (5e-3)        real-axis broadening (meV) passed to invz_chi_realaxis
%     .parallel (false)       distribute the field points over a parallel pool (parfor)
%     .Jnu, .info             precomputed coupling branches / info struct (skips the
%                             qVec_generator + invz_jq_modes step; used by the tests)
%     .verbose  (true)        print one progress line per field
%
%   Cost is one 1/z solve per field (~15-25 min for a full 61-point sweep on the 16^3 grid,
%   single core; roughly 2x that near/below Bc where the FM solve runs). Compute S once, then
%   replot as many times as you like (see invz_plot_spectra_map / invz_run_spectra).

if nargin < 5, opts = struct(); end
hyp      = getf(opts, 'hyp', true);
grid     = getf(opts, 'grid', [16 16 16]);
dpRng    = getf(opts, 'dpRng', 30);
eta      = getf(opts, 'eta', 5e-3);
verbose  = getf(opts, 'verbose', true);
parallel = getf(opts, 'parallel', false);

fields = fields(:).';   w = w(:);
nB = numel(fields);     nw = numel(w);

if isfield(opts, 'Jnu') && isfield(opts, 'info')
    Jnu = opts.Jnu(:);   info = opts.info;
else
    [qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5]);
    qc = qc(any(abs(qc) > 1e-12, 2), :);
    [Jnu, info] = invz_jq_modes(ion, qc, struct('dpRng', dpRng, 'cache', true));
    Jnu = Jnu(:);
end
Jcc0 = info.Jcc0;

% Sliced plain-array outputs (parfor-safe); packed into S after the sweep.
chizM   = nan(nw, nB);   chirpaM = nan(nw, nB);
Sig0    = nan(1, nB);    phaseC  = zeros(1, nB);

% nWorkers = 0 forces serial execution even inside a parfor, and works without the
% Parallel Computing Toolbox; Inf lets parfor use (and auto-create) the pool.
nWorkers = 0;
if parallel && ~isempty(ver('parallel')), nWorkers = Inf; end

parfor (k = 1:nB, nWorkers)
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k)] = ...
        one_field(ion, T, fields(k), Jnu, Jcc0, w, eta, hyp);
    if verbose
        ph = {'ordered/degenerate (masked)', 'ferromagnet', 'paramagnet'};
        fprintf('  B = %5.2f T : %-28s Sigma0 = %s\n', fields(k), ph{phaseC(k)+1}, num2str(Sig0(k)));
    end
end

S = struct();
S.fields = fields;  S.w = w;  S.T = T;  S.info = info;
S.chiz = chizM;  S.chirpa = chirpaM;
S.Sigma0 = Sig0;  S.phase = phaseC;
end

% -------------------------------------------------------------------------------------------
function [chiz, chirpa, Sigma0, phase] = one_field(ion, T, B, Jnu, Jcc0, w, eta, hyp)
%ONE_FIELD chi''_cc(omega) at one field. Tries the ferromagnetic (ordered) solve first, then
% the paramagnet; returns phase = 1 (FM), 2 (PM), or 0 (no solution -> NaN columns).
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;

% --- ferromagnetic (ordered) branch ---
try
    pto = invz_solve_point_ordered(ion, T, B, Jnu, struct('hyp', hyp, 'J0eff', Jcc0));
    if pto.is_ordered && pto.converged && isfinite(pto.Sigma0)
        o  = invz_chi_realaxis(ion, T, B, pto, w, struct('Jsel', Jcc0, 'eta', eta));
        chiz = imag(o.chi_cc_q(1, :)).';
        pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pto.tl, ...
                     'K', [], 'is_ordered', true, 'si', pto.si);
        o0  = invz_chi_realaxis(ion, T, B, pt0, w, struct('Jsel', Jcc0, 'npass', 1, 'eta', eta));
        chirpa = imag(o0.chi_cc_q(1, :)).';
        Sigma0 = pto.Sigma0;  phase = 1;
        return;
    end
catch
    % fall through to the paramagnetic branch
end

% --- paramagnetic branch: bare-RPA overlay (needs only the two-level params) ---
try
    tl0 = invz_twolevel(ion, T, B);
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, struct('Jsel', Jcc0, 'npass', 1, 'eta', eta));
    chirpa = imag(o0.chi_cc_q(1, :)).';
catch
end
% --- paramagnetic 1/z solve ---
try
    pt = invz_solve_point(ion, T, B, Jnu, struct('hyp', hyp, 'J0eff', Jcc0));
    Sigma0 = pt.Sigma0;
    if pt.converged && isfinite(pt.Sigma0)
        o = invz_chi_realaxis(ion, T, B, pt, w, struct('Jsel', Jcc0, 'eta', eta));
        chiz = imag(o.chi_cc_q(1, :)).';
        phase = 2;
    end
catch
end
end

% -------------------------------------------------------------------------------------------
function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
