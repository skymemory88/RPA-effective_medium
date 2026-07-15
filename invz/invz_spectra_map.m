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
%     When ion.demag ~= 0, both maps are the strict-uniform MEASURED observable
%     (demag-corrected via info.Jshape_cc; the soft mode saturates at 1/Jshape_cc
%     instead of diverging); with demag = 0 they are the intrinsic response.
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
    [qc, ~, ~] = qVec_generator(ion.a, 'mode', 'grid', 'grid', grid, 'range', [-0.5 0.5], 'verbose', false);
    qc = qc(any(abs(qc) > 1e-12, 2), :);
    [Jnu, info] = invz_jq_modes(ion, qc, struct('dpRng', dpRng, 'cache', true));
    Jnu = Jnu(:);
end
Jcc0 = info.Jcc0;
Jaa0   = ion.Jxx0;  if isfield(info, 'Jaa0'),      Jaa0   = info.Jaa0;      end
Jshape = 0;         if isfield(info, 'Jshape_cc'), Jshape = info.Jshape_cc; end

% Sliced plain-array outputs (parfor-safe); packed into S after the sweep.
chizM   = nan(nw, nB);   chirpaM = nan(nw, nB);
Sig0    = nan(1, nB);    phaseC  = zeros(1, nB);

% nWorkers = 0 forces serial execution even inside a parfor, and works without the
% Parallel Computing Toolbox; Inf lets parfor use (and auto-create) the pool.
nWorkers = 0;
if parallel && ~isempty(ver('parallel')), nWorkers = Inf; end

parfor (k = 1:nB, nWorkers)
    [chizM(:, k), chirpaM(:, k), Sig0(k), phaseC(k)] = ...
        one_field(ion, T, fields(k), Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp);
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
function [chiz, chirpa, Sigma0, phase] = one_field(ion, T, B, Jnu, Jcc0, Jaa0, Jshape, w, eta, hyp)
%ONE_FIELD chi''_cc(omega) at one field via the shared ordered-first solve
% (invz_solve_auto); returns phase = 1 (FM), 2 (PM), or 0 (no solution -> NaN columns).
% Jsel = Jcc0 is the strict-uniform observable, so the demag correction Jshape applies.
nw = numel(w);
chiz = nan(nw, 1);  chirpa = nan(nw, 1);  Sigma0 = NaN;  phase = 0;
sopts = struct('hyp', hyp, 'J0eff', Jcc0, 'Jxx0', Jaa0);
copts = struct('Jsel', Jcc0, 'eta', eta, 'Jxx0', Jaa0, 'Jshape', Jshape, 'hyp', hyp);

[pt, phase] = invz_solve_auto(ion, T, B, Jnu, sopts);

if phase == 1                                     % --- ferromagnetic (ordered) branch ---
    o  = invz_chi_realaxis(ion, T, B, pt, w, copts);   % reuses pt.si (ordered eigenstates)
    chiz = imag(o.chi_cc_q(1, :)).';
    pt0 = struct('alpha', 0, 'alpha_m', 0, 'lambda', [0; 0; 0], 'tl', pt.tl, ...
                 'K', [], 'is_ordered', true, 'si', pt.si);
    c0opts = copts;  c0opts.npass = 1;  c0opts.chi0cc_w = o.chi0cc_w;   % share bare cc
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
    Sigma0 = pt.Sigma0;
    return;
end

% --- paramagnetic side: bare-RPA overlay first (needs only the two-level params), so a
% non-converged 1/z point still gets its RPA column. invz:degenerateDoublet is the one
% expected condition here (Bx -> 0); anything else is a defect and propagates.
% When the 1/z point converged (phase 2) reuse its two-level / single-ion state and the
% bare chi0cc, so the overlay and the 1/z call don't each rebuild the diagonalization.
reuse = phase == 2 && ~isempty(pt);
if reuse, tl0 = pt.tl;  si0 = pt.si;  else, tl0 = invz_twolevel(ion, T, B, struct('Jxx0', Jaa0));  si0 = []; end
chi0cc = [];
try
    pt0 = struct('alpha', 0, 'lambda', [0; 0], 'tl', tl0, 'K', []);
    c0opts = copts;  c0opts.npass = 1;  c0opts.si = si0;
    o0  = invz_chi_realaxis(ion, T, B, pt0, w, c0opts);
    chirpa = imag(o0.chi_cc_q(1, :)).';
    chi0cc = o0.chi0cc_w;                          % share the bare cc with the 1/z call
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
end
if ~isempty(pt) && isfield(pt, 'Sigma0'), Sigma0 = pt.Sigma0; end
if phase == 2                                     % --- converged paramagnetic 1/z ---
    copts1 = copts;
    if ~isempty(chi0cc), copts1.chi0cc_w = chi0cc; end
    o = invz_chi_realaxis(ion, T, B, pt, w, copts1);
    chiz = imag(o.chi_cc_q(1, :)).';
end
end
