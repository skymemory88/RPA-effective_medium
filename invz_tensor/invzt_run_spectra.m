%INVZT_RUN_SPECTRA  Full-tensor 1/z chi''_cc spectra vs transverse field, or along a q-path.
%
% SCOPE: BOTH PHASES, ACROSS Bc. Each point is solved with invzt_solve_auto
% (Task 5), which runs the PM and ordered legs and assigns the phase by
% STABILITY, not moment onset: the PM leg's three-part validity rule
% (converged, finite crit > 0, Sigma0 >= floor) decides FIRST -- its crit > 0
% IS the tensor QPT criterion -- and the ordered leg is consulted only when
% the PM sample is invalid. Unlike its PM-only predecessor, this driver
% therefore sweeps straight across Bc and shows the soft mode hardening on
% BOTH sides. Fields where neither leg is accepted (near-Bc Option-A window,
% or a solver failure) are masked to NaN with a structured console note.
%
% TRANSVERSE ONLY: theta_c (the field's tilt out of the ab-plane) must stay at
% 0 -- a nonzero Bz routes into invzt:orderedLongitudinal, rethrown by
% invzt_solve_auto rather than absorbed, since no tensor forced-moment route
% exists here (2026-07-19 scope decision).
%
% solve_opts.mode must be 'a1' (enforced below AND by invzt_chi_realaxis's own
% invzt:realaxisMode guard): the real-axis continuation exists only for the A1
% scalar Sigma. There is consequently no 'dress' knob here -- it is A3-only.
%
% ERROR POLICY: field sweep -- per-field selective masking (physics signals and
% invalid samples become NaN columns with a structured per-leg console note;
% everything else rethrows). q-path view -- FAIL LOUD by design: the whole
% product hinges on the single Bq point, so a point where both legs fail
% raises invzt:qpathInvalid (with both legs' converged/crit/m0/err) instead of
% drawing an empty map, and any physics throw from either solve propagates for
% the same reason.

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs (mirroring invz_run_spectra's names where the concept carries over) --
T = 0.1;                            % K
fields = linspace(4.0, 6.0, 101);   % T -- SPANS the QPT (corrected QPT lies between 4.64 and 4.65 T at
                                    % 0.1 K): stability-based auto solve assigns FM below,
                                    % PM above (invzt_solve_auto; NOT the bare-MF moment,
                                    % which persists to ~5.0 T -- review P0-1)
w = linspace(0, 0.018, 401).';      % meV, uniform (invz_peak_energy assumes it). This LOW-energy
                                    % window targets the experimentally observed electro-nuclear
                                    % soft branch (roughly 0-4.35 GHz). The predominantly electronic
                                    % excitation at O(0.3 meV) is a separate branch and is deliberately
                                    % outside this view; it need not soften to zero after hybridization.
eta = 5e-5;                         % meV, real-axis Lorentzian HWHM. MUST stay >= the w step
                                    % (4.5e-5 here) or peaks alias between grid points (measured
                                    % 10x peak-height error at eta/step = 0.16); guarded below.
sliceMax = 6;                       % field count at/below which line slices are drawn
peak_wmin = 0;                      % meV -- include the electro-nuclear mode all the way to Bc.
                                    % invz_peak_energy censors a maximum in the first bin, so an
                                    % exactly critical zero-energy pole is NaN rather than a fake gap.
theta_c = 0.0;  phi_ab = 0.0;       % tilt knobs (deg). theta_c ~= 0 gives Bz ~= 0 ->
                                    % invzt:orderedLongitudinal (no tensor forced-moment
                                    % route; 2026-07-19 scope).
transverse_mf = 'legacy_x';         % 'legacy_x' | 'none' | 'vector_ab'
gridN = 16; gridConv = 'halfopen';
dipoleBackend = 'ewald';            % certified production default; 'bruteforce' remains diagnostic
ewaldOpts = invzt_ewald_defaults(ion);
dpRng = 30;                         % used only by explicit 'bruteforce'
useParallel = true;
solve_opts = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', ...
                    'transverse_mf', transverse_mf);

if eta < (w(2) - w(1))
    error('invzt_run_spectra:etaStep', ['eta = %.3g meV < w step = %.3g meV: a ' ...
        'Lorentzian narrower than the sampling aliases between grid points ' ...
        '(measured 10x peak-height error at eta/step = 0.16, 2026-07-19). Raise ' ...
        'eta or refine w.'], eta, w(2) - w(1));
end

% ---- q-path view: set qpath non-empty to switch views ------------------------
qpath = [];                            % [] = field sweep; [nq x 3] r.l.u. = q-path view.
                                       % Keep q-paths GAMMA-EXCLUDED: a strict q = [0 0 0]
                                       % is assembled with the Lorentz cavity (the strict-
                                       % uniform observable), NOT the q->0+ intrinsic limit
                                       % a dispersion plot wants -- start at finite q.
% qpath = [linspace(0.1, 0.5, 41).' zeros(41, 2)];   % example: toward (0.5,0,0)
Bq    = 2.0;                           % T -- fixed field for the q-path view
wq = w;
% wq    = linspace(0, 0.6, 400).';       % meV -- own frequency grid for the q-path view
% -----------------------------------------------------------------------------

if eta < (wq(2) - wq(1))
    error('invzt_run_spectra:etaStep', ['eta = %.3g meV < wq step = %.3g meV: a ' ...
        'Lorentzian narrower than the sampling aliases between grid points ' ...
        '(measured 10x peak-height error at eta/step = 0.16, 2026-07-19). Raise ' ...
        'eta or refine wq.'], eta, wq(2) - wq(1));
end

if ~strcmp(getf(solve_opts, 'mode', 'a1'), 'a1')
    error('invzt_run_spectra:mode', ['invzt_run_spectra requires solve_opts.mode = ' ...
        '''a1'' (got ''%s''): invzt_chi_realaxis is the A1 scalar-Sigma continuation ' ...
        'ONLY. A2/A3 real-axis continuation is an open item (README section 10).'], ...
        char(getf(solve_opts, 'mode', 'a1')));
end
if phi_ab ~= 0 && strcmp(transverse_mf, 'legacy_x')
    error('invzt_run_spectra:transverseMF', ...
        ['phi_ab = %.3g deg needs the vector transverse mean field: set transverse_mf ' ...
         'to ''vector_ab'' (or ''none'' for a bare CF+Zeeman diagnostic). legacy_x is ' ...
         'x-only and C4-inconsistent for rotated fields.'], phi_ab);
end
if ~isempty(qpath)
    % F8 convention preflight (R4, 2026-07-18 second Codex review): dispersion
    % q-paths must exclude strict Gamma (and any Gamma-equivalent row) -- see
    % the qpath knob comment above. Runs BEFORE any lattice/solve work.
    if ~(isnumeric(qpath) && isreal(qpath) && ismatrix(qpath) && size(qpath, 2) == 3 ...
            && all(isfinite(qpath(:))))
        error('invzt_run_spectra:qpath', ['qpath must be an [nq,3] real finite numeric ' ...
            'array of r.l.u. coordinates (see the qpath knob comment above); got a %s ' ...
            '%s.'], class(qpath), mat2str(size(qpath)));
    end
    if any(invz_is_gamma_equiv(qpath, ion.tau))
        error('invzt_run_spectra:qpathGamma', ['qpath must exclude Gamma-equivalent rows: ' ...
            'a strict q = [0 0 0] (or any Gamma-equivalent point) is assembled with the ' ...
            'Lorentz cavity -- the strict-uniform observable, NOT the q->0+ intrinsic limit ' ...
            'a dispersion plot wants (see the qpath knob comment above). Start the path at ' ...
            'finite q.']);
    end
end
dhat   = [cosd(theta_c)*cosd(phi_ab), cosd(theta_c)*sind(phi_ab), sind(theta_c)];
sfloor = getf(solve_opts, 'sigma_floor', -0.5);   % single-sourced with invzt_critical

g   = invzt_qgrid(gridN, gridConv);
latOpts = spectra_lattice_opts(dipoleBackend, ewaldOpts, dpRng);
lat = invzt_jq_tensor(ion, g, latOpts);

if ~isempty(qpath)
    % ---------------- q-path dispersion at one fixed field --------------------
    [pt, phaseq, diq] = invzt_solve_auto(ion, T, Bq*dhat, lat, solve_opts);
    if phaseq == 0
        error('invzt:qpathInvalid', ['q-path point B = %.2f T, T = %.2f K failed both ' ...
            'legs [PM: conv=%d crit=%.4g err=%s | ORD: conv=%d m0=%.4g crit=%.4g err=%s]: ' ...
            'the whole q-path product hinges on this one point, so failing loudly beats ' ...
            'drawing an empty map.'], Bq, T, diq.para.converged, diq.para.crit, ...
            diq.para.err, diq.ordered.converged, diq.ordered.m0, diq.ordered.crit, ...
            diq.ordered.err);
    end
    out = invzt_chi_realaxis(ion, T, Bq*dhat, pt, wq, ...
            struct('qsel', qpath, 'eta', eta));
    chipp_q = imag(out.chi_cc_q);                 % [nq, nw] positive chi'' (Component 0)
    finiteMask = isfinite(chipp_q);
    Z = log10(max(chipp_q, realmin));
    figure;
    im = imagesc(wq, 1:size(chipp_q, 1), Z);
    set(im, 'AlphaData', double(finiteMask));
    set(gca, 'YDir', 'normal', 'Color', [0.8 0.8 0.8], 'Layer', 'top');
    pos = chipp_q(finiteMask & chipp_q > 0);
    if ~isempty(pos)
        ps = sort(pos(:));
        hi = ps(max(1, min(numel(ps), ceil(0.995*numel(ps)))));
        clim([log10(hi/1e3) log10(hi)]);
    end
    colormap(turbo);
    xlabel('\omega (meV)'); ylabel('q index along path');
    cb = colorbar;  cb.Label.String = 'log_{10} \chi''''_{cc}';
    title(sprintf('tensor 1/z \\chi''''_{cc}(q,\\omega), T = %.2f K, B = %.2f T (phase %d)', ...
        T, Bq, phaseq));
    Epeak_q = invz_peak_energy(chipp_q.', wq, peak_wmin);   % columns must be per-q
    figure; plot(1:numel(Epeak_q), Epeak_q, 'o-');
    xlabel('q index along path'); ylabel('E_{peak} (meV)');
    title(sprintf('q-path E_{peak}, T = %.2f K, B = %.2f T (phase %d)', T, Bq, phaseq));
else
    % ---------------- field sweep at the uniform mode, ACROSS Bc ---------------
    nWorkers = 0;
    if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end
    nf = numel(fields);
    chipp  = nan(numel(w), nf);
    phasev = zeros(1, nf);          % 1 = ordered, 2 = PM, 0 = masked
    critv  = nan(1, nf);            % branch stability of the ACCEPTED leg (review re-plan #4)
    parfor (k = 1:nf, nWorkers)
        col = nan(numel(w), 1);  ph = 0;  ck = NaN;
        try
            [pt, ph, di] = invzt_solve_auto(ion, T, fields(k)*dhat, lat, solve_opts);
            if ph > 0
                ck = pt.crit;
                o = invzt_chi_realaxis(ion, T, fields(k)*dhat, pt, w, ...
                        struct('qsel', 'gamma_uniform', 'eta', eta));
                col = squeeze(imag(o.chi_uniform(3, 3, :)));
            else
                % Structured per-leg outcomes (P2-1): a leg that RETURNED but was
                % rejected is described by its numbers, not a blank error string.
                fprintf(['  B = %.2f T: masked -- PM(att=%d conv=%d crit=%.4g err=%s) ' ...
                         'ORD(att=%d conv=%d m0=%.4g crit=%.4g err=%s)\n'], fields(k), ...
                    di.para.attempted, di.para.converged, di.para.crit, di.para.err, ...
                    di.ordered.attempted, di.ordered.converged, di.ordered.m0, ...
                    di.ordered.crit, di.ordered.err);
            end
        catch err
            switch err.identifier
                case {'invz:degenerateDoublet', 'invzt:a1ZeroField'}
                    fprintf('  B = %.2f T: %s (masked)\n', fields(k), err.identifier);
                otherwise
                    rethrow(err);
            end
        end
        chipp(:, k) = col;  phasev(k) = ph;  critv(k) = ck;
    end
    fprintf('phase summary: %d ordered / %d PM / %d masked of %d fields\n', ...
        sum(phasev == 1), sum(phasev == 2), sum(phasev == 0), nf);

    if nf <= sliceMax
        figure; hold on;  co = lines(nf);
        for k = 1:nf
            plot(w, chipp(:, k), '-', 'Color', co(k, :), ...
                 'DisplayName', sprintf('%.2f T', fields(k)));
        end
        xlabel('\omega (meV)'); ylabel('\chi''''_{cc}');
        title(sprintf('tensor 1/z, T = %.2f K', T)); legend show;
    else
        % Match invz_projected's spectrum-map display: log10(chi'') over the
        % three decades below the robust 99.5th-percentile peak. The full
        % tensor data contains weak electronuclear satellites at 1e-3 and
        % below relative weight; a raw linear colour scale hides them behind
        % the bright collective soft pole even though the poles are present.
        finiteMask = isfinite(chipp);
        Z = log10(max(chipp, realmin));
        figure;
        im = imagesc(fields, w, Z);
        set(im, 'AlphaData', double(finiteMask));
        set(gca, 'YDir', 'normal', 'Color', [0.8 0.8 0.8], 'Layer', 'top');
        pos = chipp(finiteMask & chipp > 0);
        if ~isempty(pos)
            ps = sort(pos(:));                     % toolbox-free robust percentile
            hi = ps(max(1, min(numel(ps), ceil(0.995*numel(ps)))));
            lo = hi / 1e3;
            clim([log10(lo) log10(hi)]);
        end
        colormap(turbo);
        xlabel('B (T)'); ylabel('\omega (meV)');
        title(sprintf('tensor 1/z \\chi''''_{cc}(B,\\omega), T = %.2f K, across B_c', T));
        cb = colorbar;  cb.Label.String = 'log_{10} \chi''''_{cc}';
    end

    Epeak = invz_peak_energy(chipp, w, peak_wmin);
    figure; plot(fields, Epeak, 'o-');
    xlabel('B (T)'); ylabel('E_{peak} (meV)');
    title(sprintf('\\chi''''_{cc} peak energy vs field, T = %.2f K', T));
    if any(phasev == 1) && any(phasev == 2), xline(fields(find(phasev == 2, 1)), '--', 'B_c'); end
    % Gaps are CENSORED peaks (boundary max / non-positive column) or masked
    % ordered points -- do not interpolate over them.
end

function opts = spectra_lattice_opts(backend, ewaldOpts, dpRng)
opts = struct('dipole', backend, 'cache', true);
if strcmp(backend, 'ewald')
    opts.ewald = ewaldOpts;
elseif strcmp(backend, 'bruteforce')
    opts.dpRng = dpRng;
else
    error('invzt_run_spectra:dipoleBackend', ...
        'dipoleBackend must be ''ewald'' or ''bruteforce'' (got ''%s'').', backend);
end
end
