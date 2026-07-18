%INVZT_RUN_SPECTRA  Full-tensor 1/z chi''_cc spectra vs transverse field, or along a q-path.
%
% SCOPE: PARAMAGNETIC SIDE ONLY. invzt_solve_point has no ordered-phase branch,
% so unlike invz_run_spectra this driver cannot sweep ACROSS Bc to show the soft
% mode hardening on both sides. Fields that land on the ordered side (or fail the
% sample-validity rule) are masked to NaN with a console note.
%
% solve_opts.mode must be 'a1' (enforced below AND by invzt_chi_realaxis's own
% invzt:realaxisMode guard): the real-axis continuation exists only for the A1
% scalar Sigma. There is consequently no 'dress' knob here -- it is A3-only.
%
% ERROR POLICY: field sweep -- per-field selective masking (physics signals and
% invalid samples become NaN columns with a console note; everything else
% rethrows). q-path view -- FAIL LOUD by design: the whole product hinges on the
% single Bq point, so an invalid point raises invzt:qpathNotPM (with converged/
% crit/Sigma0) instead of drawing an empty map, and any physics throw from the
% one solve propagates for the same reason.

addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'invz_common'));

ion = invz_ion();

% ---- knobs (mirroring invz_run_spectra's names where the concept carries over) --
T      = 1.6;                          % K
fields = linspace(0.3, 4.5, 40);       % T -- keep on the PM side (see SCOPE above)
w      = linspace(0, 0.6, 400).';      % meV, uniform spacing (invz_peak_energy assumes it)
eta    = 5e-3;                         % meV, real-axis Lorentzian HWHM
sliceMax  = 6;                         % field count at/below which line slices are drawn
peak_wmin = 0.02;                      % meV -- lower bound for the peak pick. NON-ZERO by
                                       % default to exclude the low-frequency hyperfine
                                       % line, mirroring invz_spectra_map's opts.peak_wmin;
                                       % set 0 to pick over the whole w grid.
theta_c = 0.0;  phi_ab = 0.0;          % field-direction tilt knobs (deg). Ported as-is:
                                       % invzt_solve_point already takes a full [Bx By Bz]
                                       % and forwards transverse_mf.
transverse_mf = 'legacy_x';            % 'legacy_x' | 'none' | 'vector_ab'
gridN = 16; gridConv = 'halfopen'; dpRng = 30;
useParallel = true;
solve_opts = struct('mode', 'a1', 'odd', true, 'nlevels', 'std', ...
                    'transverse_mf', transverse_mf);

% ---- q-path view: set qpath non-empty to switch views ------------------------
qpath = [];                            % [] = field sweep; [nq x 3] r.l.u. = q-path view.
                                       % Keep q-paths GAMMA-EXCLUDED: a strict q = [0 0 0]
                                       % is assembled with the Lorentz cavity (the strict-
                                       % uniform observable), NOT the q->0+ intrinsic limit
                                       % a dispersion plot wants -- start at finite q.
% qpath = [linspace(0.1, 0.5, 41).' zeros(41, 2)];   % example: toward (0.5,0,0)
Bq    = 2.0;                           % T -- fixed field for the q-path view
wq    = linspace(0, 0.6, 400).';       % meV -- own frequency grid for the q-path view
% -----------------------------------------------------------------------------

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
lat = invzt_jq_tensor(ion, g, struct('dpRng', dpRng, 'cache', true));

if ~isempty(qpath)
    % ---------------- q-path dispersion at one fixed field --------------------
    pt = invzt_solve_point(ion, T, Bq*dhat, lat, solve_opts);
    if ~(pt.converged && isfinite(pt.crit) && pt.crit > 0 && pt.Sigma0 >= sfloor)
        error('invzt:qpathNotPM', ['q-path point B = %.2f T, T = %.2f K is not a ' ...
            'valid PM sample (converged = %d, crit = %.4g, Sigma0 = %.4g): the whole ' ...
            'q-path product hinges on this one point, so failing loudly beats ' ...
            'drawing an empty map. Raise Bq, raise T, or check the knobs.'], ...
            Bq, T, pt.converged, pt.crit, pt.Sigma0);
    end
    out = invzt_chi_realaxis(ion, T, Bq*dhat, pt, wq, ...
            struct('qsel', qpath, 'dpRng', dpRng, 'eta', eta));
    chipp_q = imag(out.chi_cc_q);                 % [nq, nw] positive chi'' (Component 0)
    figure; imagesc(wq, 1:size(chipp_q, 1), chipp_q); set(gca, 'YDir', 'normal');
    xlabel('\omega (meV)'); ylabel('q index along path'); colorbar;
    title(sprintf('tensor 1/z \\chi''''_{cc}(q,\\omega), T = %.2f K, B = %.2f T', T, Bq));
    Epeak_q = invz_peak_energy(chipp_q.', wq, peak_wmin);   % columns must be per-q
    figure; plot(1:numel(Epeak_q), Epeak_q, 'o-');
    xlabel('q index along path'); ylabel('E_{peak} (meV)');
else
    % ---------------- field sweep at the uniform mode -------------------------
    nWorkers = 0;
    if useParallel && ~isempty(ver('parallel')), nWorkers = Inf; end
    nf = numel(fields);
    chipp = nan(numel(w), nf);
    parfor (k = 1:nf, nWorkers)
        col = nan(numel(w), 1);
        try
            pt = invzt_solve_point(ion, T, fields(k)*dhat, lat, solve_opts);
            % Same three-part sample-validity rule as invzt_critical (converged,
            % finite positive crit, single-sourced Sigma0 floor) -- so the
            % spurious negative-Sigma fixed point invzt_critical warns about
            % never reaches the plot.
            if pt.converged && isfinite(pt.crit) && pt.crit > 0 && pt.Sigma0 >= sfloor
                o = invzt_chi_realaxis(ion, T, fields(k)*dhat, pt, w, ...
                        struct('qsel', 'gamma_uniform', 'eta', eta));
                col = squeeze(imag(o.chi_uniform(3, 3, :)));   % already positive chi''
            else
                fprintf('  B = %.2f T: ordered / non-converged / invalid (masked)\n', fields(k));
            end
        catch err
            switch err.identifier
                case {'invz:degenerateDoublet', 'invz:orderedPhase', 'invzt:a1ZeroField'}
                    fprintf('  B = %.2f T: %s (masked)\n', fields(k), err.identifier);
                otherwise
                    rethrow(err);
            end
        end
        chipp(:, k) = col;
    end

    if nf <= sliceMax
        figure; hold on;  co = lines(nf);
        for k = 1:nf
            plot(w, chipp(:, k), '-', 'Color', co(k, :), ...
                 'DisplayName', sprintf('%.2f T', fields(k)));
        end
        xlabel('\omega (meV)'); ylabel('\chi''''_{cc}');
        title(sprintf('tensor 1/z, T = %.2f K', T)); legend show;
    else
        figure; imagesc(fields, w, chipp); set(gca, 'YDir', 'normal');
        xlabel('B (T)'); ylabel('\omega (meV)'); colorbar;
        title(sprintf('tensor 1/z \\chi''''_{cc}(B,\\omega), T = %.2f K', T));
    end

    Epeak = invz_peak_energy(chipp, w, peak_wmin);
    figure; plot(fields, Epeak, 'o-');
    xlabel('B (T)'); ylabel('E_{peak} (meV)');
    title(sprintf('\\chi''''_{cc} peak energy vs field, T = %.2f K', T));
    % Gaps are CENSORED peaks (boundary max / non-positive column) or masked
    % ordered points -- do not interpolate over them.
end
