function scan = invz_static_domain_scan(ion, T, Bx, Jnu_flat, opts)
%INVZ_STATIC_DOMAIN_SCAN Prospective (stage-a) half of Gate 0: does the strict-order static
% medium have a valid reference along the H_MF profile grid, and how large are the omitted
% moment terms there? G = -chi, ferromagnetic positive J.
%
% BARE DIAGONALIZATIONS ONLY -- no Sigma<->medium solve, no iteration, no convergence
% dependence. Sigma0 is taken as 0 (the kinematic proxy): Sigma0 is O(1/z), so this bounds the
% reference denominator 1 + Sigma0 to leading order. The SOLVED-path margins are the stage-b
% half of Gate 0 and are read off prof.medium_status / prof.Delta after the real solve -- this
% function does NOT predict them (spec SS7.2).
%
% GRID OWNERSHIP (prereg SS7.2): opts.hgrid is honoured verbatim when supplied. Otherwise the
% INITIAL geometric grid is built from the shared invz_hmf_grid helper with the same hmax rule
% invz_hmf_ordered uses (hmax_fac * |bare ordered hz|, or the exact opts.hmax_abs override).
% This function deliberately does NOT reproduce invz_hmf_ordered's adaptive extension or
% redensification: two implementations that agree in one test are not grid identity.
%
% opts: J0eff (required), Jxx0, transverse_mf, hyp, nH, hmax_fac, hmax_abs, hmin_frac, hgrid,
%       static_medium (must be a strict scheme), ref_margin, Jmom, Ecut.
% scan: hgrid, Delta, valid, G0bare, Gref, ref_denom, ref_status{}, omit_mu3, omit_cubic,
%       omit_max, predictor, n_nodes, n_required, n_valid, n_skipped, n_out_of_domain,
%       n_degenerate, scheme, grid_source.
% The h=0 predictor is required by Jensen but is not part of the nonzero geometric hgrid, so
% it is evaluated and accounted separately. Counters satisfy
% n_valid+n_degenerate+n_out_of_domain+n_skipped == n_required.
if nargin < 5, opts = struct(); end
J0eff  = opts.J0eff;
Jxx0   = getf(opts, 'Jxx0', ion.Jxx0);
tmf    = getf(opts, 'transverse_mf', 'legacy_x');
hyp    = getf(opts, 'hyp', true);
nH     = getf(opts, 'nH', 33);
hfrac  = getf(opts, 'hmin_frac', 1e-3);
Ecut   = getf(opts, 'Ecut', 40);
scheme = getf(opts, 'static_medium', '');
if ~any(strcmp(scheme, {'strict_1z_dyson_ref', 'strict_1z_bare_ref'}))
    error('invz:staticMedium', ['invz_static_domain_scan requires a strict scheme ' ...
        '(''strict_1z_dyson_ref'' or ''strict_1z_bare_ref''); got ''%s''. The resummed medium ' ...
        'has no reference denominator to scan.'], scheme);
end
refopt = struct('ref_margin', getf(opts, 'ref_margin', 1e-6));
Jmom = getf(opts, 'Jmom', []);
if isempty(Jmom), Jmom = invz_coupling_moments(Jnu_flat); end
if ~isvector(Jnu_flat)
    error('invz:staticMedium', 'ordered-domain scan requires a static vector Jnu_flat.');
end

sibase = struct('hyp', hyp, 'Jxx0', Jxx0, 'transverse_mf', tmf);
if isfield(opts, 'hgrid') && ~isempty(opts.hgrid)
    hgrid = opts.hgrid(:).';  grid_source = 'explicit';
else
    sibo = sibase;  sibo.order = true;  sibo.J0z = J0eff;
    sib = invz_single_ion(ion, T, Bx, sibo);
    if isfield(opts, 'hmax_abs'), hmax = opts.hmax_abs;
    else,                          hmax = getf(opts, 'hmax_fac', 1.25) * abs(sib.hz);  end
    if ~(isfinite(hmax) && hmax > 0)
        error('invz:hmfGrid', ['no positive bracket ceiling: the bare set does not order at ' ...
            'T = %g, |Bx| = %g. Supply opts.hgrid to scan anyway.'], T, norm(Bx));
    end
    hgrid = invz_hmf_grid(hmax, nH, hfrac);  grid_source = 'invz_hmf_grid';
end

n = numel(hgrid);
empty = struct('status','not_evaluated','Delta',NaN,'G0bare',NaN,'Gref',NaN, ...
    'ref_denom',NaN,'ref_margin',NaN,'omit_mu3',NaN,'omit_cubic',NaN,'omit_max',NaN);
nodes = repmat(empty, 1, n);
for k = 1:n
    nodes(k) = eval_proxy(hgrid(k));
end
predictor = eval_proxy(0);                 % the same routine, not a second implementation

node_status = {nodes.status};
statuses = [{predictor.status}, node_status];
known = strcmp(statuses, 'ok') | strcmp(statuses, 'degenerate_doublet') | ...
        strcmp(statuses, 'ref_denom_nonpositive') | strcmp(statuses, 'ref_denom_small') | ...
        strcmp(statuses, 'nonfinite');
n_valid         = nnz(strcmp(statuses, 'ok'));
n_degenerate    = nnz(strcmp(statuses, 'degenerate_doublet'));
n_skipped       = nnz(strcmp(statuses, 'not_evaluated'));
n_out_of_domain = nnz(known) - n_valid - n_degenerate;
n_required      = n + 1;
if nnz(known) + n_skipped ~= n_required
    error('invz:staticDomainScan', 'unclassified prospective node status.');
end

scan = struct('hgrid', hgrid, 'Delta', [nodes.Delta], ...
    'valid', strcmp({nodes.status}, 'ok'), 'G0bare', [nodes.G0bare], ...
    'Gref', [nodes.Gref], 'ref_denom', [nodes.ref_denom], ...
    'ref_margin', [nodes.ref_margin], 'ref_status', {node_status}, ...
    'predictor', predictor, 'omit_mu3', [nodes.omit_mu3], ...
    'omit_cubic', [nodes.omit_cubic], 'omit_max', [nodes.omit_max], ...
    'n_nodes', n, 'n_required', n_required, 'n_valid', n_valid, ...
    'n_skipped', n_skipped, 'n_out_of_domain', n_out_of_domain, ...
    'n_degenerate', n_degenerate, 'scheme', scheme, 'grid_source', grid_source);

    function rec = eval_proxy(hp)
    % One bare, non-iterative kinematic proxy. Every return preserves the fixed schema.
    rec = empty;
    sio = sibase;  sio.hz_fixed = hp;
    si = invz_single_ion(ion, T, Bx, sio);
    tl = invz_twolevel_ordered(ion, T, Bx, hp, ...
        struct('Jxx0', Jxx0, 'transverse_mf', tmf, 'domain_policy', 'return'));
    rec.Delta = tl.Delta;
    if ~tl.valid
        rec.status = 'degenerate_doublet';
        return;
    end
    wn0 = invz_matsubara(T, Ecut);
    c0 = invz_chi0z(si, T, 1i*wn0(1), struct('elastic', true));
    X = real(c0(:, :, 1));
    % THIRD VERBATIM COPY of the transverse-MF backaction block already at
    % invz_projected/invz_hmf_ordered.m:490-505 and
    % invz_projected/invz_solve_point_ordered.m:265-278 (switch tmf / fb /
    % G0bare = -(X(3,3)+fb)). Arithmetically identical to both existing copies. This is exactly
    % the quantity Gate 0 compares, so if either of those two sites is ever edited or extracted
    % into a shared helper, this copy MUST move with it -- all three stay in lockstep.
    switch tmf
        case 'none',      fb = 0;
        case 'legacy_x',  fb = X(3,1) * (Jxx0/(1-Jxx0*X(1,1))) * X(1,3);
        case 'vector_ab', t = [1 2];
                          fb = X(3,t) * (Jxx0*((eye(2)-Jxx0*X(t,t))\X(t,3)));
        otherwise, error('invz:transverseMF', 'unknown transverse_mf ''%s''', tmf);
    end
    rec.G0bare = -(X(3,3) + fb);
    % Sigma0 = 0 is hard-wired below (the kinematic proxy, see the file header). For
    % 'strict_1z_dyson_ref' that makes denom = 1 + 0 = 1; 'strict_1z_bare_ref' is denom = 1 by
    % construction regardless. So ref_denom is IDENTICALLY 1 for BOTH strict schemes here,
    % Gref == G0bare identically, and 'ref_denom_nonpositive'/'ref_denom_small' are UNREACHABLE
    % from this proxy (only a non-finite G0bare/denom can route status away from 'ok'). This scan
    % therefore cannot validate the Dyson reference denominator -- its real content is the
    % Delta-domain map (rec.Delta/tl.valid above) and the omitted-moment ratios below.
    [rec.Gref, ref] = invz_static_medium_reference(rec.G0bare, 0, scheme, refopt);
    rec.ref_denom = ref.denom;  rec.ref_margin = ref.margin;  rec.status = ref.status;
    if ~strcmp(ref.status, 'ok'), return; end
    [~, clo] = invz_medium_moment_closure(rec.Gref, Jmom, scheme);
    rec.omit_mu3 = clo.omit_mu3;  rec.omit_cubic = clo.omit_cubic;
    rec.omit_max = clo.omit_max;  rec.status = clo.status;
    end
end
