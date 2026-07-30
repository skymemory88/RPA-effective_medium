function [Tc, out] = invz_odd_zero_field(ion, opts)
%INVZ_ODD_ZERO_FIELD Projected zero-field ordering temperature.
% opts.mode   'full' (static ODD) or 'off'
% opts.grids  grid sizes, default {12,24}; two sizes use linear 1/N extrapolation
% opts.dpRng  brute-force ODD cutoff, default 30
% opts.cache  optional geometric-block cache, default false
% opts.blocks optional precomputed single-grid Vca/Vcb/Vcc/Jcc0/Jaa0 struct
%
% Every internally built grid is zero-offset, half-open, with Gamma dropped.
if nargin < 2, opts = struct(); end
mode = char(getf(opts, 'mode', 'full'));
if ~ismember(mode, {'full','off'})
    error('invz:oddZeroFieldMode', 'mode must be ''full'' or ''off''.');
end
grids = getf(opts, 'grids', {12,24});
if iscell(grids), grids = cell2mat(grids); end
grids = grids(:).';
if isempty(grids) || any(grids < 2 | grids ~= round(grids))
    error('invz:oddZeroFieldGrid', 'grids must contain integer sizes >= 2.');
end
dpRng = getf(opts, 'dpRng', 30);
cache = getf(opts, 'cache', false);
prebuilt = getf(opts, 'blocks', []);
if ~isempty(prebuilt) && numel(grids) ~= 1
    error('invz:oddZeroFieldBlocks', ...
        'Precomputed blocks require exactly one grid size.');
end

ng = numel(grids);
Tcg = nan(1,ng);
Scg = nan(1,ng);
dg = zeros(1,ng);
tbuild = zeros(1,ng);
for k = 1:ng
    if isempty(prebuilt)
        g = invz_phase1_qgrid(ion, grids(k), [0 0 0], 'halfopen', 'P_drop');
        tb = tic;
        [Vca,Vcb,Vcc,info] = invz_odd_blocks( ...
            ion, g.qvec, struct('dpRng',dpRng,'cache',cache));
        tbuild(k) = toc(tb);
        B = struct('Vca',Vca,'Vcb',Vcb,'Vcc',Vcc, ...
            'Jcc0',info.Jcc0,'Jaa0',info.Jaa0);
    else
        B = prebuilt;
        if ~(isstruct(B) && isscalar(B) && ...
                all(isfield(B, {'Vca','Vcb','Vcc','Jcc0','Jaa0'})))
            error('invz:oddZeroFieldBlocks', ...
                'blocks must contain Vca/Vcb/Vcc/Jcc0/Jaa0.');
        end
    end

    if strcmp(mode, 'off')
        Jf = reshape(invz_odd_modes(B.Vcc, zeros(size(B.Vcc))), [], 1);
        Scg(k) = sigma_crit_quiet(B.Jcc0, Jf);
        Tcg(k) = invz_critical_T0field(ion, Scg(k), B.Jcc0);
    else
        [Tcg(k), Scg(k), dg(k)] = solve_full(ion, B);
    end
end

Tc = extrapolate(grids, Tcg);
out = struct('mode',mode, 'grids',grids, 'Tc',Tcg, 'Tc_rich',Tc, ...
    'Sc_at_Tc',Scg, 'Sc_rich',extrapolate(grids,Scg), ...
    'd_at_Tc',dg(end), 'timings',struct('build',tbuild));
end

function [Tc, Sc, d] = solve_full(ion, B)
memoT = NaN;
memoSc = NaN;
memoJ0 = NaN;
memoD = NaN;
Tc = invz_critical_T0field(ion, @sc_at, @j0_at);
refresh(Tc);
Sc = memoSc;
d = memoD;

    function value = sc_at(T)
        refresh(T);
        value = memoSc;
    end

    function value = j0_at(T)
        refresh(T);
        value = memoJ0;
    end

    function refresh(T)
        if T == memoT, return; end
        Xp = invz_chiperp(ion, T, 0, struct('Jxx0',B.Jaa0));
        [dJ,memoD] = invz_odd_deltaJ(B.Vca,B.Vcb,Xp);
        Jf = reshape(invz_odd_modes(B.Vcc,dJ), [], 1);
        memoJ0 = B.Jcc0 - memoD;
        memoSc = sigma_crit_quiet(memoJ0,Jf);
        memoT = T;
    end
end

function Sc = sigma_crit_quiet(J0,Jf)
state = warning('off','invz:sigmaCritExcluded');
cleanup = onCleanup(@() warning(state)); %#ok<NASGU>
Sc = invz_sigma_crit(J0,Jf);
end

function value = extrapolate(grids, values)
if isscalar(grids)
    value = values(1);
    return;
end
[g,idx] = sort(grids);
nc = g(1);
nf = g(end);
value = (nf*values(idx(end)) - nc*values(idx(1))) / (nf-nc);
end
