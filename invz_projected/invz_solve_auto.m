function [pt, phase, diag] = invz_solve_auto(ion, T, Bx, Jnu_flat, opts)
%INVZ_SOLVE_AUTO Jensen-ordered solve with a paramagnetic fallback.
% phase: 1 = accepted ordered state, 2 = accepted stable paramagnet,
%        0 = no accepted state.
% opts.ordered_state_mode='linearized_pm_handoff' is an explicit hyp=false
% approximation: the PM 1/z mass selects the phase and its self-energy sets
% the modified ordered molecular field.  The default ordered-first Jensen
% route below is intentionally unchanged.
if nargin < 5, opts = struct(); end
ordered_state_mode = char(getf(opts, 'ordered_state_mode', 'jensen_integral'));
if strcmp(ordered_state_mode, 'linearized_pm_handoff')
    [pt, phase, diag] = linearized_pm_dispatch(ion, T, Bx, Jnu_flat, opts);
    return;
elseif ~strcmp(ordered_state_mode, 'jensen_integral')
    error('invz:orderedStateMode', ...
        ['ordered_state_mode must be ''jensen_integral'' or ' ...
         '''linearized_pm_handoff''.']);
end
pt = [];
phase = 0;
diag = struct('ordered_err', '', 'para_err', '');

try
    ordered = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts);
    if ordered.is_ordered && ordered.converged && isfinite(ordered.Sigma0)
        pt = ordered;
        phase = 1;
        return;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    diag.ordered_err = err.identifier;
end

try
    para = invz_solve_point(ion, T, Bx, Jnu_flat, opts);
    if para.converged && isfinite(para.Sigma0) && isfinite(para.crit) && para.crit > 0
        pt = para;
        phase = 2;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    diag.para_err = err.identifier;
end
end

function [pt, phase, diag] = linearized_pm_dispatch(ion, T, Bx, Jnu_flat, opts)
% PM-first stability dispatch matched to the projected 1/z critical mass.
if getf(opts, 'hyp', true)
    error('invz:hmfLinearizedHyperfine', ...
        ['linearized_pm_handoff is restricted to hyp=false so the existing ' ...
         'electro-nuclear Jensen calculation remains unchanged.']);
end
pt = [];
phase = 0;
diag = struct('ordered_err', '', 'para_err', '', ...
    'mode', 'linearized_pm_handoff', 'pm_accepted', false, ...
    'ordered_accepted', false, 'handoff_ratio', NaN, 'hmf_J0z', NaN);

% A field-map precomputes one boundary-exact coupling and a converged PM
% Sigma seed.  Fixed dispatch avoids asking the unstable PM iteration to
% reconverge independently at every ordered-side field.
fixed_handoff = all(isfield(opts, ...
    {'hmf_J0z', 'hmf_sigma0', 'hmf_boundary_field'}));
if fixed_handoff
    Bc = opts.hmf_boundary_field;
    if ~(isnumeric(Bc) && isreal(Bc) && isscalar(Bc) && isfinite(Bc))
        error('invz:hmfLinearizedBoundary', ...
            'hmf_boundary_field must be a finite real scalar.');
    end
    if Bx < Bc
        diag.handoff_ratio = opts.hmf_J0z/getf(opts, 'J0eff', ion.J0eff);
        diag.hmf_J0z = opts.hmf_J0z;
        try
            ordered = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, opts);
            pt = ordered;
            if ordered.is_ordered && ordered.converged && ...
                    isfinite(ordered.Sigma0) && isfinite(ordered.D_uni) && ...
                    ordered.D_uni > 0
                phase = 1;
                diag.ordered_accepted = true;
            end
        catch err
            if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
            diag.ordered_err = err.identifier;
        end
        return;
    end
    try
        para = invz_solve_point(ion, T, Bx, Jnu_flat, opts);
        pt = para;
        if para.converged && isfinite(para.Sigma0) && ...
                isfinite(para.crit) && para.crit > 0
            phase = 2;
            diag.pm_accepted = true;
        else
            diag.para_err = 'invz:unstableOrUnconvergedPm';
        end
    catch err
        if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
        diag.para_err = err.identifier;
    end
    return;
end

try
    para = invz_solve_point(ion, T, Bx, Jnu_flat, opts);
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    diag.para_err = err.identifier;
    return;
end
pt = para; % retain the failed/rejected producer for diagnostics
if ~(para.converged && isfinite(para.Sigma0) && isfinite(para.crit) && ...
        (1 + para.Sigma0) > 0)
    diag.para_err = 'invz:invalidPmHandoff';
    return;
end
if para.crit > 0
    phase = 2;
    diag.pm_accepted = true;
    return;
end

ratio = 1/(1 + para.Sigma0);
J0eff = getf(opts, 'J0eff', ion.J0eff);
J0z_mf = J0eff*ratio;
if ~(isfinite(ratio) && ratio > 0 && isfinite(J0z_mf) && J0z_mf > 0)
    diag.ordered_err = 'invz:invalidPmHandoff';
    return;
end
diag.handoff_ratio = ratio;
diag.hmf_J0z = J0z_mf;

oo = opts;
oo.hmf_J0z = J0z_mf;
oo.hmf_sigma0 = para.Sigma0;
if isfield(para, 'Sigma') && all(isfinite(para.Sigma))
    oo.Sigma_seed = para.Sigma;
end
try
    ordered = invz_solve_point_ordered(ion, T, Bx, Jnu_flat, oo);
    if ordered.is_ordered && ordered.converged && ...
            isfinite(ordered.Sigma0) && isfinite(ordered.D_uni) && ...
            ordered.D_uni > 0
        pt = ordered;
        phase = 1;
        diag.ordered_accepted = true;
    end
catch err
    if ~strncmp(err.identifier, 'invz:', 5), rethrow(err); end
    diag.ordered_err = err.identifier;
end
end

function v = getf(s, name, default)
if isfield(s, name), v = s.(name); else, v = default; end
end
