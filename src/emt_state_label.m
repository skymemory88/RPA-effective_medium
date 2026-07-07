function label = emt_state_label(scan_mode, cvar_value, dscrt_value)
% EMT_STATE_LABEL Build a readable state descriptor for logs.

dscrt = NaN;
if ~isempty(dscrt_value)
    dscrt = dscrt_value(1);
end

if strcmp(scan_mode, 'field')
    label = sprintf('B = %.4f T, T = %.4f K', cvar_value, dscrt);
elseif strcmp(scan_mode, 'temp')
    label = sprintf('T = %.4f K, B = %.4f T', cvar_value, dscrt);
else
    label = sprintf('state = %.4f', cvar_value);
end

end
