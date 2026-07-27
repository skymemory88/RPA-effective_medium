ROOT = ['/Users/yikaiyang/Library/CloudStorage/OneDrive-Nexus365/' ...
        'Programming scripts/Matlab/Simulation/invZ expansion'];
cd(ROOT);
addpath(genpath(fullfile(ROOT, 'invz_projected')));
addpath(fullfile(ROOT, 'invz_common'));

ion = invz_ion();
T = 0.31;

% --- Invalid path: degenerate doublet, domain_policy = 'return' ---
o_invalid = struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x', 'domain_policy', 'return');
tl_invalid = invz_twolevel_ordered(ion, T, [0 0 0], 0, o_invalid);
fn_invalid = sort(fieldnames(tl_invalid));
fprintf('REVIEWER invalid tl fields = %s\n', strjoin(fn_invalid, ','));
fprintf('REVIEWER invalid tl.valid = %d, tl.Delta = %.15e\n', tl_invalid.valid, tl_invalid.Delta);

% Does a consumer that ignores .valid and reads a two-level field fail loudly?
try
    x = tl_invalid.g0; %#ok<NASGU>
    fprintf('REVIEWER UNEXPECTED: tl.g0 accessible on invalid tl, value=%g\n', x);
catch err
    fprintf('REVIEWER g0-access-on-invalid-tl THROWS: id=%s msg=%s\n', err.identifier, err.message);
end

% --- Valid path: same point used by test_valid_flag_present_on_the_normal_path ---
o_valid = struct('Jxx0', ion.Jxx0, 'transverse_mf', 'legacy_x');
tl_valid = invz_twolevel_ordered(ion, T, [2.85 0 0], 0.02, o_valid);
fn_valid = sort(fieldnames(tl_valid));
fprintf('REVIEWER valid tl fields = %s\n', strjoin(fn_valid, ','));
fprintf('REVIEWER valid tl.valid = %d, tl.Delta = %.15e\n', tl_valid.valid, tl_valid.Delta);

% --- domain_policy validation on an otherwise-valid node: does it fire BEFORE any
%     use of tl.Delta, i.e. is it unconditional? ---
try
    invz_twolevel_ordered(ion, T, [2.85 0 0], 0.02, ...
        struct('Jxx0', ion.Jxx0, 'domain_policy', 123));
    fprintf('REVIEWER UNEXPECTED: numeric domain_policy silently accepted\n');
catch err
    fprintf('REVIEWER numeric-policy THROWS: id=%s msg=%s\n', err.identifier, err.message);
end

% --- single-diagonalization sanity: wrap invz_single_ion call count via a global counter
%     is intrusive; instead just confirm invz_single_ion itself does not recurse into
%     invz_twolevel_ordered (grep-verified already) and that both policy branches return
%     promptly (no retry loop) by checking elapsed calls is trivially fast -- informational only.
fprintf('REVIEWER DONE\n');
