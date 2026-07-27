% Task 13 notes SS5: measure (not assume) the status at Bx = [9 0 0], T = 0.31.
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root, 'invz_projected'))); addpath(fullfile(root, 'invz_common'));
ion = invz_ion();
Jnu = linspace(-2e-3, 6.0e-3, 24).';
base = struct('J0eff', 6.42e-3, 'Jxx0', ion.Jxx0, 'hyp', true);
[h, p] = invz_hmf_ordered(ion, 0.31, [9 0 0], Jnu, base);
fprintf('Bx=[9 0 0]: hstar=%.6g status=%s\n', h, p.status);

% Also print the other two anchors for the same report, all three in one run.
[h0, p0] = invz_hmf_ordered(ion, 0.31, [0 0 0], Jnu, base);
fprintf('Bx=[0 0 0]: hstar=%.6g status=%s Delta(predictor-ish node1)=%.6g\n', h0, p0.status, p0.Delta(1));

o2 = base; o2.static_medium = 'strict_1z_dyson_ref'; o2.ref_margin = 1e9;
[h2, p2] = invz_hmf_ordered(ion, 0.31, [2.85 0 0], Jnu, o2);
fprintf('Bx=[2.85 0 0] ref_margin=1e9: hstar=%.6g status=%s (medium_status all ref_denom_small? %d)\n', ...
    h2, p2.status, all(strcmp(p2.medium_status, 'ref_denom_small')));
