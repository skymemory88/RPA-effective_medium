function report = validate_jq_input(Jq)
% VALIDATE_JQ_INPUT Validate shape and symmetry quality of J(q).

report = struct();
report.valid = true;
report.shape = size(Jq);

if ndims(Jq) ~= 3 && ndims(Jq) ~= 4
    report.valid = false;
    report.message = 'Jq must be 3-D or 4-D.';
    return;
end

if size(Jq, 1) ~= 3 || size(Jq, 2) ~= 3
    report.valid = false;
    report.message = 'First two dimensions of Jq must be 3x3.';
    return;
end

if ndims(Jq) == 3
    n_q = size(Jq, 3);
    n_c = 1;
else
    s = size(Jq);
    n_q = max(s(3), s(4));
    n_c = min(s(3), s(4));
end

report.n_q = n_q;
report.n_cvar_like = n_c;

max_asym = 0;
max_imag_diag = 0;

if ndims(Jq) == 3
    for iq = 1:size(Jq, 3)
        J = Jq(:,:,iq);
        max_asym = max(max_asym, norm(J - J', 'fro'));
        max_imag_diag = max(max_imag_diag, max(abs(imag(diag(J)))));
    end
else
    for i3 = 1:size(Jq, 3)
        for i4 = 1:size(Jq, 4)
            J = Jq(:,:,i3,i4);
            max_asym = max(max_asym, norm(J - J', 'fro'));
            max_imag_diag = max(max_imag_diag, max(abs(imag(diag(J)))));
        end
    end
end

report.max_asymmetry = max_asym;
report.max_imag_diag = max_imag_diag;
report.message = 'ok';

end
