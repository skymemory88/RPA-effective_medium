function ion = invz_cfrot(ion, r_deg)
%INVZ_CFROT Rotate the m=4 crystal-field coefficient pairs by r_deg about c.
% Matches the external LiReF4 cf.m 'coefficient' convention (Br = [c4 s4; -s4 c4]
% on each [Bc; Bs] pair, angle 4*r): the rotated ion gains ion.B44s (sine partner
% of B44, honoured by invz_single_ion; invz_ion itself carries no B44s -> 0).
%
% Equivalence (pinned by test_invz_cfrot_equiv): rotating the CF by r with the
% field FIXED equals rotating the FIELD by phi_ab = +r with the CF fixed (SAME
% sign), through the whole scalar-cc pipeline PROVIDED transverse_mf =
% 'vector_ab' -- the cc coupling channel and the hyperfine term are rotation-
% invariant about c and the vector transverse MF is in-plane isotropic, so the
% CF is the only in-plane-anisotropic ingredient. legacy_x (x-only MF) breaks
% this covariance by construction.
c4 = cosd(4*r_deg);  s4 = sind(4*r_deg);
B44s0 = 0;  if isfield(ion, 'B44s'), B44s0 = ion.B44s; end
b44c =  c4*ion.B44  + s4*B44s0;
b44s = -s4*ion.B44  + c4*B44s0;
b64c =  c4*ion.B64c + s4*ion.B64s;
b64s = -s4*ion.B64c + c4*ion.B64s;
ion.B44 = b44c;  ion.B44s = b44s;  ion.B64c = b64c;  ion.B64s = b64s;
end
