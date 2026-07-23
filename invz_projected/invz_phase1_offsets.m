function offs = invz_phase1_offsets()
%INVZ_PHASE1_OFFSETS The eight {0,1/2}^3 BZ phase offsets (Phase-1 pre-registration
% docs/invzp_phase1_quadrature_prereg.md, "Offsets" -- "the eight {0,1/2}^3 phase offsets (eight
% including the 000 baseline)"), stage-2c coupling-only quadrature audit (additive, no field
% axis; see invz_phase1_qgrid.m for how these flags are turned into an actual q-point set).
%
% Returns a 1x8 struct array, deterministic order, with fields:
%   .tag    human-readable label: '000' (baseline, no shift) then one letter per shifted axis
%           ('h','k','l','hk','hl','kl','hkl') -- h/k/l name the FIRST/SECOND/THIRD reduced-q
%           component, matching qVec_generator's column order (qvec(:,1)=h, (:,2)=k, (:,3)=l).
%   .flags  [1x3] logical [oh ok ol]. flags(a)==false: axis a sits at phase 0 (unshifted).
%           flags(a)==true: axis a sits at phase 1/2, i.e. shifted by HALF of this convention's
%           own N-point axis step (see invz_phase1_qgrid.m header for the exact per-axis
%           construction -- offset [0 0 0], all flags false, needs no special-casing there: it
%           falls out of the SAME formula as a direct qVec_generator call, bit-for-bit).
%
% No RNG/Date/randperm anywhere -- calling this twice yields the identical struct array.
tags  = {'000', 'h', 'k', 'l', 'hk', 'hl', 'kl', 'hkl'};
flags = logical([0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1]);
offs = repmat(struct('tag', '', 'flags', false(1,3)), 1, 8);
for i = 1:8
    offs(i).tag   = tags{i};
    offs(i).flags = flags(i,:);
end
end
