function d = invz_exact_numeric_digest(x)
%INVZ_EXACT_NUMERIC_DIGEST Exact-byte SHA-256 of a numeric array's class, shape and data.
% ONE definition, THREE consumers: Task 0 freezes the coupling fingerprint, Task 17's G11
% anchor asserts it, and Task 18's Gate-0 driver re-checks it before any solve. Reproducing
% this algorithm from prose in three places would let a silent mismatch destroy the
% fingerprint's whole purpose and surface as a spurious Gate-0 abort.
%
% Digest input, in this exact order: the class name, the size vector as doubles, then the
% ORIGINAL class-specific bytes of the data in column-major order. Deliberately NOT
% order-invariant and NOT shape-invariant: a reshaped or reordered multiset is a different
% input and must digest differently.
%
% The data bytes are taken with typecast(x(:),'uint8'), NEVER via double(x). Converting first
% would make the digest exact only for doubles: int64(2^53) and int64(2^53+1) are distinct
% inputs that collapse to the same double, so a double-converting digest would call them equal
% (pinned by test_digest_hashes_original_class_bytes_not_doubles). The size vector is genuinely
% double because size() returns doubles; that is a canonical encoding, not a data conversion.
% Contract: FULL, REAL, numeric x. Sparse and complex inputs are rejected rather than silently
% densified or split into re/im, since either choice would be an unstated encoding decision.
%
% JVM REQUIREMENT: uses java.security.MessageDigest, so it does not work under -nojvm. This
% matches the existing precedent and constraint documented for invz_bz_couplings.m's local
% qhash helper; this repository is plotting-oriented and -nojvm is not a supported mode.
if ~isnumeric(x) || ~isreal(x) || issparse(x)
    error('invz:exactDigest', ['x must be a full real numeric array; got %s%s.'], ...
          class(x), repmat(' (sparse)', 1, issparse(x)));
end
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(uint8(class(x)));
md.update(typecast(double(size(x)), 'uint8'));
md.update(typecast(x(:), 'uint8'));      % ORIGINAL class bytes -- never via double(x)
h = typecast(md.digest(), 'uint8');
d = lower(reshape(dec2hex(h, 2).', 1, []));
end
