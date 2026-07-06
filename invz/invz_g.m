function g = invz_g(tl, z)
%INVZ_G Two-level inelastic response g(z) = 2*n01*Delta/(Delta^2 - z^2)  (HTML eq 11).
g = 2*tl.n01*tl.Delta ./ (tl.Delta^2 - z.^2);
end
