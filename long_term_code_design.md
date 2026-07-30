# Long-Term Code Design

This document records design decisions and implementation proposals that should
survive individual debugging sessions. A proposal is not necessarily implemented;
each section states its current status.

## Hybrid RPA Mean-Field Branch Selection: Continuation with Independent Stability Validation

### Status and recommendation

The projected-spin spectra currently use an independent, pointwise bare-MF/RPA
selector implemented by the local functions `bare_rpa_point`,
`solve_bare_ordered`, and `bare_mf_residual` in
`invz_projected/invz_chi_realaxis.m`. A spectra-driver run after this repair
was reported to produce a physically correct RPA soft mode. Preserve that
implementation until the validation gates below have been run. The hybrid scheme
described here is a long-term robustness improvement, not an immediate
replacement demanded by the present result.

The recommended design combines:

1. ordered-to-paramagnetic continuation, which provides a smooth and efficient
   initial guess along a field sweep; and
2. an independent static-stability test, which remains the authority for the
   phase and prevents scan direction or solver history from selecting an
   unstable branch.

Continuation is therefore a numerical aid, not the physical phase criterion.

### Conceptual separation: MF state selection versus RPA response

For each temperature \(T\) and transverse field \(B_x\), the bare single-ion
susceptibility \(\chi_0(B_x,\omega)\) must be evaluated in a self-consistent MF
state. The RPA response then applies the algebraic resummation

\[
\chi_{\mathrm{RPA}}(q,\omega)
=
\frac{\chi_0(\omega)}
     {1-J(q)\chi_0(\omega)}
\]

in the scalar projected longitudinal channel. RPA itself does not perform the
nonlinear order-parameter iteration and does not select a phase. It inherits the
phase through the eigenstates used to construct \(\chi_0\).

At zero longitudinal bias, the symmetric \(m_z=0\) state is an exact MF fixed
point on both sides of a continuous transition. Below the critical field it is
self-consistent but unstable. Consequently, the statement “the MF solver
converged” is insufficient: the code must also establish that the converged
state is the stable branch required by the RPA calculation.

### How the reference LiReF4 workflow selects a branch

The working reference implementation is located outside this repository under:

`Mean Field/LiReF4/`

Its relevant data flow is:

1. `LiReF4_MF_Yikai.m` initializes Ho with a nonzero positive \(J_z\) moment,
   scaled to 5.51.
2. The normal field scan is ordered from low to high field.
3. `LiIonsF4.m` resets the initial moment once before the field loop, not at
   every field. Each call to `remf.m` mutates `ion.mom`, so the converged state
   at field \(B_j\) seeds field \(B_{j+1}\).
4. `remf.m` iterates the nonlinear MF equation. The ordered solution is followed
   until its moment collapses into the paramagnetic solution at the continuous
   transition.
5. The resulting eigenenergies and eigenvectors are saved.
6. `MF_Chi.m` loads those saved MF eigenstates and constructs
   \(\chi_0(B,\omega)\).
7. `RPA.m` only evaluates
   \([I-\chi_0J(q)]^{-1}\chi_0\); it contains no phase selector or MF iteration.

Thus the reference code uses implicit branch following through a
symmetry-breaking seed and sequential continuation. It does not compare free
energies or explicitly reject an unstable self-consistent state.

### Existing projected implementation

`invz_projected/invz_chi_realaxis.m` currently supports
`opts.bare_rpa = true`. Its local `bare_rpa_point` performs the following
pointwise procedure:

1. Construct the paramagnetic single-ion state `siPm`.
2. Calculate the zero-frequency longitudinal susceptibility
   \(\chi_{0,cc}^{\mathrm{PM}}(0)\).
3. Calculate the paramagnetic RPA mass

   \[
   r_{\mathrm{PM}} = 1-J_{0,\mathrm{eff}}
                         \chi_{0,cc}^{\mathrm{PM}}(0).
   \]

4. If \(r_{\mathrm{PM}}\ge 0\), use the PM state.
5. If \(r_{\mathrm{PM}}<0\), solve the nonzero ordered MF root

   \[
   f(h_z)=h_z-J_{0,\mathrm{eff}}\langle J_z\rangle_{h_z}=0
   \]

   with a sign bracket and `fzero`.
6. Build \(\chi_0(\omega)\) from the selected state and set
   \(\Sigma(\omega)=0\) for the bare RPA leg.

This is deliberately independent of the phase and state selected by the
\(1/z\) solver. `opts.bare_rpa` rejects a caller-supplied `chi0cc_w`, preventing
accidental reuse of a susceptibility constructed from the wrong state.

### Interaction-ownership invariant

The same scalar uniform coupling must be used in all three operations:

1. the PM stability mass,
   \(1-J_{0,\mathrm{eff}}\chi_{0,cc}^{\mathrm{PM}}(0)\);
2. the ordered MF equation,
   \(h_z=J_{0,\mathrm{eff}}\langle J_z\rangle\); and
3. the intrinsic Gamma-point RPA denominator,
   \(1-J_{0,\mathrm{eff}}\chi_0(\omega)\).

This is the central invariant guarding against the original suspected failure:
a mismatch between the interaction that constructs the MF state and the
interaction in the RPA denominator.

`Jsel` may differ from `J0eff` only when evaluating a genuinely dispersive
finite-\(q\) response. `Jsel` selects the response momentum; it must not select
the MF phase. At the intrinsic Gamma point, assert to near machine precision
that `Jsel == J0eff`, allowing only a documented shape/demagnetization correction
applied after the intrinsic susceptibility has been constructed. `Jxx0` is the
transverse MF coupling and is a distinct physical parameter; it must remain
consistent across PM and ordered single-ion construction.

### Proposed hybrid algorithm

For a field sweep at fixed temperature, use the following algorithm.

1. Define a deterministic traversal. The preferred traversal for a
   ferromagnetic transverse-field transition is increasing \(B_x\), beginning
   in the ordered phase. If the caller supplies fields in another order, compute
   in a sorted work order and scatter results back to the requested output
   order. Do not silently alter the returned ordering.
2. At every field, construct the PM state and calculate \(r_{\mathrm{PM}}\).
   This test is mandatory even when a continuation state is available.
3. If \(r_{\mathrm{PM}}\) is clearly positive, select the PM state. Do not retain
   an ordered state merely because continuation converged to one.
4. If \(r_{\mathrm{PM}}\) is clearly negative, solve the nonzero ordered MF
   equation. Use the previous ordered field \(h_z\) as the first root/bracket
   hint. Verify a sign-changing bracket and the final residual. If the local
   continuation bracket fails, fall back to the existing global bracket search.
5. In a narrow critical tolerance band around \(r_{\mathrm{PM}}=0\), evaluate
   both the PM stability and the existence/residual of the nonzero root. Treat
   the transition point as the merger of the two roots; do not impose an
   arbitrary nonzero moment floor.
6. Build the dynamic \(\chi_0(\omega)\) from the state selected by steps 2--5,
   then calculate RPA for every requested \(q\). Reuse this selected state across
   momenta at the same \(T,B_x\), but never across different fields or between
   the \(1/z\) and bare-RPA theories.
7. Record diagnostics before advancing the continuation state.

MATLAB-style pseudocode:

```matlab
previous = [];
for kk = traversal_order(fields)
    Bx = fields(kk);

    pm = solve_paramagnetic_state(ion, T, Bx, Jxx0, hyp);
    chi0cc0 = static_longitudinal_chi0(pm, T);
    mass_pm = 1 - J0eff * chi0cc0;

    if mass_pm > mass_tol
        selected = pm;
        phase = 2;
    else
        hint = ordered_hint(previous);
        ordered = solve_nonzero_root_with_fallback( ...
            ion, T, Bx, J0eff, Jxx0, hyp, hint);

        if mass_pm < -mass_tol
            require_valid_ordered_root(ordered);
            selected = ordered;
            phase = 1;
        else
            [selected, phase] = resolve_root_merger(pm, ordered, mass_pm);
        end
    end

    require_interaction_consistency(J0eff, Jgamma);
    rpa_state{kk} = package_state(selected, phase, mass_pm);
    previous = rpa_state{kk};
end

% Dynamic response may be parallelized after the serial/static prepass.
parfor kk = 1:numel(fields)
    output{kk} = build_rpa_response(rpa_state{kk}, w, Jq);
end
```

The serial portion should contain only the static MF state construction. The
frequency- and momentum-dependent response remains embarrassingly parallel.

### Numerical decision policy at the root merger

The implementation must define the critical tolerance band; it must not leave
`resolve_root_merger` as an informal phase choice. Use the following policy:

1. Start with `mass_tol = 1e-10` for the dimensionless PM mass, but make it a
   named option. Increase it only when a measured error bound from the static
   susceptibility calculation requires doing so.
2. Outside the band, the sign of the PM mass is authoritative:
   `mass_pm < -mass_tol` requires the ordered root and
   `mass_pm > mass_tol` requires the PM state.
3. Inside the band, recompute the PM state and static susceptibility using the
   tightest supported MF tolerance. If the refined mass exits the band, use its
   sign.
4. If the refined mass remains inside the band, attempt the nonzero ordered
   root. Accept it only if its residual satisfies
   `abs(hz - J0eff*Jz) <= 1e-10 meV` and the solver reports convergence.
5. Treat an accepted root with
   `abs(hz) <= 1e-10 meV` or
   `abs(Jz) <= 1e-10` as merged with the PM root. Use the PM representation and
   record a `critical_merged` diagnostic rather than manufacturing a finite
   ordered moment.
6. If a finite accepted root exists and the best resolved mass is nonpositive,
   use the ordered state. If the best resolved mass is positive, use the PM
   state. If neither condition can be certified, mark the point invalid and
   request field-grid refinement.

The numerical values above match the scale of the existing ordered-root
residual check and are initial defaults, not universal physical constants. They
should be scaled if the Hamiltonian energy unit changes.

For the Gamma-point coupling invariant, use a symmetric floating-point check,
for example

```matlab
coupling_tol = 64*eps(max([1, abs(J0eff), abs(Jgamma)]));
assert(abs(Jgamma - J0eff) <= coupling_tol);
```

The assertion must be applied to the intrinsic coupling before any separately
labelled sample-shape correction.

### Ownership and minimal file growth

There must be one owner of the bare-RPA MF equations. Do not copy the PM mass,
ordered residual, or coupling-selection logic into a spectra driver.

The clean implementation is one helper boundary, for example
`invz_bare_rpa_state`, supporting both scalar state construction and an optional
field-sweep continuation hint. `invz_chi_realaxis` consumes the returned state;
`invz_spectra_map` coordinates traversal and parallel dynamic evaluation. This
would add at most one focused helper file.

If avoiding even that file is a strict requirement, keep the existing scalar
helpers local to `invz_chi_realaxis.m` and retain the current independent
pointwise algorithm. Do not implement continuation by duplicating those
equations as local helpers in `invz_spectra_map.m`; duplicated ownership is a
larger long-term risk than the absence of continuation.

### Required diagnostics

Each field point should retain:

- `phase_rpa`: ordered, PM, or invalid/masked;
- `rpa_mass_pm`;
- selected \(h_z\) and \(\langle J_z\rangle\);
- ordered MF residual
  \(h_z-J_{0,\mathrm{eff}}\langle J_z\rangle\);
- MF convergence status and iteration count where available;
- the exact `J0eff`, `Jxx0`, and Gamma-point `Jsel` used;
- whether the ordered solve used continuation or the global fallback bracket;
- the minimum magnitude or relevant eigenvalue of the static RPA denominator.

A concrete returned state record should contain at least `si`, `tl`,
`is_ordered`, `phase_rpa`, `rpa_mass_pm`, `hz`, `Jz`, `mf_residual`,
`mf_converged`, `selection_source`, and a `couplings` structure containing
`J0eff`, `Jxx0`, and `Jgamma`. The dynamic response must consume this record
without rebuilding or changing its MF state.

Failures to converge, bracket a root, or satisfy the coupling invariant must be
reported as invalid points. They must not be converted silently to the PM phase.

### Validation gates

Do not adopt the hybrid implementation solely because it produces a smoother
plot. Require all of the following:

1. **Current-baseline preservation:** away from the QCP, hybrid and current
   pointwise selectors agree on phase, moment, \(\chi_0\), and RPA poles within
   numerical tolerance.
2. **Continuous QCP behavior:** the ordered \(m_z\) approaches zero, the PM mass
   approaches zero, and the RPA soft pole reaches zero at a common field within
   field-grid resolution.
3. **No missing interval:** the soft pole does not disappear before the critical
   point and reappear after it.
4. **Grid refinement:** refining the field grid changes the inferred critical
   field and minimum pole by amounts consistent with discretization error.
5. **Order independence:** ascending, descending, and randomly ordered input
   arrays produce the same scattered-back physical result. Continuation may
   change runtime, not the selected phase.
6. **Cold-start comparison:** selected checkpoints agree with independent
   pointwise solves that do not use the preceding field.
7. **Interaction audit:** at Gamma, the MF/stability coupling and intrinsic RPA
   denominator coupling agree to near machine precision.
8. **Adversarial failure test:** deliberately corrupting the continuation hint
   must trigger fallback or rejection, not alter the final physical branch.

For a first-order transition or coexistence region, PM stability alone does not
identify the globally stable phase. In that case this design must be extended
with a consistently normalized MF free-energy comparison; continuation plus a
local stability test is not sufficient.
