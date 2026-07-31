## 0. The reframe that makes most of this optional

In the bare doublet reduction, using your production splittings:

$$\frac{\hmf^*}{\hmf_c}=\frac{\sqrt{\rho^2-1}}{\sqrt{\rho^{2/3}-1}},\qquad \rho=b_c/b\;\Longrightarrow\;\frac{\hmf^*}{\hmf_c}\ \ge\ \sqrt3\quad\text{for all }\rho ,$$

with the bound approached only as \(B_x\to B_c\). At 1 T: fold 0.00548, root 0.03511 — a factor 6.4. At 4 T: your \(h_c=0.0064876\), root 0.011451 — a factor 1.77. **The root is never in the forbidden region, and never closer than \(\sqrt3\) to it.** So the fold, the two sheets at \(h=0.006\), and the mass edge all live in territory the answer never needs to visit — *provided the anchor comes from above rather than below*. That is the first suggestion, and it makes suggestions 5–7 diagnostics rather than blockers.

## 1. Anchor high, transport down

Above the root the theory is weakly coupled: \(\Jc_0\tilde\chi_0=1-\text{mass}=(b/b_c)^2\), which is **0.016 at the 1 T root**. So:

1. Pick \(h_{\rm anchor}\) where \(\Jc_0\tilde\chi_0\le10^{-2}\). At 1 T that is \(\hmf=0.0411\) meV — only \(1.17\times h^*\).
2. Compute \(\delta H(h_{\rm anchor})\) **perturbatively**, at second order in \(\Jc\), where \(K=-\langle\Jc^2\rangle_qG_0\) with no self-consistency and no sheet ambiguity. Relative error \(O((\Jc\chi)^2)\sim10^{-4}\) — an order below the truncation error.
3. Transport *downward* with (45)'s ODE along the certified branch to \(h^*\). Finite interval, entirely above both obstructions.

The anchor can be computed two independent ways — the perturbative limit of the \(S_3\) formula, and the perturbative limit of the \(\lambda\)-integral — which **must agree**, because at second order both reduce to the same textbook calculation. That agreement is the normalization check the review correctly demands of \(S_3\), obtained without settling the full ordered algebra.

Regime map: at 4 T the same criterion puts the anchor at \(\hmf=0.155\) meV, i.e. \(\hmf J=0.85\) meV against the 11 K first excited crystal-field level at 0.95 meV — CF admixture is no longer negligible, so the high anchor degrades exactly where \(B_x\to B_c\). But that is where \(h_c\to0\) and the framework's own continuous-transition closure applies. Two anchors, complementary domains, overlapping around 2–3 T where they can be cross-validated.

## 2. Sheet selection: \(\lambda\)-homotopy, using the arclength code you already wrote

At \(\lambda=0\) the solution is unique with \(\Sigma=0\). The \(1/z\) solution is the one analytically connected to it; sheets that exist only after the Dyson resummation (30) have no order-by-order counterpart and are roots of the resummed algebra, not solutions of the theory. So: continue in \(\lambda\) at fixed \(h\), with **pseudo-arclength in \(\lambda\) instead of \(h\)** — same bordered system, different continuation parameter. Two sharp outcomes:

- lands on the contracting sheet (\(\lambda_{\rm dom}=0.4173\)) ⟹ selection settled, and the by-product is \(\delta\Phi_c(m)\);
- **folds before reaching \(\lambda=1\)** ⟹ the target state is not a \(1/z\) solution at all, which is a far stronger statement than "Picard did not converge".

Run it at \(h\) just below the 1 T fold (0.0054) as well: that decides existence there constructively rather than by budgeted search.

## 3. A \(\lambda\to0\) unit test — the first *external* check in the pipeline

Everything you currently validate is a self-consistency residual. Perturbation theory supplies an independent target: with \(\Jc\to\lambda\Jc\) and \(\langle\Jc_q\rangle=\Jc(ii)=0\),

$$K(\wn)=-\lambda^2\langle\Jc_q^2\rangle_q\,G_0(\wn)+O(\lambda^3)\;\Longrightarrow\;\Sigma=c_2\lambda^2+O(\lambda^3),$$

with \(c_2\) fully explicit from \((M,m,\Delta,n_{01})\) and the elementary moments \(\frac1\beta\sum g=1\), \(\frac1\beta\sum g^2=\frac12[g(0)+\beta(1-n_{01}^2)]\), \(\frac1\beta\sum g^3=3/2\Delta^2\), \(\frac1\beta\sum g^4=5/2\Delta^3\). Note \(\xi\to1\) and \(G(0)\to G_0(0)\) in this limit, so the test covers the elastic sector and the \(\tanh\) resummation too. Plot \(\Sigma(\lambda)/\lambda^2\) and check it approaches \(c_2\). Any failure localizes a bug in (39)–(42) that no residual test can see. Same sweep also measures how much work the \(\tanh\) is doing: track its argument \(m^2n_{01}^2\beta K(0)-M^2\beta\lambda_1\) along \(\lambda\); if \(|{\rm arg}|\gg1\) at \(\lambda=1\), the elastic sector is being carried entirely by a resummation Jensen introduced for boundedness, not accuracy.

## 4. Reproduce Jensen's HoF₃ ordered phase end to end

This is the decisive discriminator between an implementation defect and a regime problem, and I don't see it in the diagnosis. HoF₃ is a genuine singlet–singlet system — no hybrid, the framework's native scalar case — and Jensen publishes ordered-phase numbers to hit: \(T_c=0.53\) K, \(1+\Sigma(0)=1.140\) at \(T\simeq2\) K and 1.098 at 4.2 K, 1.128 just above \(T_c\), \(\delta U/N=-9.06\;\mu\)eV at low \(T\), the moment deficit improving 16%→9% at \(T=0\), his Fig. 7 \(\langle J_x\rangle(T)\) and \(\langle I_x\rangle(T)\) below \(T_c\), and (2.34) satisfied to 2–3%. Crucially, HoF₃ has \(m^2/M^2\lesssim2\) at saturation, against your measured \(1.37\times10^3\): **the ordered formulae have never been exercised anywhere near your regime.** If Fig. 7 reproduces, the algebra is sound and the LiHoF₄ failure is a regime/anchor issue; if it does not, everything above is premature.

## 5. What creates the 1 T fold — a testable mechanism

Under a uniform shift of \(K\), the ordered self-energy responds with gain

$$\frac{\partial\Sigma(0)}{\partial K}\bigg|_{\rm unif}=g(0)\big[M^2-3m^2\big]
=\frac{J^2n_{01}\,[\,b^2-3h_L^2\,]}{(b^2+h_L^2)^{3/2}},$$

(the \(M^2\) cancels in the \(-2m^2\gamma(0)g(0)/M^2\) term of (39); \(\alpha\) is insensitive by the sum-rule structure). At your 1 T fold this is \(-1092\) meV\(^{-1}\). It peaks at \(h_L=\sqrt3\,b\), i.e. \(\hmf=0.007744\) meV at 1 T with magnitude 1424 meV\(^{-1}\) — versus 190 meV\(^{-1}\) at 4 T, peaking at \(\hmf=0.058\) meV, far above the whole profile. **That is exactly the observed dichotomy: a gain-driven fold at 1 T, a mass-edge termination at 4 T.** (Incidentally, your accepted node 0.0077454658 sits at \(h_L/b=1.73225\) against \(\sqrt3=1.732051\); probably coincidence, but worth one glance.)

Three cheap consequences:
- Measure the loop gain \(|\partial\Sigma/\partial K|\cdot|\partial K/\partial\Sigma|\) from the Jacobians you already form; it should cross 1 at the fold.
- **Predicted second fold.** The gain rises through the critical value at \(h_L/b\simeq1.22\) and falls back through it at \(h_L/b\simeq2.65\), so the two sheets should rejoin near \(\hmf\simeq0.0118\) meV at 1 T. Your second solution at \(h=0.006\) (\(h_L/b=1.34\)) lies inside that window. Continue the \(\lambda_{\rm dom}=2.64\) sheet upward and look for the closure.
- **Artifact test.** Replace the \(\xi\) \(\tanh\) by its linear term, and vary the static-frequency replacement \(g(\omega_{n'}\pm\omega_n)\to g(\omega_{n'})\). If the fold moves substantially, it is generated by the least-controlled part of the ordered algebra, amplified by \(m^2/M^2\), rather than by physics.

## 6. Why your triangle-inequality bound failed — and how to fix it

At \(h=0\) the elastic term vanishes, so \(\widetilde G_0=G_0/(1+\Sigma)\) exactly and admissibility is *identically* \(1+\Sigma(0)>\Jc_0\chi_0(0)\). Using \(K\le\max_q\Jc_q\) inside the positive-mass domain, \(\lambda_1\le K_{\max}\), \(\lambda_2\le K_{\max}\cdot\frac12[g(0)+\beta(1-n_{01}^2)]\), the \(\beta(1-n_{01}^2)\) terms cancel and

$$\Sigma(0)\big|_{h=0}\;\le\;\frac{M^2g(0)}{n_{01}^2}\,K_{\max}\;=\;\frac{\Jc_{\max}\,\chi_0(0)}{n_{01}^2}.$$

Since \(\Jc_{\max}=\Jc_0\) for a \(\mathbf q\to0\) ferromagnet, the bound is \(\Jc_0\chi_0/n_{01}^2\) while the requirement is \(\Jc_0\chi_0-1\). **The bound fails by exactly one unit of \(\Sigma\), universally** — which is why it failed at every field from 0.5 to 2.2 T, and why no sharpening of the \(q\)-side can ever close it. The slack has to come from the frequency side: \(K(\omega_n)\) is a positive-weight average of \(\Jc_q\) with tilt \(\propto|G_0(\omega_n)|\), hence maximal at \(n=0\) and decaying like \(g(\omega_n)\), so a single \(x\) cannot push \(d_q\) near zero at every Matsubara frequency. The \(n=0\) share of the moment weights is only \(g(0)/\beta=0.347\) for \(\lambda_1\); imposing that structure gives roughly \(\Sigma\le640\,K_{\max}\approx4.1\) against a required 6.9 at 1 T — an exclusion. Making that rigorous is a bounded, one-dimensional piece of analysis in \(x\), and it is exactly the "explicit boxes rather than budgeted starts" your item 2 asks for. Two free necessary conditions to run first: \(\chi_{\rm rest}(h=0)<1/J_{\rm sup}=155.72\) meV\(^{-1}\) (else no \(\Sigma\) can fix the hybrid at all), and the solver-free ledger \(\Sigma_{\rm req}(0)=\Jc_0\chi_0(0)-1\) tabulated across \(B_x\), which also locates \(B_c(0.1\,{\rm K})\) as the field where it meets the actual \(\Sigma(0)\).

## 7. The hybrid tail, stated exactly

With \(w\) the dominant-sector share (your "weak-transition remainder" is \(1-w\)):

$$r-1=\frac{w\,\Sigma}{1+(1-w)\Sigma}.$$

So the fitted \(|r-1|\) exponents are not \(\Sigma\)'s: with \(1-w\) falling 0.216→0.092 at 1 T, \(\Sigma\sim h^{-0.65}\) rather than \(h^{-0.59}\). More usefully, \(w\) is bounded away from zero, so **the mixing cannot rescue integrability** — the question is entirely about \(\Sigma\). Combined with the retraction, that leaves the additive decomposition (\(\alpha\), \(-\alpha_m\), \([\gamma-\frac{2m^2}{M^2}\gamma_0]g\), and the elastic \(m^2\xi h(0)\)) as the right diagnostic, with the \(\frac{2m^2}{M^2}\gamma(0)g(0)\) term the first suspect given §5.

---

If I were sequencing: §4 and §3 first as gates (a day, and they can invalidate everything else), then §1 for the deliverable, §2 for selection, §5–6 as understanding. I can reissue the HTML with all corrections plus these — say the word.