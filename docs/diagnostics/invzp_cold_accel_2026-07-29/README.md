# Safeguarded cold-start acceleration — 2026-07-29

## Scope and frozen safeguards

This packet tests the narrow solver claim identified in
`invzp_convg_fix.md`: the 0.10 K, 4.400 T edge mask is a locally
contractive, signed oscillatory outer mode that can be removed without a
physical-field warm seed. It does not propose a low-field branch selector.

The fixture is the legacy brute-force `16^3`, `dpRng=30` coupling multiset,
the resummed static medium, and `transverse_mf='legacy_x'`. In the current
worktree its exact digest is
`499922e6c9f8c44d51b5db06486aac345b6226d1f8096713d20916ca78612cb5`;
cached and freshly regenerated vectors are bit-identical. Older handoff
material registers
`ddb9532d11326458554b93b2ce09c80a3299cca9baa69202dc311f138b4fae17`
for the nominally same fixture. That historical byte digest was not
reproduced and remains a provenance conflict; the current vector does
reproduce the documented extrema and all pre-change solver anchors used
here. The
accelerated 4.400 T solve uses `mix_outer=0.70`, `max_outer=200`; its
references are ordinary cold `0.70/1000` and predictor-only QCP-down warm
`0.50/1000`.

`cold_acceleration='signed_aitken1'` is default off, resummed-only, and
applies only to a genuinely cold attempt (including the existing fresh
cold retry after a failed warm attempt). A proposal must independently
pass all of these fixed gates:

- four successive full `Sigma`-vector increments and three signed ratios;
- at least eight ordinary iterations;
- pooled `lambda` in `[-0.99,-0.50]`, ratio spread at most `0.02`, and
  relative scalar-mode fit error at most `0.10`;
- one unchanged coupling interval, finite nonzero pole distance, and
  closed inner static solve over the history;
- fresh unmixed full-vector residuals at both the ordinary next iterate and
  proposal, with the proposal strictly smaller than both the ordinary
  candidate and the current residual.

Every accepted or rejected proposal clears the history. A later proposal
must re-earn every gate. Ordinary iteration resumes between proposals, and
the unchanged A--D/all-node audit remains the only final acceptance rule.

## Results

At 4.400 T, a first one-shot implementation accepted a proposal at outer
143 (`lambda=-0.908610599`) and reduced the fresh residual from
`1.17e-3` to `9.20e-5`, but did not close within 200 iterations. This
falsified the extra one-proposal restriction; the documented
restarted-Anderson-1 alternative requires restarted proposals, not a
one-shot extrapolation.

With restart semantics, three independently gated proposals were accepted:

| outer | signed `lambda` | ratio spread | mode-fit error | ordinary residual | proposal residual | interval rank |
|---:|---:|---:|---:|---:|---:|---:|
| 143 | -0.908610599 | 0.0194 | 0.0148 | 1.17e-3 | 9.20e-5 | 16384 |
| 147 | -0.911828939 | 0.00866 | 0.00840 | 6.47e-5 | 4.20e-7 | 16384 |
| 152 | -0.921293084 | 0.0124 | 0.0117 | 2.71e-7 | 8.45e-10 | 16384 |

Ordinary iteration then closed at outer 153, and the complete HMF solve
was accepted and stable. Its full final state is bit-identical to the
pre-change ordinary `0.70/1000` cold reference:

```text
max|Delta Sigma|  = 0
max|Delta K|      = 0
max|Delta lambda| = 0
Delta hstar       = 0
Delta m           = 0
```

The predictor-only `0.50/1000` warm state has the same `hstar` and moment;
relative to the accelerated state its maximum absolute differences are
`2.07011224e-10` in `Sigma`, `2.69215769e-12 meV` in `K`, and
`1.87829691e-11` in `lambda`. All three states pass the unchanged final
audit.

The 4.300 T trace is a negative mechanism control. Its cold `h=0`
predictor accepts normally in 13 iterations and never activates the
accelerator. The sole failed profile node is instead at
`h=0.0028118554912421501 meV`; after 1000 iterations it has
Block-A residual `0.00833` despite static residual `7.57e-11`. Its late
iterations switch between interval ranks 16384 and 16376, with pooled
factors/mode-fit errors that fail the registered scalar-mode and
same-interval gates. No proposal is attempted. Thus 4.300 T is not the
4.400 T budget-limited predictor mechanism.

The user additionally reported that the visible sliver near 3.825 T
converged under neither cold nor QCP-down warm seeding. The cold
acceleration trace agrees and makes the failure sharper: 17 of 34 evaluated
nodes fail at the 1000-step cap. The predictor itself ends with Block-A
residual about `1.05` on interval rank 0; numerous low-`h` nodes occupy
different interior ranks. No predictor proposal qualifies. One later
fresh-cold retry at `h=0.0057150624450985578 meV` accepts a proposal at
outer 9 but still fails the unchanged node audit after 1000 iterations
(Block A `1.82e-5`). This is broad path nonconvergence, not a recoverable
single signed edge mode.

## Verdict

Safeguarded cold acceleration is a successful local remedy at 4.400 T and
removes warm-seed direction/basin dependence for that column: it reaches
the exact long-cold state within the original 200-iteration budget.
It correctly declines to extrapolate the non-scalar, interval-switching
4.300 T failure and does not repair the multi-node 3.825 T failure.
Nothing here selects among low-field thermodynamic branches or authorizes a
default flip. The historical/current coupling-digest conflict must also be
resolved before treating this packet as a reproduction of the older frozen
Gate-0 fixture rather than of the exact current QCP-edge fixture.

`invzp_cold_acceleration_gate.m` reproduces the three-field discriminator
and saves a compact result if given a path. `cold_acceleration_results.mat`
contains the retained compact output; full temporary traces were not
retained because the 3.825 T trace alone is about 256 MB.
