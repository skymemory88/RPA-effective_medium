# Display-only energy unit knob (`eUnit`) for `invz_run_spectra`

Date: 2026-07-08. Branch: `invz-1z-lihof4`. Status: design approved by user
(conversation of 2026-07-08).

## Problem

`invz_run_spectra.m` computes chi''_cc(omega) on a frequency grid `w` in meV
(required by the solver layer) and plots it with the same meV values on the
axis. The user wants to view the same spectra in GHz without touching the
computation. The driver has two mutually exclusive plot branches selected by
field count (`numel(fields) <= sliceMax`): a line-slice overlay drawn inline,
and a 2D field-vs-frequency colormap rendered via `invz_plot_spectra_map.m`.
A unit knob that only covers one branch would be a silent no-op under the
script's own current settings (`fields = linspace(0,9,201)`, `sliceMax = 6`
puts it in the colormap branch).

## Solution overview

Add a knob `eUnit = 'GHz' | 'meV'` near the top of `invz_run_spectra.m`.
Conversion is display-only and lives entirely in the driver, reusing the
existing (currently unused) `invz_const().Gh2mV` constant (GHz -> meV). The
grid `w` fed into `invz_spectra_map` for the actual solve is never touched.
`invz_plot_spectra_map.m` gains one optional trailing parameter used only to
build the y-axis label string — it does no conversion arithmetic itself, so
callers that already pass numeric data in whatever unit they like (e.g. the
docstring's ad hoc REPL replot example) keep working unchanged.

## Components

### 1. Driver changes `invz/invz_run_spectra.m`

Knob with the other user-facing settings near the top:

```matlab
eUnit = 'GHz';   % 'meV' or 'GHz' -- plotting only; computation always runs in meV
```

Immediately after the knobs, validate and derive a scale factor / label
(fails fast, before the sweep runs):

```matlab
C = invz_const();
switch eUnit
    case 'meV', eScale = 1;         eLabel = '\omega (meV)';
    case 'GHz', eScale = 1/C.Gh2mV; eLabel = '\omega (GHz)';
    otherwise, error('invz_run_spectra:eUnit', 'eUnit must be ''meV'' or ''GHz''.');
end
```

- Line-slice branch: plot `w*eScale` in place of `w`; `xlabel(eLabel)`.
- Colormap branch: build a display-only shallow copy `Splot = S; Splot.w =
  S.w*eScale;` and call `invz_plot_spectra_map(ax, Splot, Splot.chiz, ttl,
  eUnit)` (new 5th argument).
- `w` passed into `invz_spectra_map(ion, T, fields, w, ...)` is unchanged --
  the solve always runs on the meV grid.

### 2. Plot helper changes `invz/invz_plot_spectra_map.m`

- New optional 5th parameter `eUnit`, default `'meV'`.
- `ylabel(ax, '\omega (meV)')` becomes `ylabel(ax, sprintf('\omega (%s)',
  eUnit))`.
- No numeric conversion added -- the function keeps its existing contract of
  plotting exactly the data it is given (`S.w`, already display-scaled by the
  caller).

## Error handling

- Invalid `eUnit` value (typo, wrong case, etc.): explicit `error(...)`
  immediately after the knob, before the potentially 15-25 minute sweep runs.
- No other new failure modes -- the change touches display code only.

## Testing

Display-only labeling knob with no physics/numerics to unit test, and
`invz_plot_spectra_map.m` is not under automated test today (confirmed: no
existing test file references it). Verify manually by running the line-slice
branch (small `fields` list, fast) under both `eUnit` values and confirming
axis numbers scale by `1/C.Gh2mV` (~241.8x) with the label updating; no new
automated test is added.

## Success criteria

- Fast test suite still passes unchanged (this change touches no file it
  covers).
- Running the driver with `eUnit = 'meV'` reproduces today's plots exactly
  (byte-for-byte same axis values/labels).
- Running the driver with `eUnit = 'GHz'` shows the same spectral shape with
  frequency values scaled by `1/C.Gh2mV` and axis labels reading `(GHz)`, in
  both the line-slice and colormap branches.
- Invalid `eUnit` raises a clear error before the sweep starts.

## Out of scope

- Converting any non-plotted quantity (e.g. `eta`'s doc comment, printed
  `Sigma0` diagnostics) -- computation and console output stay in meV.
- Changing `invz_spectra_map.m` or any solver-layer file.
- Adding automated tests for MATLAB figure labeling.
