# F-16 high angle-of-attack model: file guide

Put all of these in the same folder as your existing `F16aerodata.m`.

| File | What it is |
|---|---|
| `f16_constants.m` | Geometry, mass, inertia, control limits. Every value sourced in `SOURCES.md`. |
| `f16_aero.m` | The full nonlinear coefficient buildup. Uses all ~45 tables. This is the core. |
| `f16_engine.m` | Thrust tables, throttle gearing, ISA atmosphere. |
| `f16_trim.m` | Real trim: solves force *and* moment balance for throttle, δh, α. |
| `run_pitch_equilibria.m` | Corrected replacement for your `Trim.m`. |
| `plot_coefficients.m` | Corrected replacement for your `f16_parameters.m`. |
| `SOURCES.md` | Citations, BibTeX, and what you still need to verify. |

Your original `Trim.m` and `f16_parameters.m` are superseded but keep them:
the difference between the old and new results is itself worth a paragraph in
the methodology section.

## Run this first

```matlab
[thtl, dh, alpha, res] = f16_trim(152.4, 0, 0, 0.35)
```

Expected: `alpha ≈ 2.0` deg, `dh ≈ -0.2` deg, `thtl ≈ 0.13`, residual at
machine precision. That is level flight at 500 ft/s at sea level with the c.g.
at the reference position. If you get this, the aerodynamic buildup, the
constants and the sign conventions are all consistent and you can trust
everything downstream. If you do not, stop and debug before running anything
else.

Then:

```matlab
run_pitch_equilibria     % equilibrium map and the deep-stall branch
plot_coefficients        % coefficient figures and the departure criterion
```

## What `run_pitch_equilibria` shows

It sweeps the stabilator from full nose-up to full nose-down at three c.g.
positions and finds every angle of attack at which the pitching moment
vanishes, marking each as statically stable or unstable.

The result to look for: at x_cg = 0.35 c̄ with the stabilator at full nose-down
(δh = +25°), the model still has a statically stable equilibrium at high angle
of attack. That is the locked-in deep stall. Full nose-down control is applied
and the aircraft sits there anyway. It is the single most useful figure you can
put in the results section, because it is a purely nonlinear phenomenon with no
counterpart in linear theory, and it is created by exactly the flow physics your
Section 2 describes.

Two things to notice while writing it up:

- Move the c.g. forward to 0.30 c̄ and the picture changes qualitatively, not
  just quantitatively. The number of equilibria changes. That is a bifurcation
  in the c.g. parameter, and it connects directly to the "relaxed longitudinal
  static stability" in the title of NASA TP-1538.
- At x_cg = 0.35 the low-α equilibrium is statically *unstable*. That is not a
  bug. The F-16 is deliberately unstable in pitch at subsonic speeds and relies
  on the flight control system. Worth a sentence, because a reader who does not
  know that will assume you made an error.

## Known limitations of what is here

- Static aerodynamics only. No unsteady or hysteresis effects. See `SOURCES.md`
  section 4.
- `run_pitch_equilibria` maps pitch-moment equilibrium at a fixed flight
  condition, not full six-degree-of-freedom equilibrium. It is robust and it
  answers the deep-stall question. A genuine continuation with airspeed and
  flight path angle free is the next step, and it is fiddly: the equilibrium set
  has several branches and a naive solver hops between them. When you build it,
  step the continuation parameter in small increments, seed each solve with the
  previous solution, and reject any step where α jumps by more than a degree or
  two.
- The engine thrust tables need checking. See `SOURCES.md` section 3.

## Suggested next three steps

1. Wrap `f16_aero` and `f16_engine` in a six-degree-of-freedom
   `xdot = f16_dynamics(x,u)` and verify that a trimmed state stays trimmed.
2. Numerically linearise about trim points across the α range, take
   eigenvalues, and track short period, dutch roll and spiral. This gives the
   quantitative stability answer your research question asks for.
3. Then the continuation and bifurcation study, and time simulations of
   departure and of an attempted deep-stall recovery.
