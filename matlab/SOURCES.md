# Where every number and equation came from

Read this before citing anything. It is split into what is certain, what you
should verify, and what is deliberately not included.

---

## 1. The four sources

**[1] Nguyen, Ogburn, Gilbert, Kibler, Brown, Deal (1979).** *Simulator Study
of Stall/Post-Stall Characteristics of a Fighter Airplane With Relaxed
Longitudinal Static Stability.* NASA Technical Paper 1538.

The primary source for everything: the wind-tunnel aerodynamic tables, the
engine thrust tables, the mass and geometry data, the control surface limits,
and the deep-stall pitching moment increment. You already cite this as
`nguyen1979`. Every other source below is a restatement of it.

**[2] Stevens, B.L. and Lewis, F.L.** *Aircraft Control and Simulation.*
Wiley. 2nd edition 2003, or 3rd edition 2016 with Johnson.

Appendix A and Chapter 3 reproduce the F-16 model in a form that is easier to
work from than the NASA report: reference geometry, moments of inertia, the
engine model including the throttle gearing and the three thrust tables, and a
worked trim example. This is the standard citation in the flight dynamics
literature for "the F-16 model" and the one your examiner will expect to see.

**[3] van Oort, E. and Sonneveldt, L. (2006).** *F16 aerodata tables from NASA
report 1538.* Delft University of Technology, Control and Simulation Division.

This is the `F16aerodata.m` file itself. Cite it for the data, not the model.
You already cite it as `vanoort2006`.

**[4] Sonneveldt, L. (2006).** *Nonlinear F-16 Model Description.* Delft
University of Technology, Control and Simulation Division.

The companion document to [3]. This is the source of the coefficient buildup
equations implemented in `f16_aero.m`: how the leading-edge-flap tables are
blended, how the aileron and rudder increments are scaled, how the damping
terms are non-dimensionalised, and how the c.g. transfer terms enter Cm and Cn.
If you cite one document for the equations, cite this one. The same equations
appear in Russell, R.S. (2003), *Non-linear F-16 Simulation using Simulink and
Matlab*, University of Minnesota, if you want a second reference.

---

## 2. Values I put in `f16_constants.m`

These are the standard model values from [1], as tabulated in [2]. The imperial
original is given so you can check the conversion.

| Quantity | Value used (SI) | Original |
|---|---|---|
| Wing reference area S | 27.87 m² | 300 ft² |
| Reference span b | 9.144 m | 30 ft |
| Mean aerodynamic chord c̄ | 3.450336 m | 11.32 ft |
| Mass | 9295.44 kg | 636.94 slug (20 500 lb) |
| Ixx | 12874.8 kg·m² | 9 496 slug·ft² |
| Iyy | 75673.6 kg·m² | 55 814 slug·ft² |
| Izz | 85552.1 kg·m² | 63 100 slug·ft² |
| Ixz | 1331.4 kg·m² | 982 slug·ft² |
| Engine angular momentum | 216.9 kg·m²/s | 160 slug·ft²/s |
| Reference c.g. | 0.35 c̄ | 0.35 c̄ |

Conversion factors: 1 ft = 0.3048 m, 1 ft² = 0.09290304 m², 1 slug = 14.5939 kg,
1 slug·ft² = 1.35582 kg·m², 1 lbf = 4.44822 N.

Two corrections to your original `f16_parameters.m`:

- `empty_weight = 9.207` and `Max_Togw = 21.772` were labelled kg but are
  tonnes. Corrected to 9207 kg and 21772 kg.
- `b = 9.449` m is the physical span over the wingtip launchers. The
  *aerodynamic reference* span is 9.144 m (30 ft). The coefficients in
  `F16aerodata.m` are non-dimensionalised with the reference value, so using
  9.449 m would introduce a 3% error in every rolling and yawing moment and in
  every roll and yaw damping term. Both values are now in `f16_constants`,
  clearly separated.

---

## 3. What you must verify yourself

**The three engine thrust tables in `f16_engine.m`.** These are the one part of
the model I could not check against a primary document. They are the idle,
military and maximum tables from [1] / [2], indexed by Mach number (rows, 0 to
1.0) and altitude (columns, 0 to 50 000 ft), in pounds force. Open your copy of
[2] Appendix A and check every entry against the table before you cite any
throttle setting or thrust-limited result in the report.

The good news: trim angle of attack, trim stabilator, static stability, the
departure criterion and the equilibrium map are all essentially insensitive to
thrust. Only the trim throttle value depends on these numbers. So an error here
does not invalidate the main analysis, but it would embarrass you in a viva.

**The throttle gearing.** `tgear(δt) = 64.94·δt` below 0.77 and
`217.38·δt − 117.38` above it. The break at 0.77 is the afterburner detent.
Also from [2], also worth checking.

**Page and figure numbers.** I have given document titles and years, which are
reliable. I have not given page, table or equation numbers, because I would be
guessing at them. Fill those in from your own copies.

---

## 4. What is deliberately not in the model

**The stabilator effectiveness factor η(δh).** The full NASA formulation
multiplies the tabulated Cm by a tail-effectiveness factor that depends on
stabilator deflection. That table is not in `F16aerodata.m`, so `f16_aero.m`
uses η = 1. Mention this as a modelling assumption in your methodology section.
If you find the table in [1], it drops straight into the Cm line.

**Speedbrake.** The ΔCX contribution of the speedbrake is not in the table set
either. Not needed for your analysis.

**Leading-edge-flap scheduling.** The real F-16 schedules the LEF as a function
of angle of attack and Mach number. `f16_aero` takes δlef as an input so you can
either hold it fixed, which is what the scripts do, or add the schedule later.
Holding it at 0 is the honest choice for a stability study, because it isolates
the airframe from the control system. Say so in the report.

**Unsteady aerodynamics.** Your Section 3 argues correctly that a static lookup
table becomes inaccurate at high alpha. Nothing in this code addresses that: it
is a static model. If you want to close that gap, the Goman–Khrabrov model adds
one internal state for the separation point, with a relaxation time constant and
an alpha-dot lag, and blends attached and separated coefficient values with it.
See Goman, M. and Khrabrov, A. (1994), "State-Space Representation of
Aerodynamic Characteristics of an Aircraft at High Angles of Attack", *Journal
of Aircraft*, 31(5), 1109–1115. That is the paper your Section 3 is describing.

---

## 5. BibTeX

```bibtex
@techreport{nguyen1979,
  author      = {Nguyen, Luat T. and Ogburn, Marilyn E. and Gilbert, William P.
                 and Kibler, Kemper S. and Brown, Phillip W. and Deal, Perry L.},
  title       = {Simulator Study of Stall/Post-Stall Characteristics of a
                 Fighter Airplane With Relaxed Longitudinal Static Stability},
  institution = {NASA Langley Research Center},
  number      = {NASA TP-1538},
  year        = {1979}
}

@book{stevenslewis,
  author    = {Stevens, Brian L. and Lewis, Frank L.},
  title     = {Aircraft Control and Simulation},
  publisher = {John Wiley \& Sons},
  edition   = {2},
  year      = {2003}
}

@misc{vanoort2006,
  author = {van Oort, Ewoud and Sonneveldt, Lars},
  title  = {{F-16} Aerodata Tables from {NASA} Report 1538},
  howpublished = {Control and Simulation Division, Delft University of Technology},
  year   = {2006}
}

@techreport{sonneveldt2006,
  author      = {Sonneveldt, Lars},
  title       = {Nonlinear {F-16} Model Description},
  institution = {Control and Simulation Division, Delft University of Technology},
  year        = {2006}
}

@article{goman1994,
  author  = {Goman, Mikhail and Khrabrov, Alexander},
  title   = {State-Space Representation of Aerodynamic Characteristics of an
             Aircraft at High Angles of Attack},
  journal = {Journal of Aircraft},
  volume  = {31},
  number  = {5},
  pages   = {1109--1115},
  year    = {1994}
}
```

---

## 6. Sign conventions, so the report says the right thing

Checked directly against the tables. At α = 0, β = 0:

| δh | −25° | −10° | 0° | +10° | +25° |
|---|---|---|---|---|---|
| Cm | +0.141 | +0.043 | −0.060 | −0.161 | −0.253 |

So **positive δh is trailing edge down and produces a nose-down moment**, which
is the conventional sign. Full nose-up is δh = −25°, full nose-down is +25°.
Your original script's printout did not state this and it is easy to get
backwards when writing the results section.
