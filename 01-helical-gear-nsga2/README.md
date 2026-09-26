# Helical Gear Geometry Optimisation (NSGA-II)

Run `helical_gear_nsga2.m`. Writes `pareto_front.png` and prints a feasibility
check for the baseline design before the front is plotted.

## Design variables

| | Symbol | Meaning | Range |
|---|---|---|---|
| 1 | z₁ | pinion tooth count (integer) | 17 – 30 |
| 2 | mₙ | normal module (mm) | 2 – 4 |
| 3 | β | helix angle | 8° – 20° |
| 4 | xₙ₁ | pinion profile shift | 0.2 – 0.5 |
| 5 | xₙ₂ | wheel profile shift | −0.5 – 0.5 |

## Objectives

- `f₁` — tooth-flank contact stress (minimise)
- `f₂` — difference in maximum root bending stress between pinion and wheel (minimise)

## Constraints

Contact fatigue and root bending fatigue (both gears), centre-distance match, total
contact ratio ≥ 2, tip-land thickness (both gears) and minimum profile shift against
undercut (both gears) — ten in total.

Two settings are worth knowing about before you run it:

```matlab
saCoef = 0.25;   % tip-thickness criterion: 0.25·mₙ for through-hardened gears,
                 % 0.4·mₙ for surface-hardened. Changes which designs are feasible.
constr(i,5) = abs(a_actual - a) - 0.5;   % centre-distance tolerance, mm
```

`saCoef` applies to the pinion; the wheel's tip-thickness check is fixed at 0.4·mₙ.

The centre distance is solved through the working pressure angle:

```
inv(αt′) = inv(αt) + 2(xₙ₁ + xₙ₂)·tan(αₙ) / (z₁ + z₂)
a′       = mₙ(z₁ + z₂) / (2cos β) · cos(αt) / cos(αt′)
```

`solve_inv` inverts the involute function by Newton iteration. Using αt in place of
αt′ here removes the profile-shift terms from the equation, which makes the
constraint impossible to satisfy for a shifted pair and pushes the search to pile up
against the helix-angle bound.

## Output

The feasible non-dominated set (at most the population size, 200) plus three labelled
representatives: minimum bending-stress
gap (A), balanced trade-off (B) and minimum contact stress (C).
