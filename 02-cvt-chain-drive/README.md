# CVT Chain Drive — Multibody vs. Continuum Cross-Validation

Run `run_comparison.m`. Solves both models over one operating condition and draws a
six-panel comparison.

## Files

| File | What it is |
|---|---|
| `MultibodyModel.m` | discrete model — 78 chain plates resolved individually |
| `CMMModel.m` | continuum model — chain treated as a distributed medium |
| `run_comparison.m` | entry point: runs both, plots the comparison |

## Compared quantities

Normal force P, sliding angle ψ and chain tension F, over the full time history.

## Model parameters

Pitch radius, sheave half-angle, pin count, speed and the other operating parameters
are set in the `params` struct at the top of `run_comparison.m`; the friction
coefficient and sheave-deformation amplitude are fixed inside the model classes. Contact forces, friction and centrifugal effects are coupled
in the multibody side; the continuum side derives tension and sliding-angle
distributions analytically.

## Reading the result

Sliding-angle phase and chain-tension trend agree between the two models. Peak normal
force does not: roughly 5.3×10⁴ N from the multibody model against 2.5×10³ N from the
continuum model. The discrete model resolves the impact as an individual plate enters
mesh; the continuum model averages it away. If you need peak contact loads for a
fatigue or wear estimate, the continuum model will understate them by more than an
order of magnitude.
