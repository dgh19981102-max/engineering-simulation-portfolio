# Pickup-and-Delivery Routing with Pairing and Type Constraints

Run `pickup_delivery_routing.m`. Draws the best route found and a convergence
comparison across four algorithms.

## Problem

37 nodes: one depot, 18 pickup points and 18 matching delivery points, across 4 load
types.

- **One part at a time** — every pickup is followed directly by a delivery point
- **Type matching** — that delivery point must be of the same type as the part picked up

Node coordinates are generated procedurally by `generate_array_coordinates`
(grid layout with per-type offsets), so the instance is reproducible and easy to
replace with your own data.

## Algorithms

| Algorithm | Final tour length | vs. greedy |
|---|---|---|
| Hybrid metaheuristic | 5367.7 | −1.66% |
| Tabu Search | 5368.3 | −1.65% |
| Ant Colony Optimisation | 5396.4 | −1.14% |
| Greedy (baseline) | 5458.4 | — |

All four run on the same instance and the same iteration budget (100), so the numbers
are directly comparable. Values are from MATLAB R2023a with the script's `rng(6)`.
The hybrid and tabu results differ by 0.01%, which is well inside run-to-run noise — on this instance they are tied, and tabu search gets there
with a simpler implementation.

## Swapping in your own data

Replace `generate_array_coordinates` with your coordinate matrix (n×2) and adjust
`n_total`, `n_pairs` and `n_types` at the top of the script (`type_counts` is drawn
from them). The constraint machinery is
independent of how the coordinates are produced.

## Note

`pdist2` comes from the Statistics and Machine Learning Toolbox. If you do not have
it, replace the call with this base-MATLAB line (R2016b+):

```matlab
dist_matrix = sqrt(sum((permute(coordinates, [1 3 2]) - ...
                        permute(coordinates, [3 1 2])).^2, 3));
```

It matches `pdist2` exactly, with exact zeros on the diagonal. The faster
Gram-matrix form `sqrt(max(0, s + s' - 2*G))` does not guarantee exact diagonal
zeros, and a tiny non-zero value there stops `eta(isinf(eta)) = 0` from clearing
the self-loops — set the diagonal to 0 if you use it.
