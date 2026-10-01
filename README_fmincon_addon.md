# Local SQP optimizer add-on

Add `main_trajectory_fmincon.m` and `fmincon_time_seeds.m` to the folder
containing the current 15-file time-based trajectory package. Existing files
are untouched. MATLAB Optimization Toolbox (`fmincon`) is required.

In MATLAB:

```matlab
run_trajectory_checks
main_trajectory_fmincon
```

At the top of `main_trajectory_fmincon.m`, `CONTROL_HORIZON_S = 1200` places
the last AoA control knot at 1200 s, while `MAX_FLIGHT_TIME_S = 1800`
permits integration to continue until 1800 s. The last AoA is held after
1200 s. Raise both if you need active control beyond 1200 s. The existing
GA script reads `p.T_max` in `build_vehicle_params.m` and then
`setup_waypoints.m` sets it to at least 105% of the knot horizon. To run the
GA for longer, raise `p.T_max` there. The separate `reference_t_max` controls
only the best-L/D timing probe.

The add-on makes the same physics-guided seeds as the current GA and also
projects `p.seed_results_file` (normally `results/results_ang40.mat`) onto
the current time grid. It integrates all seeds under current conditions,
chooses the longest feasible seed that agrees at 2 s and 0.5 s ODE steps,
and invokes constrained SQP `fmincon` from it. If none passes, the script
stops with an explicit message. It never starts SQP from an infeasible
trajectory. A completed feasible GA run can provide a useful saved seed.

The optimization minimizes minus impact range subject to altitude, load,
and terminal speed inequalities. The simulator marks non-impact and
out-of-table flights with a large constraint violation. The AoA bounds are
unchanged. The objective and constraints share a trajectory cache, so SQP
uses serial finite differences; it does not need a GPU or worker pool.
Change `MAX_ITER`, `MAX_EVAL`, and `FD_STEP_DEG` at the top of the new main
file to adjust runtime and finite difference behavior.

Best feasible improvements are checkpointed in
`results/fmincon_checkpoint_ang40.mat`. The final result goes to
`results/results_fmincon_ang40.mat` and
`results/trajectory_fmincon_ang40.csv`, leaving the GA result intact. The
selected schedule is checked at 0.5 s; if its improvement fails that check,
the script retains the verified starting seed.

SQP is a local search. ODE events, maximum altitude/load operations, and a
finite control grid can make gradients noisy or nonsmooth. A feasible start
and a solver success message do not establish a global maximum. MATLAB and
Octave were unavailable when this add-on was prepared, so run the included
checks on your machine and inspect the reported replay status.
