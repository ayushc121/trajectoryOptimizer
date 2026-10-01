# Time-controlled FOSTRAD trajectory optimizer

Keep every `.m` file and `owen_tricon_v1_data.csv` together. Replace the
older unsuffixed versions in your MATLAB working folder. Run:

```matlab
which trajectory_eom -all
which alpha_command -all
run_trajectory_checks
benchmark_ode_step
main_trajectory_optimizer
```

`main_trajectory_optimizer` rebuilds the aero table from the CSV each run.
It expects six headerless columns `[AoA_deg, CL, CD, LD, CMy, Qmax_Wm2]`.
`S_ref=0.8 m^2` matches the reported FOSTRAD coefficient reference area;
mass is still a provisional 100 kg. The data have no Mach axis, so CL/CD/CMy
are repeated over Mach 2–9 as an assumption. The 2D spherical Earth model
uses surface-arc downrange and radial altitude. CMy, Qmax, trim, control
rates, wind, heating, and Earth rotation do not enter its equations.

## Fixed time controls

The GA now runs **once** and optimizes AoA values only. The control is linear
between fixed times and holds the endpoint AoA beyond the last knot. There
are 25 AoA knots by default: five evenly spaced from 0 to 15 s (including
both endpoints), then 20 evenly spaced from 15 s to the horizon (excluding
the duplicate 15 s). Every candidate uses the same times. Consequently a
longer trajectory retains the same early control density, and there is no
remapping or restart between GA runs.

Change these values in `build_vehicle_params.m`:

```matlab
p.n_time_knots = 25;       % total number of AoA commands
p.n_early_knots = 5;       % commands from 0 through early_window_s
p.early_window_s = 15;    % seconds
p.time_horizon_override_s = []; % set a positive time if the probe is unsuitable
```

The setup finds the peak positive interpolated CL/CD in the supplied sweep,
flies that AoA until the **actual Mach** reaches 3, then commands interpolated
zero-lift AoA. If that reference impacts within aerodynamic and atmosphere
coverage, the last control knot is placed at 115% of its flight time. This
is a **timing reference, not an upper bound** on optimized flight time or a
feasible trajectory. If it stops without a valid impact, setup prints a
warning and uses the larger of its elapsed time and the existing 1200 s
`T_max`; inspect the printed horizon or set `time_horizon_override_s` before
running a costly GA. `T_max` is kept beyond the final knot; trajectories may
last past the knot with its final AoA held constant.

The initial population includes zero-lift, maximum-L/D cruise, cruise to
zero-lift transitions, and nine loft/altitude-hold/late-descent profiles
based on local gravity, density, and interpolated CL/CD. Their provisional
feedback law estimates an AoA history, which is sampled at the fixed knots;
all seeds are then simulated with the actual open-loop model. If
`p.seed_results_file` exists and matches reference area and launch angle,
the previous saved AoA history is also projected onto the new time grid and
re-evaluated. It is never assumed feasible merely because it was saved.
Perturbations of these histories fill the remaining initial population.
The script prints how many seeds are feasible and the best seed's range.

A feasible trajectory is ranked by impact range, with the optional terminal
Mach preference off (`p.impact_mach_soft_max = Inf`). Infeasible trajectories
are ranked by altitude, normal load, and missing/slow impact violations.
The hard limits remain 30 km, 20 g, and at least Mach 3 at actual impact.
The GA cannot certify a global optimum, and an instantaneous AoA command
may not be trimmable by a real vehicle.

`p.ga_pop_size=60` and `p.ga_max_gen=120` allow up to 7200 GA candidate
evaluations instead of three such runs. Two CPU workers are enabled by
default; a GPU is not required. Optimization uses a 2 s ODE maximum step,
then the selected schedule is replayed with identical settings and at 0.5 s
for a convergence check. `benchmark_ode_step` now benchmarks time-based
schedules and reports why an incomplete flight stopped. Only a finished
impact validates a claimed range.

The AoA plot now has elapsed time on its horizontal axis. The altitude plot
marks the control knots reached before impact. The final summary prints the
25 time knots and their AoA values. Results are saved after the single GA
completes; files from old range-knot runs can be used as optional seeds,
but their reported range should not be compared directly if their physics
or parameters differ.
