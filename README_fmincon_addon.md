# Maximum-L/D cruise trajectory search

Place `main_trajectory_fmincon.m`, `ld_equilibrium_altitude.m`, and
`ld_equilibrium_seeds.m` beside the existing simulation files and aerodynamic
CSV. In MATLAB, run:

```matlab
clear main_trajectory_fmincon ld_equilibrium_altitude ld_equilibrium_seeds
which main_trajectory_fmincon -all
main_trajectory_fmincon
```

This version has one SQP search. The previous high-cruise, low-cruise,
zero-lift-hold, and post-peak-reclimb branches are removed. It uses the same
718 AoA knots, 1800 s control horizon, 2500 s ODE limit, and 42 smooth
correction anchors plus launch angle. Change these in `User settings` at the
top of `main_trajectory_fmincon.m`.

The CSV gives the best-L/D AoA (about -0.5 degrees). The script calculates
instantaneous equilibrium altitude by solving

    (1/2) rho(h) V^2 S_ref CL(M, alpha_LD)
      = m [g(h) - V^2/(R_earth + h)]

at a specified airspeed (or Mach). This is the level-flight lift condition
at gamma=0 for the spherical-Earth model. It is not a permanent altitude:
as drag slows an unpowered vehicle, the balance altitude falls. The terminal
prints equilibrium altitudes from Mach 8 through Mach 3.

The new seeds use guided AoA only for an initial pitch turn, then hold the
best-L/D AoA through cruise and transition toward zero-lift AoA at a range
of late descent times. Turn and switch times include near-ceiling naive
trajectories and turns that bring the vehicle close to the calculated
balance altitude at cruise entry. Every complete schedule must pass the
actual open-loop simulator and the finer ODE replay. The script prints the
selected equilibrium-entry seed and a near-ceiling naive comparison when
one is feasible. It does not impose a smoothness, cruise-altitude, or
post-peak reclimb constraint.

SQP maximizes impact range subject to the existing 30 km altitude,
20 g load, AoA table, aerodynamic Mach coverage, and Mach 3 impact-speed
conditions. During the selected seed's long max-L/D cruise, correction
anchors can change its AoA by at most `CRUISE_AOA_WINDOW_DEG = 0.5`.
Outside cruise, `MAX_CORRECTION_DEG = 5` applies. Increase the cruise
window if you want the optimizer to depart further from the L/D optimum.
The launch angle is free within its configured bounds. A failed SQP exit
can still retain a better feasible trial; the final result is replayed at
0.5 s before saving.

Results are saved as `results/results_fmincon_start40_h1800_ldcruise.mat`,
`results/trajectory_fmincon_start40_h1800_ldcruise.csv`, and
`results/ld_cruise_comparison.png`. The `.mat` includes the optimized
trajectory, selected seed, near-ceiling comparison, all generated seed
schedules and their labels, and optimizer metadata. Earlier branch results
remain in `results` but are not warm starts or final candidates.

The max-L/D AoA does not guarantee a perfectly level path or remove all
oscillation. Examine the plotted cruise history: if speed changes enough,
its equilibrium altitude moves. The point-mass model also assumes AoA can
be trimmed and achieved; pitching-moment/control authority is not enforced.
MATLAB was unavailable in this workspace, so run the script on your machine
and check the reported feasibility and finer replay status.
