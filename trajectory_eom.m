function dy = trajectory_eom(t, y, p)
% TRAJECTORY_EOM  2D point-mass equations of motion with direct AoA schedule
%
% State vector:  y = [V; gamma; x; h]
%   V     - airspeed [m/s]
%   gamma - flight path angle [rad]
%   x     - downrange distance [m]
%   h     - altitude [m]
%
% Control:
%   Alpha is a piecewise-linear function of downrange x, defined by the
%   waypoint schedule (p.wp_ranges, p.wp_alphas). For x outside the grid,
%   the nearest endpoint value is held constant (zero extrapolation slope).
%   This replaces the PN altitude-tracking controller entirely.
%
% Equations of motion (2D point mass, flight-path-angle formulation):
%   dV/dt     = (-D - m*g*sin(gamma)) / m
%   dgamma/dt = ( L - m*g*cos(gamma)) / (m*V)
%   dx/dt     = V * cos(gamma)
%   dh/dt     = V * sin(gamma)

%% Extract state
V     = y(1);
gamma = y(2);
x     = y(3);
h     = max(y(4), 0);      % floor at zero to keep atmosphere valid

%% Atmosphere
[rho, ~, ~, a, ~] = atmosphere_1976(h);
V = max(V, 1.0);           % guard against near-zero speed at end of flight
M = V / max(a, 1.0);
M = max(p.M_min_table, min(p.M_max_table, M));
q = 0.5 * rho * V^2;

%% Alpha schedule — piecewise linear in downrange x
% Clamp x to the waypoint grid before interpolating so the alpha holds
% at the endpoint value outside the scheduled range (no extrapolation).
x_sched   = max(p.wp_ranges(1), min(p.wp_ranges(end), x));
alpha_deg = interp1(p.wp_ranges, p.wp_alphas, x_sched, 'linear');
alpha_deg = max(-p.alpha_max_deg, min(p.alpha_max_deg, alpha_deg));

%% Aerodynamic forces
[CL, CD] = aero_lookup(M, alpha_deg, p.aero);

% Enforce structural g-limit as a physical saturation.
% At high dynamic pressure (Mach 8 near sea level), even small alpha
% generates far more than 15g. Clamp CL to the value that produces
% exactly n_max — the flight control system would saturate here
% regardless of commanded alpha.
% This also removes g-load from the optimizer's active constraint set,
% since it is now always satisfied by construction.
q_safe  = max(q, 1.0);
CL_lim  = p.n_max * p.m * p.g / (q_safe * p.S_ref);
CL      = max(-CL_lim, min(CL_lim, CL));

L = CL * q * p.S_ref;
D = CD * q * p.S_ref;

%% Equations of motion
dV     = (-D  - p.m * p.g * sin(gamma)) / p.m;
dgamma = ( L  - p.m * p.g * cos(gamma)) / (p.m * V);
dx     = V * cos(gamma);
dh     = V * sin(gamma);

dy = [dV; dgamma; dx; dh];
end