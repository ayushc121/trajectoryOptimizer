function dy = trajectory_eom(t, y, p)
% TRAJECTORY_EOM  Planar point mass over a non-rotating spherical Earth.
%
% State vector:  y = [V; gamma; x; h]
%   V     - airspeed [m/s]
%   gamma - flight path angle relative to the LOCAL horizontal [rad]
%   x     - surface-arc downrange R_earth*theta [m]
%   h     - radial altitude above the spherical surface [m]
%
% Control:
%   Alpha is piecewise linear in time at p.wp_times, held at the endpoints.
%
% r = R_earth + h, g(r) = mu_earth/r^2:
%   dV/dt     = -D/m - g(r)*sin(gamma)
%   dgamma/dt = (L/m - [g(r) - V^2/r]*cos(gamma))/V
%   dx/dt     = R_earth*V*cos(gamma)/r
%   dh/dt     = V*sin(gamma)
% The V*cos(gamma)/r term is the rotation rate of the local horizontal.

%% Extract state
V     = y(1);
gamma = y(2);
x     = y(3);
h     = max(y(4), 0);      % floor at zero to keep atmosphere valid
r     = p.earth_radius_m + h;
g_r   = p.mu_earth / r^2;

%% Atmosphere
[rho, ~, ~, a, ~] = atmosphere_1976(h);
V = max(V, 1.0);           % guard against near-zero speed at end of flight
M_actual = V / max(a, 1.0);
M = max(p.M_min_table, min(p.M_max_table, M_actual));
q = 0.5 * rho * V^2;

%% Alpha schedule — piecewise linear in elapsed time
alpha_deg = alpha_command(t,M_actual,p);

%% Aerodynamic forces
[CL, CD] = aero_lookup(M, alpha_deg, p.aero);

% CL and CD must describe the same aerodynamic operating point. A structural
% load limit cannot change lift by itself: exceeding p.n_max makes the
% commanded trajectory infeasible (checked in simulate_trajectory).
% This point-mass model assumes the commanded AoA can be achieved and trimmed;
% actual control authority/trim must be established from aerodynamic data.

L = CL * q * p.S_ref;
D = CD * q * p.S_ref;

%% Equations of motion
dV     = -D/p.m - g_r * sin(gamma);
dgamma = (L/p.m - (g_r - V^2/r)*cos(gamma)) / V;
dx     = p.earth_radius_m * V * cos(gamma) / r;
dh     = V * sin(gamma);

dy = [dV; dgamma; dx; dh];
end
