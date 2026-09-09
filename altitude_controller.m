function [alpha_cmd_deg, n_cmd] = altitude_controller(t, V, gamma, x, h, M, q, p)
% ALTITUDE_CONTROLLER  PN altitude-tracking controller (distance-based waypoints)
%
% Determines the commanded AoA to steer toward the current altitude waypoint
% using Proportional Navigation. Waypoints switch on DOWNRANGE DISTANCE, not
% time — which is simpler and more physically meaningful than time-based
% switching: the controller always knows how far ahead the next target is.
%
% Switching rule:
%   wp_idx increments when x (downrange) crosses p.wp_ranges(wp_idx).
%
% Virtual target:
%   The PN target is the actual waypoint position (x_wp, h_target), not a
%   time-projected estimate. dx = x_wp - x is always >= 0 because the switch
%   fires as soon as x >= x_wp. This removes the T_lookahead parameter and
%   simplifies the LOS geometry.
%
% Out-of-reach waypoints:
%   When h_target exceeds the vehicle's current kinematic reach, CL_cmd
%   exceeds the table maximum, and alpha saturates at alpha_max. This is
%   correct behaviour (maximum effort), and the optimizer will learn to avoid
%   placing unreachable waypoints because they waste energy. Saturation
%   during the ascent phase is expected and not a concern.
%
% Inputs:
%   t     - current time [s]  (not used for switching; kept for signature compat)
%   V     - airspeed [m/s]
%   gamma - flight path angle [rad]
%   x     - downrange distance [m]
%   h     - altitude [m]
%   M     - Mach number [-]
%   q     - dynamic pressure [Pa]
%   p     - parameter struct (must include p.wp_ranges, p.wp_alts)
%
% Outputs:
%   alpha_cmd_deg - commanded AoA [deg]
%   n_cmd         - commanded load factor [g] after clamping

%% Step 1 — distance-based waypoint switching
wp_idx   = min(sum(x >= p.wp_ranges) + 1, p.N_wp);
h_target = p.wp_alts(wp_idx);

%% Step 2 — target position (actual waypoint, no time projection)
x_target = p.wp_ranges(wp_idx);

%% Step 3 — LOS geometry
dx = max(x_target - x, p.R_min_threshold);  % guard against dx=0 at waypoint
dh = h_target - h;
R  = sqrt(dx^2 + dh^2);

%% Step 4 — PN normal acceleration command
if R < p.R_min_threshold || V < 1
    a_n_cmd = 0;
else
    lambda = atan2(dh, dx);            % LOS angle [rad]
    V_c    = V * cos(lambda - gamma);  % closing speed [m/s]

    % LOS rate (stationary target, corrected sign):
    %   d(lambda)/dt = V * sin(lambda - gamma) / R
    %
    % Sign: if gamma > lambda (heading above target), sin(lambda-gamma) < 0
    %   → lambda_dot < 0, V_c > 0, a_n = N*V_c*lambda_dot < 0 (pitch down) ✓
    lambda_dot = V * sin(lambda - gamma) / R;
    a_n_cmd    = p.N_nav * V_c * lambda_dot;
end

%% Step 5 — convert to CL, enforce g constraint
%
% Perpendicular force balance (generalised for arbitrary gamma):
%   L = m * (a_n_cmd + g*cos(gamma))
%   CL = L / (q * S_ref)
%
% Note: at high-q low-altitude conditions (sea-level Mach 8), the 15g
% clamp will activate at very small alpha (~1-2 deg). At low-q high-altitude
% conditions (30 km, Mach 5), the alpha limit (15 deg) activates first and
% the resulting n may be well below 15g. The two limits do NOT always
% coincide — which one is active depends entirely on dynamic pressure.
q_safe = max(q, 1.0);
CL_cmd = p.m * (a_n_cmd + p.g * cos(gamma)) / (q_safe * p.S_ref);

n_cmd  = CL_cmd * q_safe * p.S_ref / (p.m * p.g);
n_cmd  = max(-p.n_max, min(p.n_max, n_cmd));
CL_cmd = n_cmd * p.m * p.g / (q_safe * p.S_ref);

%% Step 6 — invert CL → AoA
M_clamped = max(p.M_min_table, min(p.M_max_table, M));

CL_curve = interp2(p.aero.alpha_vec, p.aero.Mach_vec, p.aero.CL_table, ...
                   p.aero.alpha_vec, ...
                   M_clamped * ones(size(p.aero.alpha_vec)), ...
                   'spline');

CL_cmd_c      = max(CL_curve(1), min(CL_curve(end), CL_cmd));
alpha_cmd_deg = interp1(CL_curve, p.aero.alpha_vec, CL_cmd_c, 'linear');
alpha_cmd_deg = max(-p.alpha_max_deg, min(p.alpha_max_deg, alpha_cmd_deg));

if ~isfinite(alpha_cmd_deg), alpha_cmd_deg = 0; end
end
