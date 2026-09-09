function traj = simulate_trajectory(wp_free_alphas, p, detailed)
% SIMULATE_TRAJECTORY  Integrate the trajectory for a given AoA waypoint set
%
% Alpha is commanded as a piecewise-linear schedule in downrange x, defined
% by the waypoint positions p.wp_ranges and the decision variables
% wp_free_alphas. There is no PN controller — alpha is the direct control.
%
% Inputs:
%   wp_free_alphas - [1 x N_wp] AoA values at each waypoint position [deg].
%                   All waypoints are free (no fixed terminal value).
%   p              - parameter struct from build_vehicle_params + setup_waypoints.
%                   p.gamma0 MUST be set by the caller.
%   detailed       - (optional, default true) if false, skips per-step
%                   alpha/n/M reconstruction to speed up optimisation.
%
% Outputs:
%   traj - struct:
%     .t          [Nt x 1]  time [s]
%     .V          [Nt x 1]  airspeed [m/s]
%     .gamma      [Nt x 1]  flight path angle [rad]
%     .x          [Nt x 1]  downrange [m]
%     .h          [Nt x 1]  altitude [m]
%     .M          [Nt x 1]  Mach (0 in fast mode)
%     .n          [Nt x 1]  load factor [g]
%     .alpha      [Nt x 1]  AoA [deg] (0 in fast mode)
%     .x_final    scalar    impact downrange [m]
%     .h_max      scalar    max altitude [m]
%     .n_max      scalar    max |load factor| [g]
%     .V_final    scalar    impact speed [m/s]
%     .M_final    scalar    impact Mach [-]
%     .t_flight   scalar    total flight time [s]
%     .feasible   logical   all three constraints satisfied
%     .c_viol     [3 x 1]   constraint violations (c <= 0 satisfied):
%                             c(1) = h_max - 30000  [m]
%                             c(2) = n_max - 15     [g]
%                             c(3) = V_min - V_final [m/s]
%     .wp_alphas  [1 x N_wp] AoA schedule used [deg]
%     .gamma0     scalar    launch angle [rad]

if nargin < 3, detailed = true; end

%% ---- Store alpha schedule in p for trajectory_eom ----
p.wp_alphas = wp_free_alphas(:)';   % [1 x N_wp]

%% ---- Initial conditions ----
[~,~,~,a0,~] = atmosphere_1976(p.h0);
V0 = p.M_launch * a0;
y0 = [V0; p.gamma0; 0.0; p.h0];    % [V; gamma; x; h]

%% ---- ODE options ----
if detailed
    ode_opts = odeset( ...
        'RelTol',  p.ode_reltol_detail, ...
        'AbsTol',  p.ode_abstol_detail, ...
        'MaxStep', p.ode_max_step_detail, ...
        'Events',  @ground_impact_event);
else
    ode_opts = odeset( ...
        'RelTol',  p.ode_reltol, ...
        'AbsTol',  p.ode_abstol, ...
        'MaxStep', p.ode_max_step, ...
        'Events',  @ground_impact_event);
end

%% ---- Integrate ----
eom_handle = @(t, y) trajectory_eom(t, y, p);
[t_out, y_out] = ode45(eom_handle, [0, p.T_max], y0, ode_opts);

%% ---- Extract states ----
V_h   = y_out(:,1);
gam_h = y_out(:,2);
x_h   = y_out(:,3);
h_h   = max(y_out(:,4), 0);
nT    = length(t_out);

%% ---- Constraint quantities ----
h_max_traj = max(h_h);

% Load factor: n = (V*dgamma/dt + g*cos(gamma)) / g
% Central-difference estimate from ODE output — valid at all interior points.
dgam_dt    = gradient(gam_h, t_out);
n_approx   = (V_h .* dgam_dt + p.g * cos(gam_h)) / p.g;
n_max_traj = max(abs(n_approx));

% Terminal speed and Mach
[~,~,~,a_fin,~] = atmosphere_1976(h_h(end));
V_fin  = V_h(end);
M_fin  = V_fin / max(a_fin, 1);

%% ---- Detailed per-step reconstruction ----
if detailed
    M_hist     = zeros(nT,1);
    n_hist     = zeros(nT,1);
    alpha_hist = zeros(nT,1);

    for i = 1:nT
        [rho_i,~,~,a_i,~] = atmosphere_1976(h_h(i));
        M_i  = V_h(i) / max(a_i, 1);
        M_i  = max(p.M_min_table, min(p.M_max_table, M_i));
        q_i  = max(0.5 * rho_i * V_h(i)^2, 1e-3);

        % Alpha from schedule (same logic as trajectory_eom)
        x_s       = max(p.wp_ranges(1), min(p.wp_ranges(end), x_h(i)));
        alpha_i   = interp1(p.wp_ranges, p.wp_alphas, x_s, 'linear');
        alpha_i   = max(-p.alpha_max_deg, min(p.alpha_max_deg, alpha_i));

        [CL_i, ~] = aero_lookup(M_i, alpha_i, p.aero);
        n_i       = CL_i * q_i * p.S_ref / (p.m * p.g);

        M_hist(i)     = M_i;
        n_hist(i)     = n_i;
        alpha_hist(i) = alpha_i;
    end
    n_max_traj = max(abs(n_hist));   % more accurate than fd estimate
else
    M_hist     = zeros(nT,1);
    n_hist     = n_approx;
    alpha_hist = zeros(nT,1);
end

%% ---- Constraint vector (c <= 0 satisfied) ----
c_viol    = zeros(3,1);
c_viol(1) = h_max_traj - p.h_max;          % altitude ceiling [m]
c_viol(2) = n_max_traj - p.n_max;          % g-load [g]
c_viol(3) = p.V_min_impact - V_fin;        % terminal speed [m/s]

% If the vehicle never lands (T_max exceeded), force infeasibility
if h_h(end) > 200
    c_viol(3) = c_viol(3) + 1e6;
end

feasible = all(c_viol <= 1e-2);

%% ---- Package output ----
traj.t        = t_out;
traj.V        = V_h;
traj.gamma    = gam_h;
traj.x        = x_h;
traj.h        = h_h;
traj.M        = M_hist;
traj.n        = n_hist;
traj.alpha    = alpha_hist;
traj.x_final  = x_h(end);
traj.h_max    = h_max_traj;
traj.n_max    = n_max_traj;
traj.V_final  = V_fin;
traj.M_final  = M_fin;
traj.t_flight = t_out(end);
traj.feasible = feasible;
traj.c_viol   = c_viol;
traj.wp_alphas = p.wp_alphas;
traj.gamma0   = p.gamma0;
end

%% ======================================================================
function [value, isterminal, direction] = ground_impact_event(~, y)
value      = y(4);
isterminal = 1;
direction  = -1;
end
