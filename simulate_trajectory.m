function traj = simulate_trajectory(wp_free_alphas, p, detailed)
% SIMULATE_TRAJECTORY  Integrate the trajectory for a given AoA waypoint set
%
% Alpha is a piecewise-linear schedule at fixed p.wp_times. In the optional
% reference-only mode, alpha_command switches at actual Mach 3 instead.
%
% Inputs:
%   wp_free_alphas - [1 x N_wp] AoA values at time knots [deg].
%                   All waypoints are free (no fixed terminal value).
%   p              - parameter struct from build_vehicle_params + setup_waypoints.
%                   p.gamma0 MUST be set by the caller.
%   detailed       - (optional, default true) if false, omits the stored
%                   alpha/M/CMy histories. Loads are still checked.
%
% Outputs:
%   traj - struct:
%     .t          [Nt x 1]  time [s]
%     .V          [Nt x 1]  airspeed [m/s]
%     .gamma      [Nt x 1]  flight path angle [rad]
%     .x          [Nt x 1]  spherical surface-arc downrange [m]
%     .h          [Nt x 1]  radial altitude above spherical Earth [m]
%     .central_angle_rad [Nt x 1]  x / Earth radius [rad]
%     .M          [Nt x 1]  actual Mach (0 in fast mode)
%     .n          [Nt x 1]  load factor [g]
%     .alpha      [Nt x 1]  AoA [deg] (0 in fast mode)
%     .CMy        [Nt x 1]  source pitching-moment coefficient (NaN in fast
%                            mode; diagnostic only, no trim constraint)
%     .x_final    scalar    final downrange [m] (impact only if .landed)
%     .h_max      scalar    max altitude [m]
%     .n_max      scalar    max |load factor| [g]
%     .V_final    scalar    impact speed [m/s]
%     .M_final    scalar    impact Mach [-]
%     .t_flight   scalar    total flight time [s]
%     .landed     logical   integration ended at the ground event
%     .aero_valid logical   trajectory remained inside the aero-table grid
%     .feasible   logical   landed, valid aero, and constraints satisfied
%     .c_viol     [3 x 1]   constraint violations (c <= 0 satisfied):
%                             c(1) = h_max - p.h_max [m]
%                             c(2) = n_max - p.n_max [g]
%                             c(3) = V_min - V_final [m/s]
%     .wp_alphas  [1 x N_wp] AoA schedule at p.wp_times [deg]
%     .gamma0     scalar    launch angle [rad]

if nargin < 3, detailed = true; end
if ~isfield(p,'earth_radius_m') || ~isfield(p,'mu_earth') || ...
        ~isfinite(p.earth_radius_m) || p.earth_radius_m <= 0 || ...
        ~isfinite(p.mu_earth) || p.mu_earth <= 0
    error('simulate_trajectory:earthModel', ...
        'A positive spherical Earth radius and gravitational parameter are required.');
end

%% ---- Store alpha schedule in p for trajectory_eom ----
p.wp_alphas = wp_free_alphas(:)';   % [1 x N_wp]
if numel(p.wp_alphas) ~= numel(p.wp_times) || ...
        numel(p.wp_times) < 2 || p.wp_times(1) ~= 0 || ...
        any(~isfinite(p.wp_times)) || any(diff(p.wp_times) <= 0)
    error('simulate_trajectory:waypoints', ...
        'AoA values and strictly increasing time knots starting at zero must match.');
end
if ~isfield(p.aero, 'S_ref') || ~isfinite(p.aero.S_ref) || ...
        abs(p.aero.S_ref - p.S_ref) > 1e-8 * p.S_ref
    error('simulate_trajectory:referenceArea', ...
        'Aero reference area must be supplied and match p.S_ref.');
end
if numel(p.aero.Mach_vec) < 2 || numel(p.aero.alpha_vec) < 2 || ...
        any(diff(p.aero.Mach_vec) <= 0) || any(diff(p.aero.alpha_vec) <= 0) || ...
        ~isequal(size(p.aero.CL_table), [numel(p.aero.Mach_vec), numel(p.aero.alpha_vec)]) || ...
        ~isequal(size(p.aero.CD_table), size(p.aero.CL_table)) || ...
        any(~isfinite(p.aero.CL_table(:))) || any(~isfinite(p.aero.CD_table(:))) || ...
        any(p.aero.CD_table(:) <= 0)
    error('simulate_trajectory:aeroTable', 'Aero grid or coefficients are invalid.');
end

%% ---- Initial conditions ----
[~,~,~,a0,~] = atmosphere_1976(p.h0);
V0 = p.M_launch * a0;
y0 = [V0; p.gamma0; 0.0; p.h0];    % [V; gamma; x; h]

%% ---- ODE options ----
% Both modes MUST integrate identical equations with identical tolerances.
% The detailed flag only requests extra output reconstruction. Previously the
% GA's 2 s maximum step and the report's 0.5 s step selected different paths.
ode_opts = odeset('RelTol', p.ode_reltol, 'AbsTol', p.ode_abstol, ...
                  'MaxStep', p.ode_max_step, ...
                  'Events', @(t,y) ground_impact_event(t,y,p));

%% ---- Integrate ----
eom_handle = @(t, y) trajectory_eom(t, y, p);
[t_out, y_out, t_event, ~, i_event] = ode45(eom_handle, [0, p.T_max], y0, ode_opts);
landed = any(i_event == 1 & t_event > 0);  % event 2 is insufficient airspeed

%% ---- Extract states ----
V_h   = y_out(:,1);
gam_h = y_out(:,2);
x_h   = y_out(:,3);
h_h   = max(y_out(:,4), 0);
nT    = length(t_out);

%% ---- Constraint quantities ----
h_max_traj = max(h_h);

% Reconstruct the same (uncapped) lift used by trajectory_eom at every node.
% Finite differencing gamma on an adaptive mesh can miss a load peak.
[rho_hist,~,~,a_hist,~] = atmosphere_1976(h_h);
M_actual = V_h ./ a_hist;
alpha_used = alpha_command(t_out,M_actual,p);
CMy_hist = nan(nT,1);
if isfield(p.aero, 'mach_invariant_assumed') && p.aero.mach_invariant_assumed
    % This CSV repeats the AoA sweep at every Mach. Evaluate the load
    % history in one vector operation rather than thousands of scalar
    % atmosphere and 2D table calls for every GA candidate.
    alpha_lookup = max(p.aero.alpha_vec(1), ...
                       min(p.aero.alpha_vec(end),alpha_used));
    CL_hist = interp1(p.aero.alpha_vec, p.aero.CL_table(1,:), ...
                      alpha_lookup, 'linear');
    if detailed && isfield(p.aero, 'CMy_table')
        CMy_hist = interp1(p.aero.alpha_vec, p.aero.CMy_table(1,:), ...
                            alpha_lookup, 'linear');
    end
else
    CL_hist = zeros(nT,1);
    for i = 1:nT
        if detailed && isfield(p.aero, 'CMy_table')
            [CL_hist(i),~,CMy_hist(i)] = aero_lookup(M_actual(i),alpha_used(i),p.aero);
        else
            [CL_hist(i),~] = aero_lookup(M_actual(i),alpha_used(i),p.aero);
        end
    end
end
q_hist = 0.5 * rho_hist .* max(V_h,1).^2;
n_used = CL_hist .* q_hist * (p.S_ref/(p.m*p.g));
n_max_traj = max(abs(n_used));
aero_valid = all(M_actual >= p.aero.Mach_vec(1) & ...
                 M_actual <= p.aero.Mach_vec(end)) && ...
             all(alpha_used >= p.aero.alpha_vec(1) & ...
                 alpha_used <= p.aero.alpha_vec(end));

% Terminal speed and Mach
[~,~,~,a_fin,~] = atmosphere_1976(h_h(end));
V_fin  = V_h(end);
M_fin  = V_fin / max(a_fin, 1);

%% ---- Detailed per-step reconstruction ----
if detailed
    M_hist     = M_actual;
    alpha_hist = alpha_used;
else
    M_hist     = zeros(nT,1);
    alpha_hist = zeros(nT,1);
end
n_hist = n_used;

%% ---- Constraint vector (c <= 0 satisfied) ----
c_viol    = zeros(3,1);
c_viol(1) = h_max_traj - p.h_max;          % altitude ceiling [m]
c_viol(2) = n_max_traj - p.n_max;          % g-load [g]
c_viol(3) = p.V_min_impact - V_fin;        % terminal speed [m/s]

% Preserve the 3-component interface expected by the existing optimiser.
% Incomplete or out-of-table runs cannot satisfy its terminal-speed test.
if ~landed || ~aero_valid
    c_viol(3) = max(c_viol(3), 1e6);
end

feasible = landed && aero_valid && all(c_viol <= 1e-2);

%% ---- Package output ----
traj.t        = t_out;
traj.V        = V_h;
traj.gamma    = gam_h;
traj.x        = x_h;
traj.central_angle_rad = x_h / p.earth_radius_m;
traj.h        = h_h;
traj.M        = M_hist;
traj.n        = n_hist;
traj.alpha    = alpha_hist;
traj.CMy      = CMy_hist;
traj.x_final  = x_h(end);
traj.h_max    = h_max_traj;
traj.n_max    = n_max_traj;
traj.V_final  = V_fin;
traj.M_final  = M_fin;
traj.t_flight = t_out(end);
traj.landed   = landed;
traj.aero_valid = aero_valid;
traj.feasible = feasible;
traj.c_viol   = c_viol;
traj.wp_alphas = p.wp_alphas;
traj.gamma0   = p.gamma0;
end

%% ======================================================================
function [value, isterminal, direction] = ground_impact_event(~, y, p)
value      = [y(4); y(1) - 1; p.abort_altitude_m - y(4)];
isterminal = [1; 1; 1]; % ground, speed guard, grossly infeasible altitude
direction  = [-1; -1; -1];
end
