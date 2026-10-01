function p = build_vehicle_params()
% BUILD_VEHICLE_PARAMS  All vehicle and simulation parameters in one place.
%
% VEHICLE ASSUMPTIONS:
%   [V1] Mass:   100 kg provisional estimate; replace after mass budget.
%   [V2] S_ref:  0.8 m^2 (user-confirmed FOSTRAD opt.SREF).
%   [V3] Non-axisymmetric waverider: length 2 m, span 0.4 m.
%   [V4] Launch: Mach 8 at h=0 m.
%   [V5] AoA:    +/-10 deg available in the FOSTRAD table.
%   [V6] Limits: h<=30 km, |n|<=20 g, Mach>=3 at actual impact.
%
% OPTIMIZER ARCHITECTURE: one GA optimizes AoA at fixed time knots.
%
% KEY TUNING: ga_* weights rank infeasible trajectories; feasible candidates
% are ranked by impact range. See README for search budget and grid growth.

%% ---- Physical constants ----
p.g       = 9.80665;
p.gam_air = 1.4;
p.earth_radius_m = 6371000;        % spherical surface radius [m]
p.mu_earth = p.g*p.earth_radius_m^2; % matches p.g at h=0 [m^3/s^2]

%% ---- Vehicle ----
p.m     = 100.0;  % provisional estimate, not derived from the STL
p.S_ref = 0.8;    % coefficient reference area, NOT frontal section area
p.L_ref = 2.0;    % provisional moment length; verify FOSTRAD opt.LREF
p.vehicle_length = 2.0;
p.max_span = 0.400000006;
p.max_height_at_widest_section = 0.049780488;
p.area_at_widest_span_station = 0.00873500611; % at x=0.5 m

%% ---- Launch ----
p.h0       = 0.0;
p.M_launch = 8.0;

%% ---- Constraints ----
p.h_max        = 30000;
p.n_max        = 20.0;
p.M_min_impact = 3.0;
[~,~,~,a_sl,~] = atmosphere_1976(0);
p.V_min_impact = p.M_min_impact * a_sl;

%% ---- AoA limits ----
p.alpha_max_deg = 10.0;  % FOSTRAD sweep covers only -10 to +10 deg
p.alpha_max_rad = deg2rad(p.alpha_max_deg);

%% ---- Aero table bounds ----
p.M_min_table = 2.0;
p.M_max_table = 9.0;

%% ---- Time grid (25 AoA knots; first five span 0--15 s) ----
p.n_time_knots = 180;         % change this to tune control resolution
p.n_early_knots = 30;         % includes knots at t=0 and t=15 s
p.early_window_s = 60;
p.reference_t_max = 3600;   % reference integration safety limit [s]
p.time_horizon_override_s = []; % optional manual grid end [s]
p.time_horizon_factor = 1.15;   % grid extends past reference impact
p.seed_results_file = fullfile('results','results_ang40.mat');
p.wp_times = [];             % assigned in setup_waypoints
p.N_wp = p.n_time_knots;
p.N_free_wp = p.N_wp;
p.ga_seed_perturb_deg = 1.2;
p.impact_mach_soft_max = Inf;      % optional preference; see objective
p.ga_w_terminal_mach = 0.2;

p.wp_lb     = -p.alpha_max_deg * ones(1, p.N_wp);
p.wp_ub     =  p.alpha_max_deg * ones(1, p.N_wp);

%% ---- ODE time limit ----
p.T_max = 1800;
p.abort_altitude_m = p.h_max + 10000; % no need to integrate grossly infeasible ascent

%% ---- Phase 1: Exploratory GA ----
% 25 AoA variables. Runtime scales with ga_pop_size * ga_max_gen.
p.ga_pop_size        = 60;
p.ga_max_gen         = 120;
p.ga_stall_gen       = 20;    % best-fitness stall window [generations]
p.ga_stall_tol       = 200;   % [m] minimum improvement to not trigger stall
p.ga_tournament_size = 4;     % low = more exploration
p.ga_crossover_frac  = 0.80;
p.ga_elite_count     = 3;     % small = don't lock in early
p.ga_use_parallel    = true;  % independent GA fitness calls on CPU workers
p.ga_num_workers     = 2;     % conservative laptop default; no GPU needed

p.ga_w_alt = 6.0;
p.ga_w_g   = 6.0;
p.ga_w_vel = 2.5;

%% ---- ODE accuracy ----
p.ode_reltol          = 1e-6;
p.ode_abstol          = 1e-7;
p.ode_max_step        = 2.0;   % faster GA; final result replayed at <=0.5 s

end
