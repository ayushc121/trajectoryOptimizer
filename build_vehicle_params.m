function p = build_vehicle_params()
% BUILD_VEHICLE_PARAMS  All vehicle and simulation parameters in one place.
%
% VEHICLE ASSUMPTIONS:
%   [V1] Mass:   100 kg (post-sabot).
%   [V2] S_ref:  pi*(0.075)^2 m^2 (D=0.15 m).
%   [V3] Launch: Mach 8 at h=0 m.
%   [V4] AoA:    +/-15 deg authority.
%   [V5] Limits: h<=30 km, n<=20 g (CL-clamped in EOM), Mach>=3 at impact.
%
% OPTIMIZER ARCHITECTURE:
%   Phase 1 — Exploratory GA (alpha + waypoint fracs, 22 vars).
%   Phase 2 — Polish GA (alpha only, 11 vars, higher penalty pressure,
%             seeded from Phase 1 final population).
%   Phase 3 — fmincon SQP (optional, skip with N_random_starts=0).
%
% KEY TUNING:
%   ga_w_alt / ga_w_vel
%     Phase 1 penalty weights. If Phase 2 consistently fails to achieve
%     feasibility, raise these first before touching polish parameters.
%
%   ga_polish_penalty_mult
%     Multiplier applied to ga_w_alt and ga_w_vel for Phase 2.
%     2x (default) means Phase 2 is twice as aggressive on constraints.
%     Raise to 3x-4x if Phase 2 still returns infeasible results.
%
%   ga_polish_alt_margin / ga_polish_vel_margin
%     Safety buffers added to constraint targets in Phase 2.
%     Phase 2 aims for h_max - alt_margin and V_min + vel_margin.
%     Ensures a comfortable margin before Phase 3 sees the solution.
%
%   range_est_factor
%     Less critical since GA optimises waypoint positions, but sets the
%     outer x-search bound. Keep near 1.0-1.2.

%% ---- Physical constants ----
p.g       = 9.80665;
p.gam_air = 1.4;

%% ---- Vehicle ----
p.m     = 100.0;
p.S_ref = pi*(0.075)^2;
p.L_ref = 0.15;

%% ---- Launch ----
p.h0       = 0.0;
p.M_launch = 8.0;

%% ---- Constraints ----
p.h_max        = 30000;
p.n_max        = 20.0;
p.M_min_impact = 3.0;
[~,~,~,a_sl,~] = atmosphere_1976(0);
p.V_min_impact = p.M_min_impact * a_sl;

% Phase 1 GA soft-constraint targets (margins allow Phase 2 to tighten further)
p.h_max_ga = p.h_max - 500;
p.V_min_ga = p.V_min_impact + 60;

%% ---- AoA limits ----
p.alpha_max_deg = 15.0;
p.alpha_max_rad = deg2rad(15.0);

%% ---- Aero table bounds ----
p.M_min_table = 2.0;
p.M_max_table = 9.0;

%% ---- Waypoint grid ----
p.range_est_factor = 1.1;

p.wp_fracs = [0.02, 0.05, 0.10, 0.20, 0.35, ...
              0.50, 0.65, 0.75, 0.85, 0.92, 0.97];
p.N_wp      = numel(p.wp_fracs);
p.N_free_wp = p.N_wp;

p.ga_frac_delta = 0.12;

% Placeholders — overwritten by setup_waypoints
p.frac_lb   = max(0.005, p.wp_fracs - p.ga_frac_delta);
p.frac_ub   = min(0.995, p.wp_fracs + p.ga_frac_delta);
p.wp_ranges = p.wp_fracs * 400e3;
p.wp_lb     = -p.alpha_max_deg * ones(1, p.N_wp);
p.wp_ub     =  p.alpha_max_deg * ones(1, p.N_wp);

%% ---- ODE time limit ----
p.T_max = 600;

%% ---- Phase 1: Exploratory GA ----
% 22 vars (alpha + fracs). Best-fitness stall via custom OutputFcn.
% pop=60, gen=120 -> ~7200 ODE evals -> 20-40 min.
% pop=40, gen=80  -> ~3200 ODE evals ->  8-15 min.
p.ga_pop_size        = 60;
p.ga_max_gen         = 120;
p.ga_stall_gen       = 20;    % best-fitness stall window [generations]
p.ga_stall_tol       = 200;   % [m] minimum improvement to not trigger stall
p.ga_tournament_size = 4;     % low = more exploration
p.ga_crossover_frac  = 0.80;
p.ga_elite_count     = 3;     % small = don't lock in early

p.ga_w_alt = 6.0;
p.ga_w_g   = 0.0;   % g-load: CL-clamped in EOM
p.ga_w_vel = 2.5;

%% ---- Phase 2: Polish GA ----
% 11 vars (alpha only). Seeded from Phase 1 final population.
% Higher penalty weights and tighter constraint margins drive feasibility
% while preserving range (range term still in objective — no degenerate
% "crash early" solution is possible).
p.ga_polish_penalty_mult    = 2.5;   % multiplier on Phase 1 penalty weights
p.ga_polish_alt_margin      = 300;   % [m]   h_max_ga_p2 = h_max - this value
p.ga_polish_vel_margin      = 30;    % [m/s] V_min_ga_p2 = V_min + this value
p.ga_polish_max_gen         = 60;    % generation budget for Polish GA
p.ga_polish_stall_gen       = 15;    % best-fitness stall window
p.ga_polish_elite_count     = 8;     % high — preserve good Phase 1 solutions
p.ga_polish_tournament_size = 6;     % more selective than Phase 1

%% ---- Phase 3: fmincon SQP (optional) ----
% Skip by setting N_random_starts = 0 in main.
p.fmin_algorithm = 'sqp';
p.fmin_max_iter  = 600;
p.fmin_tol_fun   = 1.0;
p.fmin_tol_con   = 0.005;
p.fmin_fd_type   = 'central';
p.fmin_fd_step   = 0.5;
p.phase3_restart_delta = 2.0;

%% ---- ODE accuracy ----
p.ode_reltol          = 1e-4;
p.ode_abstol          = 1e-5;
p.ode_max_step        = 2.0;

p.ode_reltol_detail   = 1e-6;
p.ode_abstol_detail   = 1e-7;
p.ode_max_step_detail = 0.5;

end
