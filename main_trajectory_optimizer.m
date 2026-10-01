%% MAIN_TRAJECTORY_OPTIMIZER.M
% ==========================================================================
% Waverider — 2D Point-Mass Trajectory Optimiser
% FY26 UCAH Hypersonic Design Challenge
% ==========================================================================
%
% WORKFLOW: one physics-seeded GA on fixed time knots. fmincon is disabled.

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));

%% ---- USER SETTINGS -------------------------------------------------------

LAUNCH_ANGLE_DEG = 40;
N_RANDOM_STARTS  = 0;   % kept for call compatibility; fmincon is disabled
RNG_SEED = 42;          % repeatable exploration; change for independent runs
rng(RNG_SEED, 'twister');

SAVE_RESULTS = true;
RESULTS_DIR  = 'results';

%% ---- INITIALISATION ------------------------------------------------------

if SAVE_RESULTS && ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

p        = build_vehicle_params();
p.gamma0 = deg2rad(LAUNCH_ANGLE_DEG);

AERO_CSV  = fullfile(fileparts(mfilename('fullpath')), 'owen_tricon_v1_data.csv');
AERO_FILE = fullfile(fileparts(mfilename('fullpath')), 'aero_tables.mat');
% Always rebuild from the CSV: an old CBAero placeholder .mat must not be
% silently reused after changing geometry and reference area.
generate_aero_tables(AERO_FILE, AERO_CSV, p.S_ref, p.L_ref);
aero_data     = load(AERO_FILE);
p.aero        = aero_data;
p.M_min_table = min(aero_data.Mach_vec);
p.M_max_table = max(aero_data.Mach_vec);

p = setup_waypoints(p);

% A small process pool uses the laptop CPU. This objective does not use a GPU.
if p.ga_use_parallel
    if isempty(ver('parallel'))
        error('main_trajectory_optimizer:parallelToolbox', ...
            'Install Parallel Computing Toolbox or set p.ga_use_parallel=false.');
    end
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('local', p.ga_num_workers);
    end
    fprintf('GA parallel pool: %d CPU workers.\n', pool.NumWorkers);
end

fprintf('==========================================================================\n');
fprintf('  FY26 UCAH — Waverider Trajectory Optimiser (FOSTRAD AoA sweep)\n');
fprintf('==========================================================================\n');
fprintf('  Launch angle    : %d deg\n',   LAUNCH_ANGLE_DEG);
fprintf('  Earth model     : spherical, non-rotating; R = %.0f km\n', ...
        p.earth_radius_m/1e3);
fprintf('  Vehicle         : %.2f m length, %.3f m max span, %.1f kg assumed mass\n', ...
        p.vehicle_length, p.max_span, p.m);
fprintf('  FOSTRAD Sref    : %.4f m^2 | AoA data: %.1f to %.1f deg\n', ...
        p.S_ref, min(p.aero.alpha_vec), max(p.aero.alpha_vec));
fprintf('  GA runs         : 1  |  seed: %d\n', RNG_SEED);
fprintf('  AoA knots       : %d  |  first %d by %.1f s\n', ...
        p.N_free_wp,p.n_early_knots,p.early_window_s);
fprintf('  Control horizon : %.1f s  |  zero-lift reference: %.0f km\n', ...
        p.wp_times(end), p.ball.x_ball/1e3);
if ~p.reference.impact_valid
    fprintf('  Reference status: no valid impact; horizon is a fallback, not a flight-time bound.\n');
end
fprintf('  V_min impact    : %.1f m/s (Mach %.1f)\n', p.V_min_impact, p.M_min_impact);
fprintf('  Feasibility weights: altitude %.1f  g-load %.1f  speed %.1f\n', ...
        p.ga_w_alt, p.ga_w_g, p.ga_w_vel);
fprintf('==========================================================================\n\n');

%% ---- OPTIMISE ------------------------------------------------------------

[wp_opt, range_opt, info] = hybrid_optimize(p, N_RANDOM_STARTS);

%% ---- DETAILED SIMULATION -------------------------------------------------

traj = simulate_trajectory(wp_opt, p, true);
if abs(range_opt - traj.x_final) > 1 || ...
        abs(info.traj_ga.h_max - traj.h_max) > 1 || ...
        abs(info.traj_ga.n_max - traj.n_max) > 0.02 || ...
        info.traj_ga.feasible ~= traj.feasible
    error('main_trajectory_optimizer:replayMismatch', ...
        'GA and detailed replay disagree; result was not saved.');
end
% Verify convergence once at a finer time step. Preserve the GA trajectory
% as the report, but reject its feasibility claim if the two disagree.
p_check = p;
p_check.ode_max_step = min(0.5, p.ode_max_step/2);
p_check.ode_reltol = p.ode_reltol/2;
check = simulate_trajectory(wp_opt, p_check, false);
traj.numerically_converged = check.feasible == traj.feasible && ...
    check.landed == traj.landed && ...
    abs(check.x_final-traj.x_final) <= max(100,0.001*traj.x_final) && ...
    abs(check.h_max-traj.h_max) <= 50 && ...
    abs(check.n_max-traj.n_max) <= 0.2;
if ~traj.numerically_converged
    warning('main_trajectory_optimizer:convergence', ...
        'Half-step replay differs; refine the ODE/waypoint grid before trusting this result.');
    traj.feasible = false;
end
fprintf('\n  *** Best evaluated range: %.1f km  |  feasible: %d  |  converged: %d ***\n\n', ...
        traj.x_final/1e3, traj.feasible, traj.numerically_converged);

fprintf('==========================================================================\n');
fprintf('  FINAL TRAJECTORY SUMMARY\n');
fprintf('==========================================================================\n');
fprintf('  Surface range  : %.1f km\n',   traj.x_final/1e3);
fprintf('  Central angle  : %.2f deg\n', rad2deg(traj.x_final/p.earth_radius_m));
fprintf('  Max altitude   : %.1f km\n',   traj.h_max/1e3);
fprintf('  Flight time    : %.1f s\n',    traj.t_flight);
fprintf('  Terminal Mach  : %.3f\n',      traj.M_final);
fprintf('  Terminal speed : %.1f m/s  (min %.1f)\n', traj.V_final, p.V_min_impact);
fprintf('  Max g-load     : %.2f g\n',    traj.n_max);
fprintf('  Feasible       : %d\n',        traj.feasible);
fprintf('  Constraints (c<=0 satisfied):\n');
fprintf('    Altitude (m) : %+.1f\n',   traj.c_viol(1));
fprintf('    G-load   (g) : %+.3f\n',   traj.c_viol(2));
fprintf('    Speed  (m/s) : %+.1f\n',  traj.c_viol(3));
fprintf('--------------------------------------------------------------------------\n');
fprintf('  AoA commands [deg]: '); fprintf('%+.1f  ', wp_opt); fprintf('\n');
fprintf('  Time knots (s)     : '); fprintf('%.1f  ',  p.wp_times); fprintf('\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('  GA time        : %.0f s (one run)\n', info.t_p1);
fprintf('  Total         : %.0f s (%.1f min)\n', info.total_time, info.total_time/60);
fprintf('==========================================================================\n\n');

if ~traj.feasible
    fprintf('  *** DIAGNOSTIC ONLY: no verified feasible trajectory was found. ***\n\n');
end

%% ---- PLOTS ---------------------------------------------------------------

results_arr(1).launch_angle_deg = LAUNCH_ANGLE_DEG;
results_arr(1).range_km         = traj.x_final / 1e3;
results_arr(1).feasible         = traj.feasible;
results_arr(1).traj             = traj;
results_arr(1).info             = info;

plot_results(traj, results_arr, p);

%% ---- SAVE ----------------------------------------------------------------

if SAVE_RESULTS
    tag     = sprintf('ang%d', LAUNCH_ANGLE_DEG);
    outfile = fullfile(RESULTS_DIR, ['results_' tag '.mat']);
    save(outfile, 'traj', 'wp_opt', 'p', 'info');

    T_exp = table(traj.t, traj.x/1e3, rad2deg(traj.central_angle_rad), ...
                  traj.h/1e3, traj.V, ...
                  traj.M, traj.n, traj.alpha, traj.CMy, rad2deg(traj.gamma), ...
                  'VariableNames', ...
                  {'t_s','x_km','central_angle_deg','h_km','V_ms', ...
                   'Mach','n_g','alpha_deg', ...
                   'CMy_source','gamma_deg'});
    writetable(T_exp, fullfile(RESULTS_DIR, ['trajectory_' tag '.csv']));
    fprintf('Saved: %s\n', outfile);
end
