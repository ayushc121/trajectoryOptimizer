%% MAIN_TRAJECTORY_OPTIMIZER.M
% ==========================================================================
% Hypersonic Maneuvering Projectile — Trajectory Optimiser
% FY26 UCAH Hypersonic Design Challenge
% ==========================================================================
%
% WORKFLOW:
%   Phase 1 — Exploratory GA (22 vars: alpha + waypoint fracs).
%   Phase 2 — Polish GA (11 vars: alpha only, waypoints fixed from Phase 1).
%             Seeded from Phase 1 final population. Higher penalty pressure.
%   Phase 3 — fmincon SQP from Phase 2 feasible result (set N_RANDOM_STARTS=0
%             to skip entirely — recommended until Phase 2 proves reliable).
%
% CRITICAL: hybrid_optimize modifies p.wp_ranges internally (local copy).
%   main must apply info.wp_ranges_opt to p before calling simulate_trajectory.
%   This is already done below — do not remove that line.

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));

%% ---- USER SETTINGS -------------------------------------------------------

LAUNCH_ANGLE_DEG = 40;
N_RANDOM_STARTS  = 0;   % Phase 3 fmincon restarts. Set 0 to skip Phase 3.

SAVE_RESULTS = true;
RESULTS_DIR  = 'results';

%% ---- INITIALISATION ------------------------------------------------------

if SAVE_RESULTS && ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end

p        = build_vehicle_params();
p.gamma0 = deg2rad(LAUNCH_ANGLE_DEG);

AERO_FILE = 'aero_tables.mat';
if ~exist(AERO_FILE, 'file'), generate_aero_tables(AERO_FILE); end
aero_data     = load(AERO_FILE);
p.aero        = aero_data;
p.M_min_table = min(aero_data.Mach_vec);
p.M_max_table = max(aero_data.Mach_vec);

p = setup_waypoints(p);

fprintf('==========================================================================\n');
fprintf('  FY26 UCAH — Hypersonic Maneuvering Projectile Trajectory Optimiser\n');
fprintf('==========================================================================\n');
fprintf('  Launch angle    : %d deg\n',   LAUNCH_ANGLE_DEG);
fprintf('  Phase 3 starts  : %d (%s)\n',  N_RANDOM_STARTS, ...
        ternary(N_RANDOM_STARTS > 0, 'fmincon enabled', 'fmincon disabled'));
fprintf('  Alpha variables : %d  |  Total GA vars (P1): %d\n', p.N_free_wp, 2*p.N_free_wp);
fprintf('  Nominal wp (km) : '); fprintf('%.0f  ', p.wp_ranges/1e3); fprintf('\n');
fprintf('  x_max estimate  : %.0f km  (ballistic x %.2f)\n', p.ball.x_max/1e3, p.range_est_factor);
fprintf('  V_min impact    : %.1f m/s (Mach %.1f)\n', p.V_min_impact, p.M_min_impact);
fprintf('  Phase 1 penalty : w_alt=%.1f  w_vel=%.1f\n', p.ga_w_alt, p.ga_w_vel);
fprintf('  Phase 2 penalty : w_alt=%.1f  w_vel=%.1f  (mult=%.1fx)\n', ...
        p.ga_w_alt * p.ga_polish_penalty_mult, p.ga_w_vel * p.ga_polish_penalty_mult, ...
        p.ga_polish_penalty_mult);
fprintf('==========================================================================\n\n');

%% ---- OPTIMISE ------------------------------------------------------------

[wp_opt, range_opt, info] = hybrid_optimize(p, N_RANDOM_STARTS);

% CRITICAL BUG FIX: hybrid_optimize updates p.wp_ranges internally (MATLAB
% pass-by-value means the outer p was not modified). Apply the GA-optimal
% waypoint positions before any further simulate_trajectory calls.
% Without this line, the final trajectory uses the wrong waypoint positions
% and produces a completely incorrect result.
p.wp_ranges = info.wp_ranges_opt;

fprintf('\n  *** Best result: %.1f km  |  feasible: %d ***\n\n', ...
        range_opt/1e3, simulate_trajectory(wp_opt, p, false).feasible);

%% ---- DETAILED SIMULATION -------------------------------------------------

traj = simulate_trajectory(wp_opt, p, true);

fprintf('==========================================================================\n');
fprintf('  FINAL TRAJECTORY SUMMARY\n');
fprintf('==========================================================================\n');
fprintf('  Range          : %.1f km\n',   traj.x_final/1e3);
fprintf('  Max altitude   : %.1f km\n',   traj.h_max/1e3);
fprintf('  Flight time    : %.1f s\n',    traj.t_flight);
fprintf('  Terminal Mach  : %.3f\n',      traj.M_final);
fprintf('  Terminal speed : %.1f m/s  (min %.1f)\n', traj.V_final, p.V_min_impact);
fprintf('  Max g-load     : %.2f g\n',    traj.n_max);
fprintf('  Feasible       : %d\n',        traj.feasible);
fprintf('  Constraints (c<=0 satisfied):\n');
fprintf('    Altitude (m) : %+.1f\n',   traj.c_viol(1));
fprintf('    G-load   (g) : %+.3f\n',   traj.c_viol(2));
fprintf('    Speed  (m/s) : %+.1f\n',  -traj.c_viol(3));
fprintf('--------------------------------------------------------------------------\n');
fprintf('  Phase 1 alpha [deg]: '); fprintf('%+.1f  ', info.wp_ga);     fprintf('\n');
fprintf('  Phase 2 alpha [deg]: '); fprintf('%+.1f  ', info.wp_polish); fprintf('\n');
fprintf('  Final   alpha [deg]: '); fprintf('%+.1f  ', wp_opt);         fprintf('\n');
fprintf('  Waypoints (km)     : '); fprintf('%.0f  ',  p.wp_ranges/1e3); fprintf('\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('  Phase 1 time  : %.0f s\n', info.t_p1);
fprintf('  Phase 2 time  : %.0f s\n', info.t_p2);
if info.phase3_ran
fprintf('  Phase 3 time  : %.0f s\n', sum(info.all_times));
end
fprintf('  Total         : %.0f s (%.1f min)\n', info.total_time, info.total_time/60);
fprintf('==========================================================================\n\n');

if ~traj.feasible
    fprintf('  *** INFEASIBLE — Remedies:\n');
    if traj.c_viol(1) > 0
        fprintf('    Altitude: raise ga_w_alt (%.1f) or ga_polish_penalty_mult (%.1f)\n', ...
                p.ga_w_alt, p.ga_polish_penalty_mult);
    end
    if traj.c_viol(3) > 0
        fprintf('    Speed:    raise ga_w_vel (%.1f) or ga_polish_vel_margin (%.1f m/s)\n', ...
                p.ga_w_vel, p.ga_polish_vel_margin);
    end
    fprintf('\n');
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

    T_exp = table(traj.t, traj.x/1e3, traj.h/1e3, traj.V, ...
                  traj.M, traj.n, traj.alpha, rad2deg(traj.gamma), ...
                  'VariableNames', ...
                  {'t_s','x_km','h_km','V_ms','Mach','n_g','alpha_deg','gamma_deg'});
    writetable(T_exp, fullfile(RESULTS_DIR, ['trajectory_' tag '.csv']));
    fprintf('Saved: %s\n', outfile);
end


%% ---- Helper ---------------------------------------------------------------
function s = ternary(cond, a, b)
if cond; s = a; else; s = b; end
end

