function p = setup_waypoints(p)
% SETUP_WAYPOINTS  Build the single non-uniform waypoint grid from a
%                  ballistic reference simulation.
%
% Procedure:
%   1. Run a zero-alpha (ballistic) trajectory. Result is independent of
%      placeholder wp_ranges since alpha=0 regardless of interpolation.
%   2. Set x_max = ballistic_range * range_est_factor.
%   3. Place N_wp nominal waypoints at wp_fracs * x_max.
%   4. Compute per-waypoint fraction bounds for GA joint optimisation.
%      Each waypoint can shift +-ga_frac_delta around its nominal fraction.
%   5. Set T_max with margin over ballistic flight time.
%   6. Store ballistic reference in p.ball for penalty scaling and plotting.
%
% NOTE on GA waypoint optimisation:
%   p.frac_lb and p.frac_ub bound the fraction variables the GA will
%   jointly optimise alongside alpha. After Phase 1, the GA-optimal
%   fractions overwrite p.wp_ranges, so Phases 2 and 3 operate at the
%   positions the GA found to be best, not the nominal positions set here.

fprintf('\n[setup_waypoints] Running ballistic reference trajectory...\n');

%% --- Step 1: Ballistic reference ---
alpha_zero = zeros(1, p.N_wp);
traj_ball  = simulate_trajectory(alpha_zero, p, false);

fprintf('[setup_waypoints] Ballistic: range = %.1f km  |  h_max = %.1f km  |  t = %.1f s\n', ...
        traj_ball.x_final/1e3, traj_ball.h_max/1e3, traj_ball.t_flight);

%% --- Step 2: x_max ---
x_max = traj_ball.x_final * p.range_est_factor;
fprintf('[setup_waypoints] x_max (est) = %.1f km  (%.1f km x %.2f)\n', ...
        x_max/1e3, traj_ball.x_final/1e3, p.range_est_factor);

%% --- Step 3: Nominal waypoint positions ---
p.wp_ranges = p.wp_fracs * x_max;
p.wp_lb     = -p.alpha_max_deg * ones(1, p.N_wp);
p.wp_ub     =  p.alpha_max_deg * ones(1, p.N_wp);

min_spacing = min(diff(p.wp_ranges));
if min_spacing < 2000
    warning('setup_waypoints: minimum waypoint spacing %.1f km — consider reducing range_est_factor.', ...
            min_spacing/1e3);
end

%% --- Step 4: Fraction bounds for GA joint optimisation ---
% Each waypoint fraction is allowed to shift by +-ga_frac_delta from its
% nominal position, clamped to [0.005, 0.995] to avoid extreme positions.
% Adjacent fraction bounds can overlap — the joint GA objective function
% sorts the fractions internally, so any permutation yields the same
% trajectory and the GA naturally avoids wasting DOF on reordering.
p.frac_lb = max(0.005, p.wp_fracs - p.ga_frac_delta);
p.frac_ub = min(0.995, p.wp_fracs + p.ga_frac_delta);

fprintf('[setup_waypoints] Fraction bounds (each waypoint +-%.2f of x_max):\n', p.ga_frac_delta);
fprintf('  Nominal fracs: '); fprintf('%.2f  ', p.wp_fracs);   fprintf('\n');
fprintf('  lb fracs     : '); fprintf('%.2f  ', p.frac_lb);    fprintf('\n');
fprintf('  ub fracs     : '); fprintf('%.2f  ', p.frac_ub);    fprintf('\n');

%% --- Step 5: T_max ---
p.T_max = 2.5 * traj_ball.t_flight;

%% --- Step 6: Ballistic reference struct ---
% p.ball.traj stores the full trajectory for plot_results.m.
% p.ball.x_max is used in ga_objective penalty scaling.
p.ball.x_max    = x_max;
p.ball.x_ball   = traj_ball.x_final;
p.ball.h_max    = traj_ball.h_max;
p.ball.t_flight = traj_ball.t_flight;
p.ball.traj     = traj_ball;   % full struct for plotting (bx, bh, etc.)

fprintf('[setup_waypoints] Nominal waypoints (km): ');
fprintf('%.1f  ', p.wp_ranges/1e3); fprintf('\n');
fprintf('[setup_waypoints] T_max = %.0f s\n\n', p.T_max);

end
