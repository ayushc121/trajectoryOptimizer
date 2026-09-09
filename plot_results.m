function plot_results(best_traj, results, p)
% PLOT_RESULTS  Trajectory optimisation visualisation
%
% Figure 1 — Altitude vs downrange (trajectory + ballistic reference + alpha schedule)
% Figure 2 — Alpha schedule (own panel — this is the primary control output)
% Figure 3 — Six-panel flight data vs downrange
% Figure 4 — Launch angle comparison bar chart (only when results has >1 entry)
%
% CHANGES vs previous version:
%   - p.ball.x / p.ball.h replaced with p.ball.traj.x / p.ball.traj.h
%     (setup_waypoints now stores the full ballistic traj in p.ball.traj)
%   - best_traj.gamma0 replaced with p.gamma0 (more robust; gamma0 is
%     always available in p but may or may not be stored in traj struct)
%   - Load factor limit label now uses p.n_max instead of hardcoded 15

t   = best_traj.t;
x   = best_traj.x / 1e3;          % [km]
h   = best_traj.h / 1e3;          % [km]
V   = best_traj.V;
gam = rad2deg(best_traj.gamma);
M   = best_traj.M;
n   = best_traj.n;
al  = best_traj.alpha;             % [deg] per-timestep alpha from schedule

x_end = best_traj.x_final / 1e3;
x_ax  = x_end * 1.08;

% Ballistic reference — stored in p.ball.traj by setup_waypoints
bx = p.ball.traj.x / 1e3;
bh = p.ball.traj.h / 1e3;

% Alpha schedule knots (piecewise-linear control profile)
wp_x_km  = p.wp_ranges / 1e3;
wp_alpha = best_traj.wp_alphas;    % [deg] at each waypoint (set by simulate_trajectory)

%% ==============================================================
%% Figure 1: Trajectory profile
%% ==============================================================
figure('Name','Trajectory Profile','Position',[80 80 1000 500]);

% Forbidden zone above 30 km ceiling
fill([0 x_ax x_ax 0], [p.h_max/1e3 p.h_max/1e3 40 40], ...
     [1 0.78 0.78], 'EdgeColor','none','FaceAlpha',0.5);
hold on; grid on; box on;

% Ceiling line
yline(p.h_max/1e3, 'r--', 'LineWidth', 1.5);

% Ballistic reference
plot(bx, bh, 'Color',[0.6 0.6 0.6], 'LineStyle','--', 'LineWidth', 1.5);

% Optimised trajectory
plot(x, h, 'b-', 'LineWidth', 2.2);

% Waypoint x-positions — vertical ticks with alpha annotation
y_tick = 0.5;
for k = 1:length(wp_x_km)
    xpos = wp_x_km(k);
    if xpos <= x_ax
        xline(xpos, ':', 'Color', [0.4 0.4 0.4], 'Alpha', 0.5);
        text(xpos, y_tick, sprintf('%.0f°', wp_alpha(k)), ...
             'HorizontalAlignment','center', 'FontSize', 7, ...
             'Color', [0.2 0.2 0.6]);
    end
end

% Launch and impact markers
plot(0,     0, 'gs', 'MarkerSize', 9, 'MarkerFaceColor','g');
plot(x_end, 0, 'rv', 'MarkerSize', 9, 'MarkerFaceColor','r');

xlabel('Downrange [km]', 'FontSize', 12);
ylabel('Altitude [km]',  'FontSize', 12);
title(sprintf('Trajectory  |  \\gamma_0 = %.0f deg  |  Range = %.1f km  |  h_{max} = %.1f km  |  M_{final} = %.2f  |  Feasible: %d', ...
      rad2deg(p.gamma0), x_end, best_traj.h_max/1e3, ...
      best_traj.M_final, best_traj.feasible), 'FontSize', 11);
legend({'Forbidden zone','30 km ceiling','Ballistic ref.','Trajectory','Launch','Impact'}, ...
       'Location','northeast','FontSize', 9);
xlim([0, x_ax]);
ylim([0, max(36, best_traj.h_max/1e3 * 1.1)]);

%% ==============================================================
%% Figure 2: Alpha schedule
%% ==============================================================
figure('Name','Alpha Schedule','Position',[90 90 800 280]);
hold on; grid on; box on;

fill([0 x_ax x_ax 0], [p.alpha_max_deg  p.alpha_max_deg ...
                       -p.alpha_max_deg -p.alpha_max_deg], ...
     [0.9 0.95 1.0], 'EdgeColor','none','FaceAlpha',0.6);

plot(x, al, 'b-', 'LineWidth', 1.5);
plot(wp_x_km, wp_alpha, 'ro-', 'MarkerSize', 7, 'MarkerFaceColor','r', 'LineWidth', 1.5);

yline( p.alpha_max_deg, 'r--', 'LineWidth', 1);
yline(-p.alpha_max_deg, 'r--', 'LineWidth', 1);
yline(0, 'k-', 'LineWidth', 0.8, 'Alpha', 0.4);

xlabel('Downrange [km]', 'FontSize', 12);
ylabel('\alpha [deg]',   'FontSize', 12);
title('AoA Schedule — trajectory (blue) and waypoint knots (red)', 'FontSize', 11);
xlim([0, x_ax]);
ylim([-p.alpha_max_deg*1.3, p.alpha_max_deg*1.3]);
legend({'±\alpha_{max} band','Schedule (trajectory)','Waypoint knots'}, ...
       'Location','best','FontSize',9);

%% ==============================================================
%% Figure 3: Six-panel flight data
%% ==============================================================
figure('Name','Flight Data','Position',[130 130 1200 700]);

% 3a — Mach
ax1 = subplot(2,3,1);
plot(x, M, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline(p.M_min_impact, 'r--', 'LineWidth', 1.5);
text(x_end*0.05, p.M_min_impact*1.08, sprintf('M = %.0f limit', p.M_min_impact), ...
     'Color','r','FontSize',8);
xlabel('Downrange [km]'); ylabel('Mach [-]'); title('Mach Number');
ylim([0, p.M_launch * 1.05]);

% 3b — Load factor
ax2 = subplot(2,3,2);
plot(x, n, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline( p.n_max, 'r--', 'LineWidth', 1.5);
yline(-p.n_max, 'r--', 'LineWidth', 1.5);
% Use p.n_max dynamically rather than hardcoding ±15
text(x_end*0.05, p.n_max*0.85, sprintf('\\pm%.0f g limit', p.n_max), ...
     'Color','r','FontSize',8);
xlabel('Downrange [km]'); ylabel('n [g]'); title('Load Factor');

% 3c — Speed
ax3 = subplot(2,3,3);
plot(x, V, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline(p.V_min_impact, 'r--', 'LineWidth', 1.5);
text(x_end*0.05, p.V_min_impact*1.05, 'V_{min}', 'Color','r','FontSize',8);
xlabel('Downrange [km]'); ylabel('Speed [m/s]'); title('Airspeed');

% 3d — Flight path angle
ax4 = subplot(2,3,4);
plot(x, gam, 'b-', 'LineWidth', 1.5); hold on; grid on;
yline(0, 'k-', 'LineWidth', 0.8, 'Alpha', 0.4);
xlabel('Downrange [km]'); ylabel('\gamma [deg]'); title('Flight Path Angle');

% 3e — Altitude
ax5 = subplot(2,3,5);
plot(bx, bh, '--', 'Color',[0.6 0.6 0.6], 'LineWidth', 1.2); hold on; grid on;
plot(x, h, 'b-', 'LineWidth', 1.5);
yline(p.h_max/1e3, 'r--', 'LineWidth', 1.5);
xlabel('Downrange [km]'); ylabel('Altitude [km]'); title('Altitude');
legend({'Ballistic ref.','Trajectory'}, 'Location','best','FontSize',8);
ylim([0, max(35, best_traj.h_max/1e3 * 1.1)]);

% 3f — L/D ratio
ax6 = subplot(2,3,6);
if any(M > 0)
    LD = zeros(size(t));
    for i = 1:length(t)
        if M(i) > 0
            [CLi, CDi] = aero_lookup(M(i), al(i), p.aero);
            if CDi > 1e-8, LD(i) = CLi / CDi; end
        end
    end
    plot(x, LD, 'b-', 'LineWidth', 1.5); grid on;
    yline(0, 'k-', 'Alpha', 0.4);
else
    text(0.5, 0.5, 'Run with detailed=true for L/D', ...
         'Units','normalized','HorizontalAlignment','center','FontSize',10);
    grid on;
end
xlabel('Downrange [km]'); ylabel('L/D [-]'); title('Lift-to-Drag Ratio');

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x');
for ax = [ax1 ax2 ax3 ax4 ax5 ax6]
    xlim(ax, [0, x_ax]);
end

%% ==============================================================
%% Figure 4: Launch angle comparison (only with multiple results)
%% ==============================================================
if length(results) > 1
    figure('Name','Launch Angle Sweep','Position',[200 200 620 380]);
    la_all   = [results.launch_angle_deg];
    rng_all  = [results.range_km];
    feas_all = logical([results.feasible]);

    bar(la_all( feas_all), rng_all( feas_all), 'b'); hold on;
    bar(la_all(~feas_all), rng_all(~feas_all), 'r');
    grid on;
    xlabel('Launch Angle [deg]', 'FontSize', 12);
    ylabel('Range [km]',         'FontSize', 12);
    title('Range vs Launch Angle', 'FontSize', 12);
    legend({'Feasible','Infeasible'}, 'Location','best');
end

end
