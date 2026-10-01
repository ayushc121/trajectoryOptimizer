function diagnose_fmincon_trajectory(result_file)
% DIAGNOSE_FMINCON_TRAJECTORY  Compare the verified SQP seed and final flight.
% In MATLAB: diagnose_fmincon_trajectory
% Or: diagnose_fmincon_trajectory('results/results_fmincon_start40_h1800_launchfree.mat')
% Requires the saved .mat file and the existing simulation functions.

folder = fileparts(mfilename('fullpath'));
addpath(folder);
results_dir = fullfile(folder,'results');
if nargin < 1 || isempty(result_file)
    files = dir(fullfile(results_dir,'results_fmincon*.mat'));
    if isempty(files)
        error('diagnose_fmincon_trajectory:noFile', ...
            'No fmincon result found in %s. Supply a .mat path.',results_dir);
    end
    [~,idx] = max([files.datenum]);
    result_file = fullfile(files(idx).folder,files(idx).name);
end
saved = load(result_file,'p','traj','wp_opt','info');
if ~all(isfield(saved,{'p','traj','wp_opt','info'})) || ...
        ~isfield(saved.info,'start_alpha') || ...
        ~isfield(saved.info,'start_launch_angle_deg')
    error('diagnose_fmincon_trajectory:fields', ...
        'The result needs p, traj, wp_opt, info.start_alpha, and info.start_launch_angle_deg.');
end
p = saved.p;
if numel(saved.wp_opt) ~= numel(p.wp_times) || ...
        numel(saved.info.start_alpha) ~= numel(p.wp_times)
    error('diagnose_fmincon_trajectory:grid', ...
        'The saved AoA schedules do not match the saved time grid.');
end

p_start = p;
p_start.gamma0 = deg2rad(saved.info.start_launch_angle_deg);
start_traj = simulate_trajectory(saved.info.start_alpha,p_start,true);
final_traj = simulate_trajectory(saved.wp_opt,p,true);
if abs(final_traj.x_final-saved.traj.x_final)> ...
        max(100,0.001*saved.traj.x_final) || ...
        final_traj.feasible ~= saved.traj.feasible
    warning('diagnose_fmincon_trajectory:replay', ...
        'The final replay differs from the saved trajectory. Check the active simulation files.');
end

fprintf('\nSource: %s\n',result_file);
fprintf('                Launch     Range    Flight   Peak h  Re-climb  Impact M   Feasible\n');
report('SQP start',start_traj,p_start);
report('Final',final_traj,p);
fprintf('Re-climb sums all positive altitude increments after the first maximum.\n');
fprintf('Impact constraint: V >= %.2f m/s, equivalent to Mach %.2f at sea level.\n', ...
    p.V_min_impact,p.M_min_impact);
fprintf('Final impact: V = %.2f m/s; Mach = %.3f; speed margin = %.2f m/s.\n', ...
    final_traj.V_final,final_traj.M_final, ...
    final_traj.V_final-p.V_min_impact);

fig = figure('Name','SQP start versus final flight','Position',[90 90 1050 790]);
subplot(3,1,1); hold on; grid on;
plot(start_traj.x/1e3,start_traj.h/1e3,'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.7);
plot(final_traj.x/1e3,final_traj.h/1e3,'b-','LineWidth',1.7);
yline(p.h_max/1e3,'r--','Altitude limit');
xlabel('Surface-arc range [km]'); ylabel('Altitude [km]');
legend('Verified SQP start','Final','Altitude limit','Location','best');
title('Did the SQP start already contain the altitude waves?');

subplot(3,1,2); hold on; grid on;
plot(start_traj.t,start_traj.h/1e3,'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.5);
plot(final_traj.t,final_traj.h/1e3,'b-','LineWidth',1.5);
yline(p.h_max/1e3,'r--');
xlabel('Time [s]'); ylabel('Altitude [km]');
legend('Verified SQP start','Final','Location','best');

subplot(3,1,3); hold on; grid on;
plot(start_traj.t,start_traj.alpha,'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.3);
plot(final_traj.t,final_traj.alpha,'b-','LineWidth',1.3);
xlabel('Time [s]'); ylabel('AoA command [deg]');
legend('Verified SQP start','Final','Location','best');

out = fullfile(results_dir,'fmincon_seed_vs_final.png');
saveas(fig,out);
fprintf('Saved comparison: %s\n',out);
end

function report(label,tr,p)
[~,first_peak] = max(tr.h);
reclimb_km = sum(max(0,diff(tr.h(first_peak:end))))/1e3;
fprintf('%-10s %8.3f %9.2f %9.1f %8.2f %9.2f %9.3f %9d\n', ...
    label,rad2deg(p.gamma0),tr.x_final/1e3,tr.t_flight, ...
    tr.h_max/1e3,reclimb_km,tr.M_final,tr.feasible);
end
