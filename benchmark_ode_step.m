% BENCHMARK_ODE_STEP  Compare ODE maximum steps using the time control grid.
clear; clc;
folder = fileparts(mfilename('fullpath'));
addpath(folder);
p = build_vehicle_params();
p.gamma0 = deg2rad(40);
aero_file = fullfile(folder,'step_benchmark_aero.mat');
generate_aero_tables(aero_file, ...
    fullfile(folder,'owen_tricon_v1_data.csv'),p.S_ref,p.L_ref);
p.aero = load(aero_file);
p.M_min_table = min(p.aero.Mach_vec);
p.M_max_table = max(p.aero.Mach_vec);
p = setup_waypoints(p);

profiles(1).name = 'Zero lift';
profiles(1).times = p.wp_times;
profiles(1).alpha = p.aero.alpha_zero_lift_deg*ones(size(p.wp_times));
profiles(1).gamma0 = 40;
profiles(2).name = 'Low-angle complete-flight probe';
profiles(2).times = p.wp_times;
profiles(2).alpha = profiles(1).alpha;
profiles(2).gamma0 = 5;
profiles(3).name = 'Constant best L/D';
profiles(3).times = p.wp_times;
profiles(3).alpha = p.reference.alpha_max_ld_deg*ones(size(p.wp_times));
profiles(3).gamma0 = 40;

saved_file = fullfile(pwd,'results','results_ang40.mat');
if exist(saved_file,'file')
    saved = load(saved_file,'wp_opt','p');
    if isfield(saved,'wp_opt') && isfield(saved,'p') && ...
            isfield(saved.p,'wp_times') && ...
            numel(saved.wp_opt)==numel(saved.p.wp_times)
        profiles(end+1).name = 'Saved time-based control schedule';
        profiles(end).times = saved.p.wp_times;
        profiles(end).alpha = saved.wp_opt;
        profiles(end).gamma0 = rad2deg(saved.p.gamma0);
    end
end

steps = [2,1,0.5,0.25];
for j=1:numel(profiles)
    p.wp_times = profiles(j).times;
    fprintf('\n%s\n',profiles(j).name);
    for k=1:numel(steps)
        pk = p;
        pk.gamma0 = deg2rad(profiles(j).gamma0);
        pk.ode_max_step = steps(k);
        timer = tic;
        tr(k) = simulate_trajectory(profiles(j).alpha,pk,false); %#ok<SAGROW>
        seconds(k) = toc(timer); %#ok<SAGROW>
    end
    reference = tr(end);
    fprintf(' step(s)  time(s)  range(km)  hmax(km)  dh(m)  dx(m)  dn(g)  feasible  landed  stop\n');
    for k=1:numel(steps)
        if tr(k).landed
            stop_reason = 'impact';
        elseif tr(k).h_max >= p.abort_altitude_m-1
            stop_reason = 'altitude';
        elseif tr(k).V_final <= 1.01
            stop_reason = 'speed';
        else
            stop_reason = 'time';
        end
        fprintf(' %6.2f  %7.2f  %9.2f  %8.2f  %+6.1f  %+6.1f  %+6.3f  %8d  %6d  %s\n', ...
            steps(k),seconds(k),tr(k).x_final/1e3,tr(k).h_max/1e3, ...
            tr(k).h_max-reference.h_max,tr(k).x_final-reference.x_final, ...
            tr(k).n_max-reference.n_max,tr(k).feasible,tr(k).landed,stop_reason);
    end
    clear tr seconds
end
delete(aero_file);
fprintf('\nThe 0.25 s row is the reference. Incomplete flights test only the path up to the stop event.\n');
