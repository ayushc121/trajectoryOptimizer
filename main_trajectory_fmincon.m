function main_trajectory_fmincon()
% MAIN_TRAJECTORY_FMINCON  Single SQP search from a best-L/D cruise seed.
% All dynamics and impact constraints come from the existing simulator.
% Requires MATLAB Optimization Toolbox.

%% User settings
LAUNCH_ANGLE_DEG = 40;
LAUNCH_ANGLE_MIN_DEG = 20;
LAUNCH_ANGLE_MAX_DEG = 60;
CONTROL_HORIZON_S = 1800;
MAX_FLIGHT_TIME_S = 2500;
BASE_CONTROL_HORIZON_S = 1200;
N_BASE_TIME_KNOTS = 180;
N_EXTRA_TIME_KNOTS = 60;
N_EARLY_KNOTS = 30;
EARLY_WINDOW_S = 60;
SUBSTEPS_PER_OLD_INTERVAL = 3;
N_OPT_EARLY = 10;
N_OPT_MID = 24;
N_OPT_LATE = 8;
MAX_CORRECTION_DEG = 5;
CRUISE_AOA_WINDOW_DEG = 0.5; % allowed departure from best-L/D AoA during cruise
EQUIL_ENTRY_TOLERANCE_KM = 2; % prefer a feasible seed near its entry equilibrium
MAX_ITER = 80;
MAX_EVAL = 2500;
FD_STEP_DEG = 0.03;
LAUNCH_FD_STEP_DEG = 0.05;
RNG_SEED = 42;

folder=fileparts(mfilename('fullpath'));
addpath(folder);
results_dir=fullfile(folder,'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end
if exist('fmincon','file')~=2
    error('main_trajectory_fmincon:toolbox','Optimization Toolbox is required.');
end
fprintf('Running: %s\n',mfilename('fullpath'));
fprintf('Dependencies: %s | %s\n',which('build_vehicle_params'), ...
    which('ld_equilibrium_seeds'));
rng(RNG_SEED,'twister');
p=build_vehicle_params();
p.gamma0=deg2rad(LAUNCH_ANGLE_DEG);
p.n_time_knots=SUBSTEPS_PER_OLD_INTERVAL* ...
    (N_BASE_TIME_KNOTS+N_EXTRA_TIME_KNOTS-1)+1;
p.n_early_knots=SUBSTEPS_PER_OLD_INTERVAL*(N_EARLY_KNOTS-1)+1;
p.early_window_s=EARLY_WINDOW_S;
p.time_horizon_override_s=CONTROL_HORIZON_S;
p.T_max=max(p.T_max,MAX_FLIGHT_TIME_S);
aero_file=fullfile(folder,'aero_tables.mat');
generate_aero_tables(aero_file, ...
    fullfile(folder,'owen_tricon_v1_data.csv'),p.S_ref,p.L_ref);
p.aero=load(aero_file);
p.M_min_table=min(p.aero.Mach_vec);
p.M_max_table=max(p.aero.Mach_vec);
p=setup_waypoints(p);
if BASE_CONTROL_HORIZON_S<=EARLY_WINDOW_S || ...
        BASE_CONTROL_HORIZON_S>=CONTROL_HORIZON_S || ...
        N_BASE_TIME_KNOTS<=N_EARLY_KNOTS || N_EXTRA_TIME_KNOTS<1 || ...
        SUBSTEPS_PER_OLD_INTERVAL<1 || ...
        SUBSTEPS_PER_OLD_INTERVAL~=floor(SUBSTEPS_PER_OLD_INTERVAL)
    error('main_trajectory_fmincon:grid','Invalid control-grid settings.');
end
old_t=[linspace(0,EARLY_WINDOW_S,N_EARLY_KNOTS), ...
    EARLY_WINDOW_S+(1:(N_BASE_TIME_KNOTS-N_EARLY_KNOTS))* ...
    (BASE_CONTROL_HORIZON_S-EARLY_WINDOW_S)/ ...
    (N_BASE_TIME_KNOTS-N_EARLY_KNOTS), ...
    BASE_CONTROL_HORIZON_S+(1:N_EXTRA_TIME_KNOTS)* ...
    (CONTROL_HORIZON_S-BASE_CONTROL_HORIZON_S)/N_EXTRA_TIME_KNOTS];
dt=diff(old_t);
segments=zeros(SUBSTEPS_PER_OLD_INTERVAL,numel(dt));
for k=0:SUBSTEPS_PER_OLD_INTERVAL-1
    segments(k+1,:)=old_t(1:end-1)+k*dt/SUBSTEPS_PER_OLD_INTERVAL;
end
p.wp_times=[reshape(segments,1,[]),old_t(end)];
if numel(p.wp_times)~=p.N_wp || ...
        abs(p.wp_times(end)-CONTROL_HORIZON_S)>1e-8 || ...
        p.T_max<=p.wp_times(end)
    error('main_trajectory_fmincon:grid','Time grid and horizon do not match.');
end
anchor_t=[linspace(0,EARLY_WINDOW_S,N_OPT_EARLY), ...
    EARLY_WINDOW_S+(1:N_OPT_MID)* ...
    (BASE_CONTROL_HORIZON_S-EARLY_WINDOW_S)/N_OPT_MID, ...
    BASE_CONTROL_HORIZON_S+(1:N_OPT_LATE)* ...
    (CONTROL_HORIZON_S-BASE_CONTROL_HORIZON_S)/N_OPT_LATE];
n_control=numel(anchor_t);
B=interp1(anchor_t,eye(n_control),p.wp_times,'linear');
if any(~isfinite(B(:))) || size(B,1)~=p.N_wp
    error('main_trajectory_fmincon:basis','Invalid AoA correction basis.');
end
fprintf('Single max-L/D search: %d AoA knots, %d correction anchors and launch angle; horizon %.0f s, T_max %.0f s.\n', ...
    p.N_wp,n_control,p.wp_times(end),p.T_max);
fprintf('Best-L/D AoA %.3f deg. Instantaneous level-flight balance:\n', ...
    p.reference.alpha_max_ld_deg);
for mach=[8,7,6,5,4,3]
    [h_eq,exists]=ld_equilibrium_altitude(p,mach,'mach');
    if exists
        fprintf('  Mach %.0f: %.2f km\n',mach,h_eq/1e3);
    else
        fprintf('  Mach %.0f: outside 0--%.0f km\n',mach,p.h_max/1e3);
    end
end

%% Screen max-L/D cruise histories, including a near-ceiling naive history.
[seeds,labels,design]=ld_equilibrium_seeds(p);
range=-inf(size(seeds,1),1);
trials=cell(size(seeds,1),1);
entry_error=inf(size(range));
for k=1:size(seeds,1)
    trials{k}=simulate_trajectory(seeds(k,:),p,false);
    if trials{k}.feasible
        range(k)=trials{k}.x_final;
        entry_error(k)=entry_mismatch(trials{k},design(k).switch_s,p);
    end
end
fprintf('Max-L/D seeds: %d/%d feasible; %d within %.1f km of balance altitude at cruise entry.\n', ...
    sum(isfinite(range)),numel(range), ...
    sum(isfinite(range)&entry_error<=EQUIL_ENTRY_TOLERANCE_KM), ...
    EQUIL_ENTRY_TOLERANCE_KM);
valid=find(isfinite(range));
if isempty(valid)
    error('main_trajectory_fmincon:noSeed', ...
        'No max-L/D seed landed with altitude/load/impact-speed constraints.');
end
near=valid(entry_error(valid)<=EQUIL_ENTRY_TOLERANCE_KM);
[~,rank_near]=sort(range(near),'descend');
[~,rank_all]=sort(entry_error(valid),'ascend');
order=unique([near(rank_near);valid(rank_all)],'stable');
seed_idx=NaN;
for k=order(:)'
    if replay_ok(seeds(k,:),p,trials{k})
        seed_idx=k; break;
    end
end
if isnan(seed_idx)
    error('main_trajectory_fmincon:replay', ...
        'No max-L/D seed passed the finer ODE replay.');
end
seed_alpha=seeds(seed_idx,:);
seed_traj=trials{seed_idx};
switch_s=design(seed_idx).switch_s;
descend_s=sscanf(labels{seed_idx},'ld_turn%d_switch%d_down%d');
descend_s=descend_s(end);
fprintf('Selected equilibrium-entry seed: %s; %.2f km, peak %.2f km, entry mismatch %.2f km, Mach %.3f.\n', ...
    labels{seed_idx},seed_traj.x_final/1e3,seed_traj.h_max/1e3, ...
    entry_error(seed_idx),seed_traj.M_final);
[~,entry_eq,entry_h,entry_V,entry_gamma]=entry_mismatch( ...
    seed_traj,switch_s,p);
fprintf('At max-L/D switch %.0f s: altitude %.2f km, balance altitude %.2f km, speed %.1f m/s, flight-path angle %.2f deg.\n', ...
    switch_s,entry_h/1e3,entry_eq/1e3,entry_V, ...
    rad2deg(entry_gamma));
peak_heights=cellfun(@(tr)tr.h_max,trials(valid));
near_ceiling=valid(peak_heights>=0.90*p.h_max);
naive_traj=[];
naive_label='';
if ~isempty(near_ceiling)
    [~,naive_order]=sort(range(near_ceiling),'descend');
    for k=near_ceiling(naive_order(:))'
        if replay_ok(seeds(k,:),p,trials{k})
            naive_traj=simulate_trajectory(seeds(k,:),p,true);
            naive_label=labels{k};
            fprintf('Near-ceiling naive max-L/D comparison: %s; %.2f km, peak %.2f km.\n', ...
                naive_label,naive_traj.x_final/1e3,naive_traj.h_max/1e3);
            break;
        end
    end
end

%% Maximize verified impact range, retaining a narrow best-L/D cruise AoA.
cruise_anchor=(anchor_t>=switch_s+30 & anchor_t<=descend_s-30);
max_change=MAX_CORRECTION_DEG*ones(1,n_control);
max_change(cruise_anchor)=CRUISE_AOA_WINDOW_DEG;
seed_x=[zeros(1,n_control),LAUNCH_ANGLE_DEG];
best_x=seed_x;
best_traj=seed_traj;
last_x=[]; last_traj=[];
eval_count=0;
checkpoint_range=seed_traj.x_final;
checkpoint_file=fullfile(results_dir,sprintf( ...
    'fmincon_checkpoint_start%d_h%d_ldcruise.mat', ...
    LAUNCH_ANGLE_DEG,CONTROL_HORIZON_S));
A=[B,zeros(p.N_wp,1);-B,zeros(p.N_wp,1)];
b=[p.wp_ub(:)-seed_alpha(:);seed_alpha(:)-p.wp_lb(:)];
fd_steps=[FD_STEP_DEG*ones(1,n_control), ...
    LAUNCH_FD_STEP_DEG/LAUNCH_ANGLE_DEG];
options=optimoptions('fmincon','Algorithm','sqp','Display','iter', ...
    'MaxIterations',MAX_ITER,'MaxFunctionEvaluations',MAX_EVAL, ...
    'FiniteDifferenceType','forward','FiniteDifferenceStepSize',fd_steps, ...
    'TypicalX',ones(1,n_control+1), ...
    'ConstraintTolerance',1e-3,'StepTolerance',1e-4, ...
    'OptimalityTolerance',1e-3,'OutputFcn',@checkpoint);
tic_opt=tic;
[returned_x,~,exitflag,output]=fmincon(@objective,seed_x,A,b,[],[], ...
    [-max_change,LAUNCH_ANGLE_MIN_DEG], ...
    [max_change,LAUNCH_ANGLE_MAX_DEG],@constraints,options);
seconds=toc(tic_opt);
returned_traj=evaluate(returned_x);
wp_opt=profile(best_x);
p.gamma0=deg2rad(best_x(end));
traj=simulate_trajectory(wp_opt,p,true);
pc=p; pc.ode_max_step=min(0.5,p.ode_max_step/2);
pc.ode_reltol=p.ode_reltol/2;
fine=simulate_trajectory(wp_opt,pc,false);
if ~traj.feasible || ~fine.feasible || ~same_trajectory(traj,fine)
    warning('main_trajectory_fmincon:replay', ...
        'SQP improvement failed finer replay; retaining the max-L/D seed.');
    wp_opt=seed_alpha;
    p.gamma0=deg2rad(LAUNCH_ANGLE_DEG);
    pc.gamma0=p.gamma0;
    traj=simulate_trajectory(wp_opt,p,true);
    fine=simulate_trajectory(wp_opt,pc,false);
    best_x=seed_x;
end
traj.numerically_converged=fine.feasible && same_trajectory(traj,fine);
if ~traj.feasible || ~traj.numerically_converged
    error('main_trajectory_fmincon:unverified', ...
        'Neither the SQP point nor the starting seed passed finer replay.');
end

% Keep the plot reference at the final launch angle.
p_ball=p; p_ball.abort_altitude_m=Inf;
p_ball.control_mode='time_knots';
p.ball.traj=simulate_trajectory( ...
    p.aero.alpha_zero_lift_deg*ones(1,p.N_wp),p_ball,false);
p.ball.x_ball=p.ball.traj.x_final;
p.ball.t_flight=p.ball.traj.t_flight;
info.method='fmincon sqp, best-L/D cruise';
info.start_label=labels{seed_idx};
info.start_alpha=seed_alpha;
info.start_range_m=seed_traj.x_final;
info.start_entry_mismatch_km=entry_error(seed_idx);
info.start_switch_s=switch_s;
info.start_descent_s=descend_s;
info.best_launch_angle_deg=rad2deg(p.gamma0);
info.cruise_aoa_window_deg=CRUISE_AOA_WINDOW_DEG;
info.correction_anchor_times_s=anchor_t;
info.best_correction_deg=best_x(1:n_control);
info.exitflag=exitflag;
info.output=output;
info.returned_feasible=returned_traj.feasible;
info.evaluations=eval_count;
info.seconds=seconds;
info.naive_label=naive_label;
outfile=fullfile(results_dir,sprintf( ...
    'results_fmincon_start%d_h%d_ldcruise.mat', ...
    LAUNCH_ANGLE_DEG,CONTROL_HORIZON_S));
save(outfile,'traj','wp_opt','p','info','seed_traj','naive_traj', ...
    'seeds','labels','design');
T_exp=table(traj.t,traj.x/1e3,traj.h/1e3,traj.V,traj.M, ...
    traj.n,traj.alpha,rad2deg(traj.gamma), ...
    'VariableNames',{'t_s','x_km','h_km','V_ms','Mach','n_g', ...
    'alpha_deg','gamma_deg'});
writetable(T_exp,fullfile(results_dir,sprintf( ...
    'trajectory_fmincon_start%d_h%d_ldcruise.csv', ...
    LAUNCH_ANGLE_DEG,CONTROL_HORIZON_S)));
summary=struct('launch_angle_deg',rad2deg(p.gamma0), ...
    'range_km',traj.x_final/1e3,'feasible',traj.feasible, ...
    'traj',traj,'info',info);
plot_results(traj,summary,p);
fig=figure('Name','Maximum-L/D cruise comparison'); hold on; grid on;
if ~isempty(naive_traj)
    plot(naive_traj.x/1e3,naive_traj.h/1e3, ...
        'Color',[0.5 0.5 0.5],'LineWidth',1.4);
end
seed_detail=simulate_trajectory(seed_alpha,p_with_angle(p,LAUNCH_ANGLE_DEG),true);
plot(seed_detail.x/1e3,seed_detail.h/1e3,'b-','LineWidth',1.5);
plot(traj.x/1e3,traj.h/1e3,'g-','LineWidth',1.8);
yline(p.h_max/1e3,'r--','Altitude ceiling');
xlabel('Surface-arc range [km]');ylabel('Altitude [km]');
if isempty(naive_traj)
    legend('Equilibrium-entry seed','SQP result','Altitude ceiling', ...
        'Location','best');
else
    legend('Near-ceiling naive seed','Equilibrium-entry seed', ...
        'SQP result','Altitude ceiling','Location','best');
end
title(sprintf('Entry mismatch %.2f km | seed %.1f km | SQP %.1f km', ...
    entry_error(seed_idx),seed_traj.x_final/1e3,traj.x_final/1e3));
saveas(fig,fullfile(results_dir,'ld_cruise_comparison.png'));
fprintf('Max-L/D SQP: start %.2f km, verified best %.2f km at %.3f deg; flight %.1f s; %d evaluations in %.1f min.\n', ...
    seed_traj.x_final/1e3,traj.x_final/1e3,rad2deg(p.gamma0), ...
    traj.t_flight,eval_count,seconds/60);
fprintf('Exitflag %d; returned point feasible %d; finer replay %d; impact Mach %.3f; peak %.2f km.\n', ...
    exitflag,returned_traj.feasible,traj.numerically_converged, ...
    traj.M_final,traj.h_max/1e3);
fprintf('Post-peak upward travel: seed %.2f km; SQP %.2f km (diagnostic, not a constraint).\n', ...
    upward_travel(seed_traj),upward_travel(traj));
fprintf('Saved: %s\n',outfile);

    function tr=evaluate(x)
        x=x(:)';
        if isempty(last_x) || ~isequal(x,last_x)
            pt=p; pt.gamma0=deg2rad(x(end));
            tr=simulate_trajectory(profile(x),pt,false);
            last_x=x; last_traj=tr;
            eval_count=eval_count+1;
            if tr.feasible && tr.x_final>best_traj.x_final
                best_traj=tr; best_x=x;
            end
        else
            tr=last_traj;
        end
    end
    function f=objective(x)
        tr=evaluate(x);
        f=-tr.x_final/1e6;
    end
    function [c,ceq]=constraints(x)
        tr=evaluate(x);
        c=[tr.c_viol(1)/1000;tr.c_viol(2);tr.c_viol(3)/1000];
        ceq=[];
    end
    function stop=checkpoint(~,~,state)
        stop=false;
        if strcmp(state,'iter') && best_traj.x_final>checkpoint_range+1e3
            cp_alpha=profile(best_x);
            cp_launch_angle_deg=best_x(end);
            cp_traj=best_traj;
            save(checkpoint_file,'cp_alpha','cp_launch_angle_deg','cp_traj');
            checkpoint_range=best_traj.x_final;
        end
    end
    function alpha=profile(x)
        x=x(:);
        alpha=seed_alpha+(B*x(1:n_control))';
    end
end

function [err,eq,h,V,gamma]=entry_mismatch(tr,switch_s,p)
if tr.t(end)<switch_s
    err=Inf; eq=NaN; h=NaN; V=NaN; gamma=NaN; return;
end
V=interp1(tr.t,tr.V,switch_s,'linear');
h=interp1(tr.t,tr.h,switch_s,'linear');
gamma=interp1(tr.t,tr.gamma,switch_s,'linear');
[eq,exists]=ld_equilibrium_altitude(p,V,'speed');
if exists, err=abs(h-eq)/1e3; else, err=Inf; end
end

function km=upward_travel(tr)
[~,i]=max(tr.h);
km=sum(max(0,diff(tr.h(i:end))))/1e3;
end

function tf=replay_ok(alpha,p,tr)
pc=p;pc.ode_max_step=min(0.5,p.ode_max_step/2);
pc.ode_reltol=p.ode_reltol/2;
fine=simulate_trajectory(alpha,pc,false);
tf=fine.feasible && same_trajectory(tr,fine);
end

function tf=same_trajectory(a,b)
tf=a.feasible==b.feasible && a.landed==b.landed && ...
    abs(a.x_final-b.x_final)<=max(100,0.001*a.x_final) && ...
    abs(a.h_max-b.h_max)<=50 && abs(a.n_max-b.n_max)<=0.2;
end

function q=p_with_angle(p,angle)
q=p;q.gamma0=deg2rad(angle);
end
