function main_trajectory_fmincon()
% MAIN_TRAJECTORY_FMINCON  Local SQP search over AoA knots and launch angle.
% Standalone alternative entry point. Existing GA and simulation files are
% used without modification. Requires Optimization Toolbox for fmincon.

%% User settings
LAUNCH_ANGLE_DEG = 40;
LAUNCH_ANGLE_MIN_DEG = 20;
LAUNCH_ANGLE_MAX_DEG = 60;
CONTROL_HORIZON_S = 1800;  % last AoA knot; final command is held after this
MAX_FLIGHT_TIME_S = 2500;  % ODE integration limit; must exceed actual impact
BASE_CONTROL_HORIZON_S = 1200; % preserve the original dense 0--1200 s grid
N_BASE_TIME_KNOTS = 180;    % knots through 1200 s
N_EXTRA_TIME_KNOTS = 60;   % append knots from 1200 to 1800 s
N_EARLY_KNOTS = 30;        % first 30 knots in the first 60 seconds
EARLY_WINDOW_S = 60;
N_OPT_EARLY = 10;           % SQP correction anchors in the early window
N_OPT_MID = 24;             % correction anchors from 60 to 1200 s
N_OPT_LATE = 8;             % correction anchors from 1200 to 1800 s
MAX_CORRECTION_DEG = 5;    % allowed change from the feasible seed AoA
ANGLE_SCAN_HALF_WIDTH_DEG = 8;
ANGLE_SCAN_STEP_DEG = 0.5;
RNG_SEED = 42;
MAX_ITER = 80;
MAX_EVAL = 2500;
FD_STEP_DEG = 0.03;        % relative AoA step (about 0.03 deg near 0 deg)
LAUNCH_FD_STEP_DEG = 0.05; % launch-angle perturbation near the starting angle
RESULTS_DIR_NAME = 'results';

folder = fileparts(mfilename('fullpath'));
addpath(folder);
RESULTS_DIR = fullfile(folder,RESULTS_DIR_NAME);
fprintf('Running: %s\n',mfilename('fullpath'));
fprintf('Dependencies: %s | %s\n', ...
    which('build_vehicle_params'),which('setup_waypoints'));
if exist('fmincon','file') ~= 2
    error('main_trajectory_fmincon:toolbox', ...
        'Optimization Toolbox with fmincon is required.');
end
if ~isscalar(LAUNCH_ANGLE_DEG) || ~isfinite(LAUNCH_ANGLE_DEG) || ...
        ~isscalar(LAUNCH_ANGLE_MIN_DEG) || ~isfinite(LAUNCH_ANGLE_MIN_DEG) || ...
        ~isscalar(LAUNCH_ANGLE_MAX_DEG) || ~isfinite(LAUNCH_ANGLE_MAX_DEG) || ...
        LAUNCH_ANGLE_MIN_DEG <= 0 || LAUNCH_ANGLE_MAX_DEG >= 90 || ...
        LAUNCH_ANGLE_MIN_DEG >= LAUNCH_ANGLE_DEG || ...
        LAUNCH_ANGLE_DEG >= LAUNCH_ANGLE_MAX_DEG
    error('main_trajectory_fmincon:launchBounds', ...
        'Set 0 < launch minimum < starting angle < launch maximum < 90 degrees.');
end
if ~exist(RESULTS_DIR,'dir'), mkdir(RESULTS_DIR); end
rng(RNG_SEED,'twister');

p = build_vehicle_params();
p.gamma0 = deg2rad(LAUNCH_ANGLE_DEG);
p.n_time_knots = N_BASE_TIME_KNOTS+N_EXTRA_TIME_KNOTS;
p.n_early_knots = N_EARLY_KNOTS;
p.early_window_s = EARLY_WINDOW_S;
p.seed_results_file = fullfile(RESULTS_DIR, ...
    sprintf('results_ang%d.mat',LAUNCH_ANGLE_DEG));
p.time_horizon_override_s = CONTROL_HORIZON_S;
p.T_max = max(p.T_max,MAX_FLIGHT_TIME_S);
aero_file = fullfile(folder,'aero_tables.mat');
generate_aero_tables(aero_file, ...
    fullfile(folder,'owen_tricon_v1_data.csv'),p.S_ref,p.L_ref);
p.aero = load(aero_file);
p.M_min_table = min(p.aero.Mach_vec);
p.M_max_table = max(p.aero.Mach_vec);
p = setup_waypoints(p);
if BASE_CONTROL_HORIZON_S <= EARLY_WINDOW_S || ...
        BASE_CONTROL_HORIZON_S >= CONTROL_HORIZON_S || ...
        N_BASE_TIME_KNOTS <= N_EARLY_KNOTS || ...
        N_EXTRA_TIME_KNOTS < 1
    error('main_trajectory_fmincon:gridSettings','Invalid extended grid settings.');
end
% Keep all the old knot times through BASE_CONTROL_HORIZON_S. Appending
% knots avoids changing an existing high-performing 1200 s schedule.
p.wp_times = [linspace(0,EARLY_WINDOW_S,N_EARLY_KNOTS), ...
    EARLY_WINDOW_S+(1:(N_BASE_TIME_KNOTS-N_EARLY_KNOTS))* ...
    (BASE_CONTROL_HORIZON_S-EARLY_WINDOW_S)/ ...
    (N_BASE_TIME_KNOTS-N_EARLY_KNOTS), ...
    BASE_CONTROL_HORIZON_S+(1:N_EXTRA_TIME_KNOTS)* ...
    (CONTROL_HORIZON_S-BASE_CONTROL_HORIZON_S)/N_EXTRA_TIME_KNOTS];
if p.N_wp ~= N_BASE_TIME_KNOTS+N_EXTRA_TIME_KNOTS || ...
        numel(p.wp_times) ~= p.N_wp || ...
        abs(p.wp_times(N_BASE_TIME_KNOTS)-BASE_CONTROL_HORIZON_S)>1e-8 || ...
        abs(p.wp_times(end)-CONTROL_HORIZON_S) > 1e-8 || ...
        p.T_max < MAX_FLIGHT_TIME_S || ...
        abs(p.early_window_s-EARLY_WINDOW_S) > 1e-8
    error('main_trajectory_fmincon:settingsMismatch', ...
        ['The active setup_waypoints or build_vehicle_params did not honor ' ...
         'the settings above. Check which main_trajectory_fmincon -all.']);
end
if p.T_max <= p.wp_times(end)
    error('main_trajectory_fmincon:time','The integration limit must exceed the last control knot.');
end
anchor_t = [linspace(0,p.early_window_s,N_OPT_EARLY), ...
    p.early_window_s+(1:N_OPT_MID)* ...
    (BASE_CONTROL_HORIZON_S-p.early_window_s)/N_OPT_MID, ...
    BASE_CONTROL_HORIZON_S+(1:N_OPT_LATE)* ...
    (p.wp_times(end)-BASE_CONTROL_HORIZON_S)/N_OPT_LATE];
n_control = numel(anchor_t);
B = interp1(anchor_t,eye(n_control),p.wp_times,'linear');
if any(~isfinite(B(:))) || size(B,1) ~= p.N_wp
    error('main_trajectory_fmincon:controlBasis','Invalid control correction grid.');
end
fprintf('SQP setup: %d full AoA knots; %d correction anchors plus launch angle [%.1f, %.1f] deg; horizon %.1f s; T_max %.1f s.\n', ...
    p.N_wp,n_control,LAUNCH_ANGLE_MIN_DEG,LAUNCH_ANGLE_MAX_DEG, ...
    p.wp_times(end),p.T_max);
fprintf('Preserved grid: %d knots through %.1f s; appended %d knots through %.1f s.\n', ...
    N_BASE_TIME_KNOTS,BASE_CONTROL_HORIZON_S, ...
    N_EXTRA_TIME_KNOTS,CONTROL_HORIZON_S);

%% Find a feasible start (old saved schedule is included when compatible)
seeds = fmincon_time_seeds(p);
seed_angles = LAUNCH_ANGLE_DEG*ones(size(seeds,1),1);
% Warm start from a prior SQP run even when its angle or knot horizon differs.
old_results = {sprintf('results_fmincon_start%d_h%d_launchfree.mat', ...
    LAUNCH_ANGLE_DEG,CONTROL_HORIZON_S), ...
    sprintf('results_fmincon_start%d_launchfree.mat',LAUNCH_ANGLE_DEG), ...
    sprintf('results_fmincon_ang%d.mat',LAUNCH_ANGLE_DEG)};
for k=1:numel(old_results)
    path = fullfile(RESULTS_DIR,old_results{k});
    if ~exist(path,'file'), continue; end
    saved = load(path,'p','wp_opt');
    if ~isfield(saved,'p') || ~isfield(saved,'wp_opt') || ...
            ~isfield(saved.p,'wp_times') || ...
            ~isfield(saved.p,'gamma0') || ~isfield(saved.p,'S_ref') || ...
            ~isfield(saved.p,'m') || ~isfield(saved.p,'aero') || ...
            abs(saved.p.S_ref-p.S_ref)>1e-8 || ...
            abs(saved.p.m-p.m)>1e-8 || ...
            ~isequal(saved.p.aero.CL_table,p.aero.CL_table) || ...
            ~isequal(saved.p.aero.CD_table,p.aero.CD_table) || ...
            numel(saved.p.wp_times)~=numel(saved.wp_opt) || ...
            numel(saved.p.wp_times)<2 || ...
            any(diff(saved.p.wp_times)<=0) || ...
            abs(saved.p.wp_times(1))>1e-8
        fprintf('Skipped incompatible saved trajectory: %s\n',path);
        continue;
    end
    angle = rad2deg(saved.p.gamma0);
    if angle < LAUNCH_ANGLE_MIN_DEG || angle > LAUNCH_ANGLE_MAX_DEG
        fprintf('Skipped saved launch angle outside bounds: %s\n',path);
        continue;
    end
    query = min(saved.p.wp_times(end),p.wp_times);
    alpha = interp1(saved.p.wp_times,saved.wp_opt(:)',query,'linear');
    if any(~isfinite(alpha)) || any(alpha<p.wp_lb-1e-8) || ...
            any(alpha>p.wp_ub+1e-8), continue; end
    seeds(end+1,:) = alpha; %#ok<AGROW>
    seed_angles(end+1,1) = angle; %#ok<AGROW>
    fprintf('Loaded saved SQP schedule at %.3f deg: %s\n',angle,path);
end
seed_range = -inf(size(seeds,1),1);
seed_trajs = cell(size(seeds,1),1);
for j=1:size(seeds,1)
    pj = p;
    pj.gamma0 = deg2rad(seed_angles(j));
    seed_trajs{j} = simulate_trajectory(seeds(j,:),pj,false);
    if seed_trajs{j}.feasible, seed_range(j) = seed_trajs{j}.x_final; end
end
fprintf('Physics and saved-history seeds: %d/%d feasible at the 2 s step.\n', ...
    sum(isfinite(seed_range)),numel(seed_range));
[~,order] = sort(seed_range,'descend');
seed_alpha = [];
seed_angle = NaN;
seed_traj = [];
for j=order(:)'
    if ~isfinite(seed_range(j)), break; end
    pc = p;
    pc.gamma0 = deg2rad(seed_angles(j));
    pc.ode_max_step = min(0.5,p.ode_max_step/2);
    pc.ode_reltol = p.ode_reltol/2;
    fine = simulate_trajectory(seeds(j,:),pc,false);
    if fine.feasible && same_trajectory(seed_trajs{j},fine)
        seed_alpha = seeds(j,:);
        seed_angle = seed_angles(j);
        seed_traj = seed_trajs{j};
        break;
    end
end
if isempty(seed_alpha)
    error('main_trajectory_fmincon:noFeasibleSeed', ...
        ['No physics or saved-history seed passed the 2 s and 0.5 s ' ...
         'feasibility/replay checks. Run main_trajectory_optimizer to ' ...
         'save a feasible result, or tune the time horizon/seed controls.']);
end
% A fixed-control angle scan shows whether angle alone can improve the
% verified starting schedule. It also supplies a better SQP starting angle.
scan_angles = unique([seed_angle,LAUNCH_ANGLE_DEG, ...
    max(LAUNCH_ANGLE_MIN_DEG,seed_angle-ANGLE_SCAN_HALF_WIDTH_DEG): ...
    ANGLE_SCAN_STEP_DEG: ...
    min(LAUNCH_ANGLE_MAX_DEG,seed_angle+ANGLE_SCAN_HALF_WIDTH_DEG)]);
scan_ranges = -inf(size(scan_angles));
for j=1:numel(scan_angles)
    ps = p;
    ps.gamma0 = deg2rad(scan_angles(j));
    ts = simulate_trajectory(seed_alpha,ps,false);
    if ts.feasible, scan_ranges(j) = ts.x_final; end
end
[~,scan_order] = sort(scan_ranges,'descend');
for j=scan_order(:)'
    if ~isfinite(scan_ranges(j)) || scan_ranges(j)<=seed_traj.x_final
        break;
    end
    ps = p;
    ps.gamma0 = deg2rad(scan_angles(j));
    ts = simulate_trajectory(seed_alpha,ps,false);
    ps.ode_max_step = min(0.5,p.ode_max_step/2);
    ps.ode_reltol = p.ode_reltol/2;
    tf = simulate_trajectory(seed_alpha,ps,false);
    if tf.feasible && same_trajectory(ts,tf)
        seed_angle = scan_angles(j);
        seed_traj = ts;
        break;
    end
end
fprintf('Fixed-control angle scan: %d/%d feasible; verified start %.3f deg, %.2f km.\n', ...
    sum(isfinite(scan_ranges)),numel(scan_angles),seed_angle,seed_traj.x_final/1e3);
fprintf('Verified SQP starting range: %.2f km, flight time %.1f s, peak %.2f km, load %.2f g.\n', ...
    seed_traj.x_final/1e3,seed_traj.t_flight, ...
    seed_traj.h_max/1e3,seed_traj.n_max);

%% Optimize smooth AoA corrections and launch angle. Zero correction
% reproduces the verified seed exactly, regardless of the full knot count.
seed_x = [zeros(1,n_control),seed_angle];
best_x = seed_x;
best_traj = seed_traj;
last_x = [];
last_traj = [];
eval_count = 0;
checkpoint_range = seed_traj.x_final;
checkpoint_file = fullfile(RESULTS_DIR,sprintf( ...
    'fmincon_checkpoint_start%d_h%d_launchfree.mat', ...
    LAUNCH_ANGLE_DEG,CONTROL_HORIZON_S));
% fmincon's finite-difference step is relative to max(abs(x),TypicalX).
% Dividing by the starting angle gives about LAUNCH_FD_STEP_DEG degrees.
fd_steps = [FD_STEP_DEG*ones(1,n_control), ...
    LAUNCH_FD_STEP_DEG/seed_angle];
options = optimoptions('fmincon', ...
    'Algorithm','sqp','Display','iter', ...
    'MaxIterations',MAX_ITER,'MaxFunctionEvaluations',MAX_EVAL, ...
    'FiniteDifferenceType','forward', ...
    'FiniteDifferenceStepSize',fd_steps, ...
    'TypicalX',ones(1,n_control+1), ...
    'ConstraintTolerance',1e-3, ...
    'StepTolerance',1e-4, ...
    'OptimalityTolerance',1e-3, ...
    'OutputFcn',@checkpoint);
% Serial finite differences preserve the shared trajectory cache and best
% candidate tracker. On a laptop this also avoids a process-pool startup.
tic_opt = tic;
% Enforce every actual AoA knot, even though SQP has fewer variables.
A = [B,zeros(p.N_wp,1); -B,zeros(p.N_wp,1)];
b = [p.wp_ub(:)-seed_alpha(:); seed_alpha(:)-p.wp_lb(:)];
[x_sqp,~,exitflag,output] = fmincon(@objective,seed_x, ...
    A,b,[],[],[-MAX_CORRECTION_DEG*ones(1,n_control),LAUNCH_ANGLE_MIN_DEG], ...
    [MAX_CORRECTION_DEG*ones(1,n_control),LAUNCH_ANGLE_MAX_DEG], ...
    @constraints,options);
seconds = toc(tic_opt);
% Evaluate the returned point as well; best_x retains a feasible seed
% even when SQP stops at an infeasible trial point.
returned_traj = evaluate(x_sqp);

%% Replay the best candidate at both ODE steps
chosen_x = best_x;
wp_opt = profile_from_x(chosen_x);
p.gamma0 = deg2rad(chosen_x(end));
traj = simulate_trajectory(wp_opt,p,true);
if ~traj.feasible || ~same_trajectory(best_traj,traj)
    error('main_trajectory_fmincon:replay', ...
        'The selected 2 s trajectory did not reproduce the optimizer evaluation.');
end
pc = p;
pc.ode_max_step = min(0.5,p.ode_max_step/2);
pc.ode_reltol = p.ode_reltol/2;
check = simulate_trajectory(wp_opt,pc,false);
traj.numerically_converged = check.feasible && same_trajectory(traj,check);
if ~traj.numerically_converged
    % The initial seed was already checked at both resolutions.
    chosen_x = seed_x;
    wp_opt = seed_alpha;
    p.gamma0 = deg2rad(seed_angle);
    pc.gamma0 = p.gamma0;
    traj = simulate_trajectory(wp_opt,p,true);
    check = simulate_trajectory(wp_opt,pc,false);
    traj.numerically_converged = check.feasible && same_trajectory(traj,check);
    warning('main_trajectory_fmincon:bestFailedReplay', ...
        'SQP improvement failed the finer replay; retained the verified seed.');
end
if ~traj.feasible || ~traj.numerically_converged
    error('main_trajectory_fmincon:unverified', ...
        'Neither the SQP result nor the starting seed passed final replay.');
end

% The reference curve in plot_results must share the chosen launch angle.
p_ball = p;
p_ball.abort_altitude_m = Inf;
p_ball.control_mode = 'time_knots';
p.ball.traj = simulate_trajectory( ...
    p.aero.alpha_zero_lift_deg*ones(1,p.N_wp),p_ball,false);
p.ball.x_ball = p.ball.traj.x_final;
p.ball.t_flight = p.ball.traj.t_flight;

info.method = 'fmincon sqp';
info.exitflag = exitflag;
info.output = output;
info.start_alpha = seed_alpha;
info.start_launch_angle_deg = seed_angle;
info.configured_launch_angle_deg = LAUNCH_ANGLE_DEG;
info.start_range_m = seed_traj.x_final;
info.sqp_returned_alpha = profile_from_x(x_sqp);
info.sqp_returned_launch_angle_deg = x_sqp(end);
info.best_launch_angle_deg = chosen_x(end);
info.correction_anchor_times_s = anchor_t;
info.best_correction_deg = chosen_x(1:n_control);
info.angle_scan_deg = scan_angles;
info.angle_scan_feasible_range_m = scan_ranges;
info.evaluations = eval_count;
info.seconds = seconds;
info.fine_check = check;
outfile = fullfile(RESULTS_DIR,sprintf( ...
    'results_fmincon_start%d_h%d_launchfree.mat', ...
    LAUNCH_ANGLE_DEG,CONTROL_HORIZON_S));
save(outfile,'traj','wp_opt','p','info');
T_exp = table(traj.t,traj.x/1e3,traj.h/1e3,traj.V, ...
    traj.M,traj.n,traj.alpha,rad2deg(traj.gamma), ...
    'VariableNames',{'t_s','x_km','h_km','V_ms','Mach','n_g', ...
                     'alpha_deg','gamma_deg'});
writetable(T_exp,fullfile(RESULTS_DIR, ...
    sprintf('trajectory_fmincon_start%d_h%d_launchfree.csv', ...
    LAUNCH_ANGLE_DEG,CONTROL_HORIZON_S)));
summary(1).launch_angle_deg = chosen_x(end);
summary(1).range_km = traj.x_final/1e3;
summary(1).feasible = traj.feasible;
summary(1).traj = traj;
summary(1).info = info;
plot_results(traj,summary,p);
fprintf('\nSQP: start %.2f km at %.3f deg; verified best %.2f km at %.3f deg; flight %.1f s; %.0f evaluations in %.1f min.\n', ...
    seed_traj.x_final/1e3,seed_angle,traj.x_final/1e3, ...
    chosen_x(end),traj.t_flight,eval_count,seconds/60);
fprintf('SQP exitflag %d; returned point feasible %d; retained best verified trajectory.\n', ...
    exitflag,returned_traj.feasible);
fprintf('Peak %.2f km / %.2f g; terminal Mach %.3f; feasible %d; finer replay %d.\n', ...
    traj.h_max/1e3,traj.n_max,traj.M_final,traj.feasible,traj.numerically_converged);
fprintf('Saved: %s\n',outfile);

    function tr = evaluate(x)
        x = x(:)';
        if isempty(last_x) || ~isequal(x,last_x)
            p_trial = p;
            p_trial.gamma0 = deg2rad(x(end));
            tr = simulate_trajectory(profile_from_x(x),p_trial,false);
            last_x = x;
            last_traj = tr;
            eval_count = eval_count + 1;
            if tr.feasible && tr.x_final > best_traj.x_final
                best_x = x;
                best_traj = tr;
            end
        else
            tr = last_traj;
        end
    end

    function f = objective(x)
        tr = evaluate(x);
        f = -tr.x_final/1e6; % range in 1000 km units
    end

    function [c,ceq] = constraints(x)
        tr = evaluate(x);
        % The simulator assigns c_viol(3)>=1e6 to non-impact or out-of-
        % table trajectories. An SQP trial must really reach the ground.
        c = [tr.c_viol(1)/1000; tr.c_viol(2); tr.c_viol(3)/1000];
        ceq = [];
    end

    function stop = checkpoint(~,~,state)
        stop = false;
        if strcmp(state,'iter') && best_traj.x_final > checkpoint_range + 1e3
            cp_alpha = profile_from_x(best_x);
            cp_launch_angle_deg = best_x(end);
            cp_traj = best_traj;
            save(checkpoint_file,'cp_alpha','cp_launch_angle_deg','cp_traj');
            checkpoint_range = best_traj.x_final;
        end
    end

    function alpha = profile_from_x(x)
        x = x(:);
        alpha = seed_alpha + (B*x(1:n_control))';
    end
end

function tf = same_trajectory(a,b)
tf = a.feasible == b.feasible && a.landed == b.landed && ...
    abs(a.x_final-b.x_final) <= max(100,0.001*a.x_final) && ...
    abs(a.h_max-b.h_max) <= 50 && abs(a.n_max-b.n_max) <= 0.2;
end
