function [pop,labels] = fmincon_time_seeds(p)
% FMINCON_TIME_SEEDS  Candidate open-loop schedules for local range search.
% The level-glide family is calculated with a short feedback rollout of the
% point-mass model, then sampled on the fixed AoA grid. It is NOT a flight
% controller: main_trajectory_fmincon replays every candidate open loop.

N = p.N_wp;
a0 = p.aero.alpha_zero_lift_deg;
ald = p.reference.alpha_max_ld_deg;
T = p.wp_times(end);
tk = p.wp_times;
pop = [a0*ones(1,N);ald*ones(1,N)];
labels = {'constant_zero_lift';'constant_best_ld'};
for frac = [0.55,0.75,0.9]
    u = max(0,min(1,(tk-frac*T)/(0.10*T)));
    pop(end+1,:) = ald*(1-u)+a0*u; %#ok<AGROW>
    labels{end+1,1} = sprintf('best_ld_switch_%.2f',frac); %#ok<AGROW>
end

% Gamma rises from the launch value to approximately level flight. A late
% descent is commanded once; no target path includes renewed climbs.
for turn_s = [24,26,28,32,36,44,52]
    for descent_s = [950,1100,1200,1300,1450]
        for sink_deg = [-2,-3,-4,-6,-9]
            pop(end+1,:) = level_glide_seed(p,turn_s,descent_s,sink_deg); %#ok<AGROW>
            labels{end+1,1} = sprintf('level_turn%d_down%d_sink%d', ...
                turn_s,descent_s,-sink_deg); %#ok<AGROW>
        end
    end
end
% Later turns reach a higher level-flight altitude. They also need earlier
% descents to retain the Mach-3 impact-speed margin.
for turn_s = [30,31,32]
    for descent_s = [500,600,700,800,900]
        for sink_deg = [-3,-4,-5,-8,-10]
            pop(end+1,:) = level_glide_seed(p,turn_s,descent_s,sink_deg); %#ok<AGROW>
            labels{end+1,1} = sprintf('high_turn%d_down%d_sink%d', ...
                turn_s,descent_s,-sink_deg); %#ok<AGROW>
        end
    end
end

% Explore a short zero-LIFT hold before the guided turn. This vehicle's
% zero-lift AoA is about -2.34 degrees, not zero degrees.
% Roll the held control through the model before calculating later commands;
% simply overwriting the first knots of an existing schedule would leave
% its later guidance based on the wrong state.
for hold_s = [1,2,4]
    for turn_s = [27,28,30]
        for descent_s = [600,800,1000,1100]
            for sink_deg = [-2,-4]
                pop(end+1,:) = level_glide_seed(p,turn_s,descent_s, ...
                    sink_deg,hold_s); %#ok<AGROW>
                labels{end+1,1} = sprintf('zero_lift_%ds_turn%d_down%d_sink%d', ...
                    hold_s,turn_s,descent_s,-sink_deg); %#ok<AGROW>
            end
        end
    end
end

% Existing GA result remains a candidate, so the new seed family cannot
% erase a previously discovered feasible schedule.
if isfield(p,'seed_results_file') && exist(p.seed_results_file,'file')
    saved = load(p.seed_results_file,'traj','p','wp_opt');
    if isfield(saved,'traj') && isfield(saved,'p') && ...
            isfield(saved,'wp_opt') && isfield(saved.traj,'t') && ...
            isfield(saved.p,'S_ref') && isfield(saved.p,'gamma0') && ...
            abs(saved.p.S_ref-p.S_ref)<1e-8 && ...
            abs(saved.p.gamma0-p.gamma0)<1e-8
        if isfield(saved.p,'wp_times') && ...
                numel(saved.p.wp_times)==numel(saved.wp_opt)
            old_t = saved.p.wp_times;
            old_a = saved.wp_opt;
        elseif isfield(saved.p,'wp_ranges') && ...
                numel(saved.p.wp_ranges)==numel(saved.wp_opt) && ...
                isfield(saved.traj,'x')
            old_t = saved.traj.t(:)';
            old_a = interp1(saved.p.wp_ranges,saved.wp_opt, ...
                max(saved.p.wp_ranges(1), ...
                min(saved.p.wp_ranges(end),saved.traj.x(:)')),'linear');
        else
            old_t = []; old_a = [];
        end
        if numel(old_t)>=2 && all(diff(old_t)>0)
            query = max(old_t(1),min(old_t(end),tk));
            pop(end+1,:) = interp1(old_t,old_a,query,'linear'); %#ok<AGROW>
            labels{end+1,1} = 'saved_ga'; %#ok<AGROW>
            fprintf('Seeded from saved control history: %s\n',p.seed_results_file);
        end
    end
end
pop = max(p.wp_lb,min(p.wp_ub,pop));
end

function ak = level_glide_seed(p,turn_s,descent_s,sink_deg,hold_s)
% Midpoint propagation uses the same V/gamma/h equations and CSV lift/drag
% as trajectory_eom. The true open-loop ODE replay still decides feasibility.
if nargin<5, hold_s=0; end
T = p.wp_times(end);
[~,~,~,a0] = atmosphere_1976(p.h0);
state = [p.M_launch*a0;p.gamma0;p.h0]; % speed, flight-path angle, altitude
ts = 0:1:T;
if ts(end)<T, ts(end+1)=T; end
aa = p.aero.alpha_zero_lift_deg*ones(size(ts));
for k=1:numel(ts)
    [aa(k),k1] = guided_rhs(ts(k),state,p,turn_s,descent_s,sink_deg,hold_s);
    if k==numel(ts), break; end
    step = ts(k+1)-ts(k);
    [~,k2] = guided_rhs(ts(k)+step/2,state+step*k1/2, ...
        p,turn_s,descent_s,sink_deg,hold_s);
    state = state+step*k2;
    if state(3)<0 || state(3)>86000 || abs(state(2))>1.5
        aa(k+1:end)=aa(k); break;
    end
end
ak = interp1(ts,aa,p.wp_times,'linear');
end

function [alpha,ds] = guided_rhs(t,s,p,turn_s,descent_s,sink_deg,hold_s)
V = max(1,s(1)); gamma = s(2); h = max(0,s(3));
r = p.earth_radius_m+h;
grav = p.mu_earth/r^2;
[rho,~,~,sound] = atmosphere_1976(h);
M = max(p.M_min_table,min(p.M_max_table,V/max(sound,1)));
q = 0.5*rho*V^2;
if t < turn_s
    u = t/turn_s;
    smooth = 3*u^2-2*u^3;
    target = p.gamma0*(1-smooth);
    target_rate = -p.gamma0*(6*u-6*u^2)/turn_s;
elseif t < descent_s
    target = 0; target_rate = 0;
elseif t < descent_s+100
    u = (t-descent_s)/100;
    target = deg2rad(sink_deg)*(3*u^2-2*u^3);
    target_rate = deg2rad(sink_deg)*(6*u-6*u^2)/100;
else
    target = deg2rad(sink_deg); target_rate = 0;
end
rate_cmd = target_rate+(target-gamma)/12;
cl_req = p.m*(V*rate_cmd+(grav-V^2/r)*cos(gamma))/ ...
    max(q*p.S_ref,1);
cl_limit = p.n_max*p.m*p.g/max(q*p.S_ref,1);
cl_req = max(-cl_limit,min(cl_limit,cl_req));
alpha = max(-p.alpha_max_deg,min(p.alpha_max_deg, ...
    invert_lift(cl_req,M,p.aero)));
if t<hold_s
    alpha = p.aero.alpha_zero_lift_deg;
elseif hold_s>0 && t<hold_s+3
    u = (t-hold_s)/3;
    w = 3*u^2-2*u^3;
    alpha = p.aero.alpha_zero_lift_deg*(1-w)+alpha*w;
end
[cl,cd] = aero_lookup(M,alpha,p.aero);
ds = [-cd*q*p.S_ref/p.m-grav*sin(gamma); ...
      (cl*q*p.S_ref/p.m-(grav-V^2/r)*cos(gamma))/V; ...
      V*sin(gamma)];
end

function alpha = invert_lift(target,M,aero)
cl = interp1(aero.Mach_vec,aero.CL_table,M,'linear');
crossings = find((cl(1:end-1)-target).* ...
                 (cl(2:end)-target)<=0);
if isempty(crossings)
    [~,idx]=min(abs(cl-target)); alpha=aero.alpha_vec(idx); return
end
possible=zeros(size(crossings));
for j=1:numel(crossings)
    i=crossings(j);
    d=cl(i+1)-cl(i);
    if abs(d)<1e-12
        possible(j)=aero.alpha_vec(i);
    else
        possible(j)=aero.alpha_vec(i)+(target-cl(i))* ...
            (aero.alpha_vec(i+1)-aero.alpha_vec(i))/d;
    end
end
[~,idx]=min(abs(possible-aero.alpha_zero_lift_deg));
alpha=possible(idx);
end
