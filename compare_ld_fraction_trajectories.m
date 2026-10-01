function compare_ld_fraction_trajectories()
% COMPARE_LD_FRACTION_TRAJECTORIES  Illustrate reduced lift-to-drag ratio.
% Runs the existing curved-Earth ODE for ten hypothetical coefficient sets.
% No GA or fmincon is used: screen a fixed family of analytical glide seeds.

%% Settings (change here without editing the shared simulation files)
LAUNCH_ANGLE_DEG = 40;
CONTROL_HORIZON_S = 1700;
MAX_FLIGHT_TIME_S = 1800;
N_TIME_KNOTS = 150;
N_EARLY_KNOTS = 20;
EARLY_WINDOW_S = 60;
FRACTIONS = (1:10)/10;
RESULTS_DIR = 'results';

folder = fileparts(mfilename('fullpath'));
addpath(folder);
if ~exist(RESULTS_DIR,'dir'), mkdir(RESULTS_DIR); end
p = build_vehicle_params();
p.gamma0 = deg2rad(LAUNCH_ANGLE_DEG);
p.n_time_knots = N_TIME_KNOTS;
p.n_early_knots = N_EARLY_KNOTS;
p.early_window_s = EARLY_WINDOW_S;
p.time_horizon_override_s = CONTROL_HORIZON_S;
p.T_max = max(p.T_max,MAX_FLIGHT_TIME_S);
aero_file = fullfile(RESULTS_DIR,'ld_comparison_aero.mat');
generate_aero_tables(aero_file, ...
    fullfile(folder,'owen_tricon_v1_data.csv'),p.S_ref,p.L_ref);
p.aero = load(aero_file);
p.M_min_table = min(p.aero.Mach_vec);
p.M_max_table = max(p.aero.Mach_vec);
p = setup_waypoints(p);
base_aero = p.aero;
base_ld = base_aero.CL_table(1,:) ./ base_aero.CD_table(1,:);
max_ld = max(base_ld);

n = numel(FRACTIONS);
fractions_col = FRACTIONS(:);
feasible = false(n,1);
range_km = nan(n,1);
diagnostic_range_km = nan(n,1);
peak_km = nan(n,1);
flight_s = nan(n,1);
final_mach = nan(n,1);
peak_g = nan(n,1);
landed = false(n,1);
aero_valid = false(n,1);
chosen = cell(n,1);
controls = cell(n,1);
status = cell(n,1);
trials = zeros(n,1);

for i=1:n
    f = FRACTIONS(i);
    pk = p;
    pk.aero = base_aero;
    pk.aero.CL_table = f * base_aero.CL_table;
    pk.aero.LD_source = f * base_aero.LD_source;
    % CD, Sref, mass, initial state, CMy, and heating metadata stay fixed.
    % Zero-lift AoA and the maximizing AoA for positive L/D are unchanged.
    assert(max(abs(pk.aero.CL_table(:)-f*base_aero.CL_table(:))) < 1e-12);
    assert(isequal(pk.aero.CD_table,base_aero.CD_table));

    bank = make_guided_bank(pk);
    fprintf('Fraction %d/10: checking %d analytical schedules...\n',i,size(bank,1));
    fraction_timer = tic;
    found_alpha = zeros(0,pk.N_wp);
    found_ranges = zeros(0,1);
    nearest = [];
    nearest_alpha = [];
    nearest_violation = Inf;
    for k=1:size(bank,1)
        tr = simulate_trajectory(bank(k,:),pk,false);
        if tr.feasible
            found_alpha(end+1,:) = bank(k,:); %#ok<AGROW>
            found_ranges(end+1,1) = tr.x_final; %#ok<AGROW>
        else
            v = violation_score(tr,pk);
            if v < nearest_violation
                nearest = tr;
                nearest_alpha = bank(k,:);
                nearest_violation = v;
            end
        end
        if mod(k,50)==0
            fprintf('  checked %d/%d; feasible at 2 s: %d\n', ...
                k,size(bank,1),numel(found_ranges));
        end
    end
    trials(i) = size(bank,1);
    [~,order] = sort(found_ranges,'descend');
    for j=order(:)'
        pc = pk;
        pc.ode_max_step = min(0.5,pk.ode_max_step/2);
        pc.ode_reltol = pk.ode_reltol/2;
        coarse = simulate_trajectory(found_alpha(j,:),pk,false);
        fine = simulate_trajectory(found_alpha(j,:),pc,false);
        if fine.feasible && same_trajectory(coarse,fine)
            controls{i} = found_alpha(j,:);
            chosen{i} = simulate_trajectory(controls{i},pk,true);
            feasible(i) = true;
            status{i} = 'FEASIBLE';
            break;
        end
    end
    if ~feasible(i)
        status{i} = 'NO VERIFIED SEED';
        if ~isempty(nearest)
            controls{i} = nearest_alpha;
            chosen{i} = simulate_trajectory(nearest_alpha,pk,true);
        end
    end
    if ~isempty(chosen{i})
        diagnostic_range_km(i) = chosen{i}.x_final/1e3;
        peak_km(i) = chosen{i}.h_max/1e3;
        flight_s(i) = chosen{i}.t_flight;
        final_mach(i) = chosen{i}.M_final;
        peak_g(i) = chosen{i}.n_max;
        landed(i) = chosen{i}.landed;
        aero_valid(i) = chosen{i}.aero_valid;
        if feasible(i), range_km(i) = chosen{i}.x_final/1e3; end
    end
    fprintf('%3.0f%% CL: %s; path %.1f km; peak %.1f km; time %.1f s; M_end %.2f; load %.2f g; %.1f s elapsed\n', ...
        100*f,status{i},diagnostic_range_km(i),peak_km(i), ...
        flight_s(i),final_mach(i),peak_g(i),toc(fraction_timer));
end

summary = table(fractions_col,fractions_col*max_ld,feasible,landed,aero_valid, ...
    range_km,diagnostic_range_km,peak_km,flight_s,final_mach,peak_g,trials,status, ...
    'VariableNames',{'CL_fraction','Peak_L_over_D','Feasible', ...
    'Landed','Aero_valid','Feasible_range_km','Path_distance_km', ...
    'Peak_altitude_km','Flight_time_s','Final_Mach', ...
    'Peak_load_g','Seeds_tested','Status'});
fprintf('\nRepresentative feasible seed trajectories (not optimized ranges):\n');
disp(summary);
fprintf('Feasible: %d of %d. A dashed curve marks an infeasible diagnostic flight.\n', ...
    sum(feasible),n);
colors = lines(n);

fig1 = figure('Name','L/D fraction trajectories','Position',[80 80 1100 650]);
hold on; grid on; box on;
for i=1:n
    if isempty(chosen{i}), continue; end
    tr = chosen{i};
    if feasible(i), style = '-'; else, style = '--'; end
    label = sprintf('%.1f of baseline L/D',FRACTIONS(i));
    if ~feasible(i), label = [label ' (infeasible)']; end
    plot(tr.x/1e3,tr.h/1e3,'Color',colors(i,:), ...
        'LineStyle',style,'LineWidth',1.8,'DisplayName',label);
end
yline(p.h_max/1e3,'k--','30 km ceiling','HandleVisibility','off');
xlabel('Surface-arc downrange [km]'); ylabel('Altitude [km]');
title('Guided seed trajectories at ten fractions of baseline C_L');
legend('Location','eastoutside');
saveas(fig1,fullfile(RESULTS_DIR,'ld_fraction_trajectories.png'));

fig2 = figure('Name','Range versus L/D','Position',[120 120 850 520]);
hold on; grid on; box on;
plot(fractions_col(feasible)*max_ld,range_km(feasible),'ko-', ...
    'MarkerFaceColor',[0.2 0.5 0.9],'LineWidth',1.5);
if any(~feasible)
    plot(fractions_col(~feasible)*max_ld,diagnostic_range_km(~feasible), ...
        'rx','MarkerSize',9,'LineWidth',1.5);
    legend({'Verified feasible seed','Diagnostic only'},'Location','best');
end
for i=find(feasible(:))'
    text(FRACTIONS(i)*max_ld,range_km(i),sprintf('  %d/10',i));
end
xlabel('Maximum positive L/D of scaled coefficient curve');
ylabel('Downrange distance [km]');
title('Heuristic seed range versus scaled L/D');
saveas(fig2,fullfile(RESULTS_DIR,'ld_fraction_range.png'));

fig3 = figure('Name','Scaled L/D curves','Position',[160 160 850 520]);
hold on; grid on; box on;
for i=1:n
    plot(base_aero.alpha_vec,FRACTIONS(i)*base_ld, ...
        'Color',colors(i,:),'LineWidth',1.4, ...
        'DisplayName',sprintf('%d/10',i));
end
xlabel('Angle of attack [deg]'); ylabel('L/D = C_L/C_D');
title('Imported L/D curve after scaling C_L');
legend('Location','eastoutside');
saveas(fig3,fullfile(RESULTS_DIR,'ld_fraction_curves.png'));

writetable(summary,fullfile(RESULTS_DIR,'ld_fraction_summary.csv'));
save(fullfile(RESULTS_DIR,'ld_fraction_comparison.mat'), ...
    'summary','chosen','controls','FRACTIONS','p');
fprintf('Saved figures, summary CSV, and trajectories in %s.\n',RESULTS_DIR);
end

function bank = make_guided_bank(p)
% Fixed candidate family for every fraction. No GA or gradient search.
N = p.N_wp;
a0 = p.aero.alpha_zero_lift_deg;
ald = p.reference.alpha_max_ld_deg;
tk = p.wp_times;
T = tk(end);
bank = [a0*ones(1,N); ald*ones(1,N)];
for switch_frac = [0.55,0.75,0.9]
    u = max(0,min(1,(tk-switch_frac*T)/(0.10*T)));
    bank(end+1,:) = ald*(1-u)+a0*u; %#ok<AGROW>
end
for peak_fraction = [0.30,0.45,0.60,0.75,0.90]
    for descent_fraction = [0.50,0.70,0.85,0.98]
        for rise_factor = [0.65,1.00,1.45]
            bank(end+1,:) = guided_seed(p, ...
                peak_fraction*p.h_max,descent_fraction,rise_factor); %#ok<AGROW>
        end
    end
end
% Analytically select an initial pitch-down AoA by desired aerodynamic
% load, then blend to a cruise AoA. This catches low-lift vehicles for
% which the longer altitude-hold seed spends too much speed in high drag.
[rho0,~,~,a0] = atmosphere_1976(p.h0);
V0 = p.M_launch*a0;
q0 = 0.5*rho0*V0^2;
for turn_load_g = [16,19]
    cl_turn = -turn_load_g*p.m*p.g/(q0*p.S_ref);
    alpha_turn = invert_lift(cl_turn,p.M_launch,p.aero);
    for switch_s = [8,12,16,20,26]
        for alpha_cruise = [-2.25,-1.75,-1.25,-0.75,0]
            for ramp_s = [6,12]
                u = max(0,min(1,(tk-switch_s)/ramp_s));
                bank(end+1,:) = alpha_turn*(1-u)+alpha_cruise*u; %#ok<AGROW>
            end
        end
    end
end
bank = max(p.wp_lb,min(p.wp_ub,bank));
end

function ak = guided_seed(p,peak,descent_fraction,rise_factor)
% Provisional altitude tracking, converted to open-loop time-knot AoA.
T = p.wp_times(end);
[~,~,~,a0] = atmosphere_1976(p.h0);
V = p.M_launch*a0;
gamma = p.gamma0;
h = p.h0;
vz0 = V*sin(gamma);
tr = max(25,min(160,rise_factor*2.2*(peak-p.h0)/max(vz0,1)));
tr = min(tr,0.45*T);
td = max(tr+15,descent_fraction*T);
td = min(td,T-10);
dt = min(1,max(0.25,T/2400));
ts = 0:dt:T;
if ts(end)<T, ts(end+1)=T; end
aa = p.aero.alpha_zero_lift_deg*ones(size(ts));
for k=1:numel(ts)
    t=ts(k);
    if t<tr
        delta=peak-p.h0;
        c2=3*delta/tr^2-2*vz0/tr;
        c3=vz0/tr^2-2*delta/tr^3;
        hd=p.h0+vz0*t+c2*t^2+c3*t^3;
        vd=vz0+2*c2*t+3*c3*t^2;
        ad=2*c2+6*c3*t;
    elseif t<td
        hd=peak; vd=0; ad=0;
    else
        u=min(1,(t-td)/(T-td));
        hd=peak*(1-3*u^2+2*u^3);
        vd=peak*(-6*u+6*u^2)/(T-td);
        ad=peak*(-6+12*u)/(T-td)^2;
    end
    r=p.earth_radius_m+max(h,0);
    gravity=p.mu_earth/r^2;
    [rho,~,~,a] = atmosphere_1976(max(0,h));
    Mach=max(p.M_min_table,min(p.M_max_table,V/max(a,1)));
    q=0.5*rho*V^2;
    gamma_des=asin(max(-0.85,min(0.85,vd/max(V,1)))) + ...
        max(-0.2,min(0.2,(hd-h)/15000));
    gamma_rate=(gamma_des-gamma)/10 + ...
        ad/(max(V,1)*max(0.3,abs(cos(gamma))));
    cl_req=p.m*(V*gamma_rate+(gravity-V^2/r)*cos(gamma))/ ...
        max(q*p.S_ref,1);
    cl_limit=p.n_max*p.m*p.g/max(q*p.S_ref,1);
    cl_req=max(-cl_limit,min(cl_limit,cl_req));
    aa(k)=max(-p.alpha_max_deg,min(p.alpha_max_deg, ...
        invert_lift(cl_req,Mach,p.aero)));
    [cl,cd]=aero_lookup(Mach,aa(k),p.aero);
    if k==numel(ts), break; end
    step=ts(k+1)-ts(k);
    dV=-cd*q*p.S_ref/p.m-gravity*sin(gamma);
    dgamma=(cl*q*p.S_ref/p.m-(gravity-V^2/r)*cos(gamma))/max(V,1);
    dh=V*sin(gamma);
    V=max(100,V+step*dV);
    gamma=gamma+step*dgamma;
    h=h+step*dh;
    if h<0 || h>86000 || abs(gamma)>1.5
        aa(k+1:end)=aa(k); break;
    end
end
ak=interp1(ts,aa,p.wp_times,'linear');
end

function alpha = invert_lift(target,M,aero)
cl=interp1(aero.Mach_vec,aero.CL_table,M,'linear');
crossings=find((cl(1:end-1)-target).* ...
    (cl(2:end)-target)<=0);
if isempty(crossings)
    [~,idx]=min(abs(cl-target));
    alpha=aero.alpha_vec(idx);
    return;
end
possible=zeros(size(crossings));
for j=1:numel(crossings)
    idx=crossings(j);
    d=cl(idx+1)-cl(idx);
    if abs(d)<1e-12
        possible(j)=aero.alpha_vec(idx);
    else
        possible(j)=aero.alpha_vec(idx)+(target-cl(idx))* ...
            (aero.alpha_vec(idx+1)-aero.alpha_vec(idx))/d;
    end
end
[~,idx]=min(abs(possible-aero.alpha_zero_lift_deg));
alpha=possible(idx);
end

function score = violation_score(tr,p)
score = max(0,tr.c_viol(1))/p.h_max + ...
    max(0,tr.c_viol(2))/p.n_max + ...
    max(0,tr.c_viol(3))/p.V_min_impact;
end

function tf = same_trajectory(a,b)
tf = a.feasible==b.feasible && a.landed==b.landed && ...
    abs(a.x_final-b.x_final) <= max(100,0.001*a.x_final) && ...
    abs(a.h_max-b.h_max) <= 50 && abs(a.n_max-b.n_max) <= 0.2;
end
