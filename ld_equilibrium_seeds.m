function [pop,labels,diagnostics] = ld_equilibrium_seeds(p)
% LD_EQUILIBRIUM_SEEDS  Turn toward cruise, hold best-L/D AoA, then descend.
% All candidates are open-loop schedules; the full ODE decides feasibility.
turns = [18,20,22,24,26,28,30,32];
switch_delays = [0,10,20];
descent_times = [1100,1250,1400,1550,1700];
T = p.wp_times;
ald = p.reference.alpha_max_ld_deg;
a0 = p.aero.alpha_zero_lift_deg;
n = numel(turns)*numel(switch_delays)*numel(descent_times);
pop = zeros(n,p.N_wp);
labels = cell(n,1);
diagnostics = repmat(struct('switch_s',0,'h_switch_m',0, ...
    'V_switch_ms',0,'h_eq_m',NaN,'mismatch_km',NaN),n,1);
k = 0;
for turn = turns
    for delay = switch_delays
        switch_s = turn+delay;
        [early_t,early_alpha,state] = pitch_over(p,turn,switch_s);
        [h_eq,exists] = ld_equilibrium_altitude(p,state(1),'speed');
        for descend = descent_times
            k = k+1;
            alpha = interp1(early_t,early_alpha,min(T,switch_s),'linear');
            alpha(T>=switch_s) = ald;
            u = max(0,min(1,(T-descend)/20));
            u = 3*u.^2-2*u.^3;
            alpha = (1-u).*alpha+u*a0;
            pop(k,:) = max(p.wp_lb,min(p.wp_ub,alpha));
            labels{k} = sprintf('ld_turn%d_switch%d_down%d', ...
                turn,switch_s,descend);
            diagnostics(k).switch_s = switch_s;
            diagnostics(k).h_switch_m = state(3);
            diagnostics(k).V_switch_ms = state(1);
            if exists
                diagnostics(k).h_eq_m = h_eq;
                diagnostics(k).mismatch_km = (state(3)-h_eq)/1e3;
            end
        end
    end
end
end

function [ts,aa,state] = pitch_over(p,turn_s,switch_s)
% Generate a smooth initial turn with the same V/gamma/h equations as the
% full simulator. The first zero-crossing of gamma can be tuned by turn_s.
[~,~,~,sound] = atmosphere_1976(p.h0);
state = [p.M_launch*sound;p.gamma0;p.h0];
ts = 0:1:switch_s;
aa = zeros(size(ts));
for j=1:numel(ts)
    [aa(j),k1] = rhs(ts(j),state,p,turn_s);
    if j==numel(ts), break; end
    [~,k2] = rhs(ts(j)+0.5,state+0.5*k1,p,turn_s);
    state = state+k2;
end
end

function [alpha,ds] = rhs(t,s,p,turn_s)
V=max(1,s(1)); gamma=s(2); h=max(0,s(3));
r=p.earth_radius_m+h; grav=p.mu_earth/r^2;
[rho,~,~,sound]=atmosphere_1976(h);
M=max(p.M_min_table,min(p.M_max_table,V/sound));
q=0.5*rho*V^2;
if t<turn_s
    u=t/turn_s;
    blend=3*u^2-2*u^3;
    target=p.gamma0*(1-blend);
    rate=-p.gamma0*(6*u-6*u^2)/turn_s;
else
    target=0; rate=0;
end
rate_cmd=rate+(target-gamma)/12;
cl_req=p.m*(V*rate_cmd+(grav-V^2/r)*cos(gamma))/max(q*p.S_ref,1);
cl_max=p.n_max*p.m*p.g/max(q*p.S_ref,1);
cl_req=max(-cl_max,min(cl_max,cl_req));
alpha=max(-p.alpha_max_deg,min(p.alpha_max_deg, ...
    lift_inverse(cl_req,M,p.aero)));
[cl,cd]=aero_lookup(M,alpha,p.aero);
ds=[-cd*q*p.S_ref/p.m-grav*sin(gamma); ...
    (cl*q*p.S_ref/p.m-(grav-V^2/r)*cos(gamma))/V; ...
    V*sin(gamma)];
end

function alpha=lift_inverse(target,M,aero)
cl=interp1(aero.Mach_vec,aero.CL_table,M,'linear');
cross=find((cl(1:end-1)-target).*(cl(2:end)-target)<=0);
if isempty(cross)
    [~,idx]=min(abs(cl-target)); alpha=aero.alpha_vec(idx); return;
end
possible=zeros(size(cross));
for k=1:numel(cross)
    i=cross(k);
    delta=cl(i+1)-cl(i);
    if abs(delta)<1e-12
        possible(k)=aero.alpha_vec(i);
    else
        possible(k)=aero.alpha_vec(i)+(target-cl(i))* ...
            (aero.alpha_vec(i+1)-aero.alpha_vec(i))/delta;
    end
end
[~,idx]=min(abs(possible-aero.alpha_zero_lift_deg));
alpha=possible(idx);
end
