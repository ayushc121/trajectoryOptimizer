% FMINCON_TIME_SEEDS  Generate the same physical seed family as the GA.
% This copy is independent so the existing optimizer files remain unchanged.
function pop = fmincon_time_seeds(p)
N = p.N_wp;
a0 = p.aero.alpha_zero_lift_deg;
ald = p.reference.alpha_max_ld_deg;
T = p.wp_times(end);
tk = p.wp_times;
bases = [a0*ones(1,N); ald*ones(1,N)];
% Positive-lift cruise followed by zero-lift descent; these complement
% closed-loop-generated loft/hold targets.
for frac = [0.55,0.75,0.9]
    a = ald*ones(1,N);
    u = max(0,min(1,(tk-frac*T)/(0.10*T)));
    a = a.*(1-u)+a0.*u;
    bases(end+1,:) = a; %#ok<AGROW>
end
for peak = [0.72,0.85,0.96]
    for descend = [0.70,0.85,0.98]
        bases(end+1,:) = guided_time_seed(p,peak*p.h_max,descend); %#ok<AGROW>
    end
end
% Optional saved control schedule: project actual control versus time onto
% the new grid, then re-evaluate it under today's model and constraints.
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
            bases(end+1,:) = interp1(old_t,old_a,query,'linear'); %#ok<AGROW>
            fprintf('Seeded from saved control history: %s\n',p.seed_results_file);
        end
    end
end
bases = max(p.wp_lb,min(p.wp_ub,bases));
pop = zeros(p.ga_pop_size,N);
n_base = min(size(bases,1),p.ga_pop_size);
pop(1:n_base,:) = bases(1:n_base,:);
for k=n_base+1:p.ga_pop_size
    idx = randi(size(bases,1));
    amp = p.ga_seed_perturb_deg*(0.5+2*rand);
    % Adjacent knots get correlated perturbations so the initial controls
    % are smoother than independent white noise, while GA may still vary
    % every control independently.
    noise = randn(1,N);
    noise = conv(noise,[0.25 0.5 0.25],'same');
    pop(k,:) = max(p.wp_lb,min(p.wp_ub,bases(idx,:)+amp*noise));
end
end

function ak = guided_time_seed(p,peak,descend_frac)
% Integrate a provisional altitude feedback law and sample its AoA commands
% on the common time grid. The actual GA evaluation is always open loop.
T = p.wp_times(end);
[~,~,~,a0] = atmosphere_1976(p.h0);
V = p.M_launch*a0;
gamma = p.gamma0;
h = p.h0;
vz0 = V*sin(gamma);
tr = max(35,min(130,2.2*(peak-p.h0)/max(vz0,1)));
tr = min(tr,0.45*T);
td = max(tr+15,descend_frac*T);
td = min(td,T-10);
dt = min(1,max(0.25,T/2400));
ts = 0:dt:T;
if ts(end)<T, ts(end+1)=T; end
aa = p.aero.alpha_zero_lift_deg*ones(size(ts));
for k=1:numel(ts)
    t = ts(k);
    if t<tr
        u=t/tr;
        % Cubic Hermite: start at actual launch vertical speed and arrive
        % at target altitude with zero vertical speed.
        dh=peak-p.h0;
        c2=3*dh/tr^2-2*vz0/tr;
        c3=vz0/tr^2-2*dh/tr^3;
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
    grav=p.mu_earth/r^2;
    [rho,~,~,a] = atmosphere_1976(max(0,h));
    M=max(p.M_min_table,min(p.M_max_table,V/max(a,1)));
    q=0.5*rho*V^2;
    gamma_des=asin(max(-0.85,min(0.85,vd/max(V,1)))) + ...
        max(-0.2,min(0.2,(hd-h)/15000));
    desired_rate=(gamma_des-gamma)/10 + ...
        ad/(max(V,1)*max(0.3,abs(cos(gamma))));
    cl_req=p.m*(V*desired_rate+(grav-V^2/r)*cos(gamma))/ ...
        max(q*p.S_ref,1);
    cl_limit=p.n_max*p.m*p.g/max(q*p.S_ref,1);
    cl_req=max(-cl_limit,min(cl_limit,cl_req));
    aa(k)=max(-p.alpha_max_deg,min(p.alpha_max_deg, ...
        invert_lift(cl_req,M,p.aero)));
    [cl,cd]=aero_lookup(M,aa(k),p.aero);
    if k==numel(ts), break; end
    step=ts(k+1)-ts(k);
    dv=-cd*q*p.S_ref/p.m-grav*sin(gamma);
    dg=(cl*q*p.S_ref/p.m-(grav-V^2/r)*cos(gamma))/max(V,1);
    dhdt=V*sin(gamma);
    V=max(100,V+step*dv);
    gamma=gamma+step*dg;
    h=h+step*dhdt;
    if h<0 || h>86000 || abs(gamma)>1.5
        aa(k+1:end)=aa(k); break
    end
end
ak=interp1(ts,aa,p.wp_times,'linear');
end

function alpha = invert_lift(target,M,aero)
cl=interp1(aero.Mach_vec,aero.CL_table,M,'linear');
crossings=find((cl(1:end-1)-target).* ...
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

