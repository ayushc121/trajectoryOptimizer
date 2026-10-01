function [h_eq,exists] = ld_equilibrium_altitude(p,value,mode)
% LD_EQUILIBRIUM_ALTITUDE  Instantaneous level-flight balance at best L/D.
% value is airspeed [m/s] in 'speed' mode or Mach in 'mach' mode. The
% equilibrium is instantaneous; drag changes airspeed, so it is not a
% constant-altitude solution for a long unpowered glide.
if nargin<3, mode='speed'; end
if ~strcmp(mode,'speed') && ~strcmp(mode,'mach')
    error('ld_equilibrium_altitude:mode','Use speed or mach.');
end
alpha = p.reference.alpha_max_ld_deg;
f = @(h) lift_excess(h,p,value,mode,alpha);
lo = f(0);
hi = f(p.h_max);
exists = isfinite(lo) && isfinite(hi) && lo*hi<=0;
if exists
    h_eq = fzero(f,[0,p.h_max]);
else
    h_eq = NaN;
end
end

function excess = lift_excess(h,p,value,mode,alpha)
[rho,~,~,sound] = atmosphere_1976(h);
if strcmp(mode,'mach')
    V = value*sound;
else
    V = value;
end
M = V/sound;
if M<p.M_min_table || M>p.M_max_table
    excess = NaN; return;
end
CL = aero_lookup(M,alpha,p.aero);
r = p.earth_radius_m+h;
g = p.mu_earth/r^2;
excess = 0.5*rho*V^2*p.S_ref*CL - ...
    p.m*(g-V^2/r);
end
