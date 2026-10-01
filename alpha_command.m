function alpha_deg = alpha_command(t, M_actual, p)
% ALPHA_COMMAND  Same AoA in the ODE and load reconstruction.
% Reference-only mode switches at actual Mach 3; GA controls use time.
if isfield(p,'control_mode') && strcmp(p.control_mode,'reference_mach_switch')
    alpha_deg = p.reference_high_alpha_deg + zeros(size(M_actual));
    alpha_deg(M_actual <= p.reference_switch_mach) = p.reference_low_alpha_deg;
else
    tq = max(p.wp_times(1),min(p.wp_times(end),t));
    alpha_deg = interp1(p.wp_times,p.wp_alphas,tq,'linear');
end
alpha_deg = max(-p.alpha_max_deg,min(p.alpha_max_deg,alpha_deg));
end
