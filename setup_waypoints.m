function p = setup_waypoints(p)
% SETUP_WAYPOINTS  Set a fixed time grid before the single GA run.
% The best-L/D reference is a timing probe, not a feasible trajectory or
% a mathematical upper bound on all possible flight times.
if ~isfield(p.aero,'alpha_zero_lift_deg')
    error('setup_waypoints:zeroLift','Aerodynamic table needs a zero-lift AoA.');
end
if p.n_time_knots < 3 || p.n_early_knots < 2 || ...
        p.n_time_knots <= p.n_early_knots || p.early_window_s <= 0
    error('setup_waypoints:grid','Need at least two early knots and one later knot.');
end

valid_alpha = p.aero.alpha_vec >= -p.alpha_max_deg & ...
              p.aero.alpha_vec <= p.alpha_max_deg;
ld = p.aero.CL_table(1,:) ./ p.aero.CD_table(1,:);
ld(~valid_alpha) = -Inf;
[best_ld, i] = max(ld);
if ~isfinite(best_ld) || best_ld <= 0
    error('setup_waypoints:aero','No positive-L/D AoA lies within the allowed range.');
end
lo = max(1,i-1); hi = min(numel(ld),i+1);
alpha_fine = linspace(p.aero.alpha_vec(lo),p.aero.alpha_vec(hi),101);
cl = interp1(p.aero.alpha_vec,p.aero.CL_table(1,:),alpha_fine);
cd = interp1(p.aero.alpha_vec,p.aero.CD_table(1,:),alpha_fine);
[~,k] = max(cl./cd);
alpha_ld = alpha_fine(k);

% The reference mode changes AoA when actual Mach crosses the threshold.
probe = p;
probe.wp_times = [0,1]; % validated by simulator, ignored in reference mode
probe.control_mode = 'reference_mach_switch';
probe.reference_high_alpha_deg = alpha_ld;
probe.reference_low_alpha_deg = p.aero.alpha_zero_lift_deg;
probe.reference_switch_mach = p.M_min_impact;
probe.abort_altitude_m = 86000; % stop before using clamped high-altitude atmosphere
probe.T_max = p.reference_t_max;
ref = simulate_trajectory([alpha_ld,alpha_ld],probe,false);
fprintf('Best-L/D reference: alpha %.3f deg, L/D %.3f, t %.1f s, impact %d, peak %.1f km.\n', ...
        alpha_ld, max(cl./cd), ref.t_flight, ref.landed, ref.h_max/1e3);

% If reference misses impact or leaves data coverage, use a conservative
% fallback and tell the operator. A manual override is available.
if ~isempty(p.time_horizon_override_s)
    horizon = p.time_horizon_override_s;
elseif ref.landed && ref.aero_valid && ref.h_max <= 86000
    horizon = p.time_horizon_factor * ref.t_flight;
else
    warning('setup_waypoints:referenceIncomplete', ...
        ['Best-L/D reference did not produce a usable impact within its ' ...
         'aerodynamic/atmospheric coverage. Using its elapsed time and ' ...
         'the configured T_max; inspect the horizon before optimizing.']);
    horizon = max(p.T_max, ref.t_flight);
end
if ~isscalar(horizon) || ~isfinite(horizon) || ...
        horizon <= p.early_window_s
    error('setup_waypoints:horizon','Choose a horizon greater than the early window.');
end
p.wp_times = [linspace(0,p.early_window_s,p.n_early_knots), ...
    p.early_window_s + (1:(p.n_time_knots-p.n_early_knots))* ...
    (horizon-p.early_window_s)/(p.n_time_knots-p.n_early_knots)];
p.N_wp = numel(p.wp_times);
p.N_free_wp = p.N_wp;
p.wp_lb = -p.alpha_max_deg*ones(1,p.N_wp);
p.wp_ub =  p.alpha_max_deg*ones(1,p.N_wp);
p.T_max = max(p.T_max,1.05*horizon);

zero = p;
zero.abort_altitude_m = Inf;
zero.control_mode = 'time_knots';
p.ball.traj = simulate_trajectory(p.aero.alpha_zero_lift_deg*ones(1,p.N_wp),zero,false);
p.ball.x_ball = p.ball.traj.x_final;
p.ball.t_flight = p.ball.traj.t_flight;
p.reference = struct('alpha_max_ld_deg',alpha_ld,'traj',ref, ...
    'time_horizon_s',horizon,'impact_valid', ...
    ref.landed && ref.aero_valid && ref.h_max <= 86000);
fprintf('Time grid: %d knots; %d knots in first %.1f s; last at %.1f s; T_max %.1f s.\n', ...
    p.N_wp,p.n_early_knots,p.early_window_s,horizon,p.T_max);
end
