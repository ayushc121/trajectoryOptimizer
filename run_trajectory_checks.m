% RUN_TRAJECTORY_CHECKS  Small consistency checks with the FOSTRAD CSV.
% Run from the folder containing the revised .m files.
clear; clc;
p = build_vehicle_params();
csv_file = fullfile(fileparts(mfilename('fullpath')), 'owen_tricon_v1_data.csv');
generate_aero_tables('check_aero_tables.mat', csv_file, p.S_ref, p.L_ref);
p.aero = load('check_aero_tables.mat');
p.M_min_table = min(p.aero.Mach_vec);
p.M_max_table = max(p.aero.Mach_vec);
p.h0 = 100;
p.gamma0 = deg2rad(-5);
p.wp_times = [0, 5];
p.T_max = 30;

assert(numel(p.aero.alpha_vec) == 81 && ...
    abs(p.aero.alpha_zero_lift_deg + 2.33831548487) < 1e-6, ...
    'CSV was not imported or zero-lift AoA is wrong.');
assert(all(all(p.aero.CL_table == repmat(p.aero.CL_table(1,:), ...
    numel(p.aero.Mach_vec), 1))), 'Mach rows are not repeated.');
% The faster one-dimensional lookup must reproduce bilinear interpolation
% on the repeated FOSTRAD Mach rows.
for M_test = [2, 4.5, 8.5]
    for alpha_test = [-9.7, -2.3, 0.1, 7.6]
        [cl,cd,cm] = aero_lookup(M_test,alpha_test,p.aero);
        cl_ref = interp2(p.aero.alpha_vec,p.aero.Mach_vec,p.aero.CL_table, ...
                         alpha_test,M_test,'linear');
        cd_ref = interp2(p.aero.alpha_vec,p.aero.Mach_vec,p.aero.CD_table, ...
                         alpha_test,M_test,'linear');
        cm_ref = interp2(p.aero.alpha_vec,p.aero.Mach_vec,p.aero.CMy_table, ...
                         alpha_test,M_test,'linear');
        assert(max(abs([cl-cl_ref,cd-cd_ref,cm-cm_ref])) < 1e-12, ...
            'The fast aero lookup changed the FOSTRAD coefficients.');
    end
end
alpha0 = p.aero.alpha_zero_lift_deg;
fast = simulate_trajectory([alpha0, alpha0], p, false);
detail = simulate_trajectory([alpha0, alpha0], p, true);
assert(fast.landed && detail.landed, 'Ground impact was not detected.');
assert(abs(fast.n_max - detail.n_max) < 0.05, ...
    'Fast and detailed load factors disagree.');
assert(abs(fast.x_final - detail.x_final) < 5, ...
    'Fast and detailed landing ranges disagree.');
assert(isequal(fast.t, detail.t) && isequal(fast.x, detail.x), ...
    'Optimization and detailed replay must use exactly the same integration.');
assert(all(detail.M > 0), 'Detailed Mach history is invalid.');

% Starting 100 m up but stopping after 0.1 s is not an impact, even though
% the final altitude is below the old 200 m heuristic.
p.T_max = 0.1;
short = simulate_trajectory([alpha0, alpha0], p, false);
assert(~short.landed && ~short.feasible, ...
    'Time-limited airborne trajectory was incorrectly declared feasible.');

% The FOSTRAD 10-degree lift exceeds the nominal 20 g constraint at launch.
p.T_max = 30;
overload = simulate_trajectory([10, 10], p, true);
assert(overload.n_max > p.n_max && ~overload.feasible, ...
    'Overload was hidden by a lift clamp.');

% Vacuum ballistic flight over a spherical Earth conserves specific energy
% and angular momentum. This detects a missing curvature term or a wrong
% surface-arc range rate without assuming any particular impact distance.
pv = p;
pv.h0 = 0; pv.gamma0 = deg2rad(40); pv.T_max = 1200;
pv.abort_altitude_m = Inf;
pv.aero.CL_table(:) = 0;
pv.aero.CD_table(:) = realmin;
vac = simulate_trajectory([0,0],pv,false);
assert(vac.landed,'Vacuum ballistic check did not impact.');
r = pv.earth_radius_m + vac.h;
energy = 0.5*vac.V.^2 - pv.mu_earth./r;
momentum = r.*vac.V.*cos(vac.gamma);
assert(max(abs(energy-energy(1)))/abs(energy(1)) < 1e-5, ...
    'Curved-Earth ballistic energy was not conserved.');
assert(max(abs(momentum-momentum(1)))/abs(momentum(1)) < 1e-5, ...
    'Curved-Earth ballistic angular momentum was not conserved.');
assert(max(abs(vac.central_angle_rad-vac.x/pv.earth_radius_m)) < 1e-12, ...
    'Surface arc and central angle disagree.');

% A nonconstant time schedule with an overload at launch must be rejected,
% while fast and detailed simulation still replay the identical states.
p.h0 = 0; p.gamma0 = deg2rad(40); p.T_max = 1200;
p.wp_times = [0 2 5 10 15 50 100 200 400 600 800];
bad_a = [10 -8.2 -5.7 -8.6 -1.3 -1.3 -2.0 -2.5 -0.3 6.0 -6.3];
fast = simulate_trajectory(bad_a,p,false);
detail = simulate_trajectory(bad_a,p,true);
assert(isequal(fast.x,detail.x) && isequal(fast.h,detail.h) && ...
    isequal(fast.c_viol,detail.c_viol), ...
    'The formerly divergent replay differs from the optimizer.');
assert(~fast.feasible, 'An infeasible ascent was accepted.');

% The new control must interpolate in time, hold its endpoint values, and
% switch the reference at actual Mach 3.
p.wp_times = [0 5 10]; p.wp_alphas = [-4 2 -1];
assert(abs(alpha_command(2.5,8,p)+1)<1e-12 && ...
       alpha_command(12,8,p)==-1,'Time interpolation or endpoint hold failed.');
pr = p; pr.control_mode = 'reference_mach_switch';
pr.reference_high_alpha_deg = 4;
pr.reference_low_alpha_deg = p.aero.alpha_zero_lift_deg;
pr.reference_switch_mach = 3;
assert(alpha_command(10,3.01,pr)==4 && ...
       abs(alpha_command(10,2.99,pr)-p.aero.alpha_zero_lift_deg)<1e-12, ...
       'Mach-triggered reference switch failed.');

p = build_vehicle_params(); p.gamma0 = deg2rad(40);
p.aero = load('check_aero_tables.mat');
p.M_min_table = min(p.aero.Mach_vec);
p.M_max_table = max(p.aero.Mach_vec);
p = setup_waypoints(p);
assert(numel(p.wp_times)==p.n_time_knots && p.wp_times(1)==0 && ...
       abs(p.wp_times(p.n_early_knots)-p.early_window_s)<1e-12 && ...
       all(diff(p.wp_times)>0) && p.T_max>p.wp_times(end), ...
       'Time grid setup failed.');
assert(max(abs(diff(p.wp_times(1:p.n_early_knots))- ...
    p.early_window_s/(p.n_early_knots-1)))<1e-12, ...
    'Early time knot density is wrong.');

fprintf('Trajectory consistency checks passed.\n');
delete('check_aero_tables.mat');
