function [wp_opt, range_opt, info] = hybrid_optimize(p, N_random_starts)
% HYBRID_OPTIMIZE  Two-phase trajectory optimiser.
%
% -----------------------------------------------------------------------
% PHASE 1 — Exploratory GA  (joint alpha + waypoint position)
%
%   Design variables: [alpha(1:N_wp), frac(1:N_wp)]  — 2*N_wp = 22 total.
%   Fracs are SORTED inside the objective so any permutation is equivalent.
%   Soft quadratic penalties on altitude and terminal-speed constraints.
%   Stopping: custom best-fitness stall OutputFcn.
%   The entire final population is captured for Phase 2 seeding.
%
% -----------------------------------------------------------------------
% PHASE 2 — Polish GA  (alpha only, waypoints fixed from Phase 1)
%
%   Initial population: alpha columns extracted from Phase 1 final population.
%   Uses the same objective form but with:
%     - Higher penalty weights  (ga_w_alt * polish_penalty_mult)
%     - Tighter constraint margins  (h_max_ga_p2, V_min_ga_p2)
%     - Higher selection pressure  (more elite, higher tournament)
%   This drives the population toward feasibility without giving up range,
%   and cannot produce the degenerate "crash-early" solution that afflicted
%   the previous fmincon Phase 2 (since range is still in the objective).
%
%   If the Polish GA produces a feasible result and N_random_starts > 0,
%   a final fmincon SQP pass (Phase 3) refines from that feasible point.
%   Set N_random_starts = 0 to skip fmincon entirely.
%
% -----------------------------------------------------------------------
% CRITICAL NOTE ON WAYPOINT POSITIONS:
%   p.wp_ranges is updated inside this function (local copy).
%   The caller MUST use info.wp_ranges_opt for any subsequent
%   simulate_trajectory calls:
%       p.wp_ranges = info.wp_ranges_opt;
%       traj = simulate_trajectory(wp_opt, p, true);

if nargin < 2, N_random_starts = 2; end

N_wp = p.N_free_wp;    % 11
N_ga = N_wp + N_wp;    % 22: alpha + fracs

%% ==========================================================================
%% Phase 1: Exploratory GA (joint alpha + fraction optimisation)
%% ==========================================================================

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('  Phase 1: Exploratory GA  (%d alpha + %d frac = %d vars)\n', N_wp, N_wp, N_ga);
fprintf('==========================================================================\n');
fprintf('  Nominal waypoints (km): '); fprintf('%.0f  ', p.wp_ranges/1e3); fprintf('\n');
fprintf('  Pop = %d  |  MaxGen = %d  |  StallGen = %d (best, tol=%.0f m)\n', ...
        p.ga_pop_size, p.ga_max_gen, p.ga_stall_gen, p.ga_stall_tol);
fprintf('  Penalty: w_alt=%.1f  w_vel=%.1f\n', p.ga_w_alt, p.ga_w_vel);

lb_ga = [p.wp_lb,  p.frac_lb];
ub_ga = [p.wp_ub,  p.frac_ub];
pop_seed = build_joint_seed_population(p);

stall_gens = p.ga_stall_gen;
stall_tol  = p.ga_stall_tol;
output_fn  = @(opts, state, flag) stall_on_best(opts, state, flag, stall_gens, stall_tol);

ga_opts_p1 = optimoptions('ga', ...
    'PopulationSize',           p.ga_pop_size, ...
    'MaxGenerations',           p.ga_max_gen, ...
    'MaxStallGenerations',      p.ga_max_gen, ...   % disable built-in; use OutputFcn
    'FunctionTolerance',        0, ...              % disable built-in stall tolerance
    'EliteCount',               p.ga_elite_count, ...
    'CrossoverFraction',        p.ga_crossover_frac, ...
    'SelectionFcn',             {@selectiontournament, p.ga_tournament_size}, ...
    'MutationFcn',              @mutationadaptfeasible, ...
    'InitialPopulationMatrix',  pop_seed, ...
    'OutputFcn',                output_fn, ...
    'Display',                  'iter', ...
    'UseParallel',              false);

t_p1 = tic;
% Request final population (6th output) — used to seed Phase 2.
[xga, ~, ~, ~, final_pop_p1, ~] = ga(@(x) ga_objective_joint(x, p), N_ga, ...
                                       [], [], [], [], lb_ga, ub_ga, [], ga_opts_p1);
t_p1 = toc(t_p1);

% Extract GA-optimal alpha and waypoint positions.
wp_ga      = xga(1:N_wp);
fracs_ga   = sort(xga(N_wp+1:end));
p.wp_ranges = fracs_ga * p.ball.x_max;   % update local p for all subsequent calls

traj_ga = simulate_trajectory(wp_ga, p, false);

fprintf('\n  Phase 1 result: %.0f s  |  range = %.1f km  |  feasible = %d\n', ...
        t_p1, traj_ga.x_final/1e3, traj_ga.feasible);
fprintf('  Violations: alt = %+.0f m  |  V = %+.1f m/s\n', ...
        traj_ga.c_viol(1), traj_ga.c_viol(3));
fprintf('  GA-optimal fracs (sorted): '); fprintf('%.3f  ', fracs_ga); fprintf('\n');
fprintf('  GA-optimal waypoints (km): '); fprintf('%.0f  ', p.wp_ranges/1e3); fprintf('\n');
fprintf('  Alpha schedule [deg]:\n    wp_x (km): ');
fprintf('%+6.0f  ', p.wp_ranges/1e3);
fprintf('\n    alpha    : ');
fprintf('%+6.1f  ', wp_ga); fprintf('\n');


%% ==========================================================================
%% Phase 2: Polish GA (alpha only, fixed waypoints, higher penalty pressure)
%% ==========================================================================

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('  Phase 2: Polish GA  (%d alpha vars, waypoints fixed)\n', N_wp);
fprintf('==========================================================================\n');

% Build Phase 2 parameter struct — same as p but with tighter penalties/margins.
p2           = p;
p2.h_max_ga  = p.h_max  - p.ga_polish_alt_margin;   % tighter altitude margin
p2.V_min_ga  = p.V_min_impact + p.ga_polish_vel_margin; % tighter speed margin
p2.ga_w_alt  = p.ga_w_alt * p.ga_polish_penalty_mult;
p2.ga_w_vel  = p.ga_w_vel * p.ga_polish_penalty_mult;

fprintf('  Penalty weights: w_alt=%.1f  w_vel=%.1f  (mult=%.1fx)\n', ...
        p2.ga_w_alt, p2.ga_w_vel, p.ga_polish_penalty_mult);
fprintf('  Constraint margins: alt buffer=%.0f m  |  V buffer=%.0f m/s\n', ...
        p.ga_polish_alt_margin, p.ga_polish_vel_margin);
fprintf('  Pop = %d  |  MaxGen = %d  |  StallGen = %d\n', ...
        p.ga_pop_size, p.ga_polish_max_gen, p.ga_polish_stall_gen);

% Seed Phase 2 from Phase 1 final population (alpha columns only).
% Phase 1's final population already lives in the high-range basin —
% the Polish GA exploits this rather than re-exploring the full space.
alpha_pop_init = final_pop_p1(:, 1:N_wp);   % extract alpha columns

% Ensure first row is the Phase 1 best individual (guaranteed in population).
[~, best_idx]    = min(cellfun(@(x) ga_objective_joint([x, fracs_ga], p), ...
                               num2cell(alpha_pop_init, 2)));
alpha_pop_init([1 best_idx], :) = alpha_pop_init([best_idx 1], :);

stall_gens_p2 = p.ga_polish_stall_gen;
stall_tol_p2  = p.ga_stall_tol;
output_fn_p2  = @(opts, state, flag) stall_on_best(opts, state, flag, stall_gens_p2, stall_tol_p2);

ga_opts_p2 = optimoptions('ga', ...
    'PopulationSize',           p.ga_pop_size, ...
    'MaxGenerations',           p.ga_polish_max_gen, ...
    'MaxStallGenerations',      p.ga_polish_max_gen, ...
    'FunctionTolerance',        0, ...
    'EliteCount',               p.ga_polish_elite_count, ...
    'CrossoverFraction',        p.ga_crossover_frac, ...
    'SelectionFcn',             {@selectiontournament, p.ga_polish_tournament_size}, ...
    'MutationFcn',              @mutationadaptfeasible, ...
    'InitialPopulationMatrix',  alpha_pop_init, ...
    'OutputFcn',                output_fn_p2, ...
    'Display',                  'iter', ...
    'UseParallel',              false);

t_p2 = tic;
[wp_polish, ~] = ga(@(x) ga_objective_alpha_only(x, p2), N_wp, ...
                    [], [], [], [], p.wp_lb, p.wp_ub, [], ga_opts_p2);
t_p2 = toc(t_p2);

traj_polish = simulate_trajectory(wp_polish, p, false);   % evaluate with original p (hard constraints)

fprintf('\n  Phase 2 result: %.0f s  |  range = %.1f km  |  feasible = %d\n', ...
        t_p2, traj_polish.x_final/1e3, traj_polish.feasible);
fprintf('  Violations: alt = %+.0f m  |  V = %+.1f m/s\n', ...
        traj_polish.c_viol(1), traj_polish.c_viol(3));
fprintf('  Alpha schedule [deg]:\n    wp_x (km): ');
fprintf('%+6.0f  ', p.wp_ranges/1e3);
fprintf('\n    alpha    : ');
fprintf('%+6.1f  ', wp_polish); fprintf('\n');

%% Pick best result from Phase 1 and Phase 2 so far ---
% Accept Phase 2 as improvement if it has lower range-equivalent cost
% (weighing feasibility alongside range).
% If neither is feasible, keep whichever has lower violation.
if traj_polish.feasible && traj_ga.feasible
    wp_best   = pick_by_range(wp_ga, traj_ga, wp_polish, traj_polish);
    is_feas   = true;
elseif traj_polish.feasible
    wp_best   = wp_polish;
    is_feas   = true;
elseif traj_ga.feasible
    wp_best   = wp_ga;
    is_feas   = true;
else
    % Neither feasible — keep whichever has lower total violation
    if violation_scalar(traj_polish, p) < violation_scalar(traj_ga, p)
        wp_best = wp_polish;
    else
        wp_best = wp_ga;
    end
    is_feas   = false;
end

traj_best = simulate_trajectory(wp_best, p, false);
fprintf('\n  Best after Phase 2: %.1f km  |  feasible = %d\n', ...
        traj_best.x_final/1e3, traj_best.feasible);

if ~is_feas
    fprintf('\n  *** Neither Phase 1 nor Phase 2 produced a feasible result. ***\n');
    fprintf('  Phase 3 skipped. Returning best result (infeasible).\n');
    fprintf('  Remedies: raise ga_w_alt (%.1f) or ga_polish_penalty_mult (%.1f),\n', ...
            p.ga_w_alt, p.ga_polish_penalty_mult);
    fprintf('            tighten ga_polish_alt_margin / ga_polish_vel_margin.\n\n');
    wp_opt    = wp_best;
    range_opt = traj_best.x_final;
    info      = pack_info(traj_ga, wp_ga, fracs_ga, traj_polish, wp_polish, ...
                          wp_best, [], [], [], [], t_p1, t_p2, false);
    info.wp_ranges_opt = p.wp_ranges;   % CRITICAL: return for main to use
    return
end


%% ==========================================================================
%% Phase 3: fmincon SQP refinement (optional, skip with N_random_starts = 0)
%% ==========================================================================

if N_random_starts <= 0
    fprintf('\n  Phase 3 skipped (N_random_starts = 0).\n');
    wp_opt    = wp_best;
    range_opt = traj_best.x_final;
    info      = pack_info(traj_ga, wp_ga, fracs_ga, traj_polish, wp_polish, ...
                          wp_best, [], [], [], [], t_p1, t_p2, false);
    info.wp_ranges_opt = p.wp_ranges;
    info.total_time    = t_p1 + t_p2;
    return
end

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('  Phase 3: fmincon SQP  (%d start(s), feasible input guaranteed)\n', ...
        1 + N_random_starts);
fprintf('==========================================================================\n');

lb = p.wp_lb;
ub = p.wp_ub;

fmin_opts = optimoptions('fmincon', ...
    'Algorithm',                'sqp', ...
    'MaxIterations',             p.fmin_max_iter, ...
    'MaxFunctionEvaluations',    p.fmin_max_iter * (2*N_wp + 2), ...
    'FunctionTolerance',         p.fmin_tol_fun, ...
    'ConstraintTolerance',       p.fmin_tol_con, ...
    'StepTolerance',             0.01, ...
    'FiniteDifferenceType',      p.fmin_fd_type, ...
    'FiniteDifferenceStepSize',  p.fmin_fd_step, ...
    'Display',                   'off');

obj_fn = @(x) objective_only(x, p);
con_fn = @(x) constraints_normalised(x, p);

best_J  = -traj_best.x_final;
best_wp = wp_best;

all_ranges = zeros(1, 1 + N_random_starts);
all_times  = zeros(1, 1 + N_random_starts);
all_exits  = zeros(1, 1 + N_random_starts);
all_feas   = false(1, 1 + N_random_starts);

fprintf('\n  [Start 1/%d]  Phase-2 result (no perturbation)...\n', 1 + N_random_starts);
t0 = tic;
[wp_k, J_k, exit_k] = run_fmincon(wp_best, obj_fn, con_fn, lb, ub, fmin_opts);
t_k    = toc(t0);
traj_k = simulate_trajectory(wp_k, p, false);
all_ranges(1) = traj_k.x_final; all_times(1) = t_k;
all_exits(1)  = exit_k;         all_feas(1)  = traj_k.feasible;
print_p3_result(traj_k, exit_k, t_k);
if exit_k > 0 && traj_k.feasible && J_k < best_J
    best_J = J_k; best_wp = wp_k;
end

for k = 1:N_random_starts
    noise    = (rand(1, N_wp) - 0.5) * 2 * p.phase3_restart_delta;
    wp0_rand = max(lb, min(ub, wp_best + noise));
    if ~simulate_trajectory(wp0_rand, p, false).feasible
        wp0_rand = wp_best;
    end
    fprintf('  [Start %d/%d]  Perturbed (+/-%.1f deg)...\n', ...
            k+1, 1+N_random_starts, p.phase3_restart_delta);
    t0 = tic;
    [wp_k, J_k, exit_k] = run_fmincon(wp0_rand, obj_fn, con_fn, lb, ub, fmin_opts);
    t_k    = toc(t0);
    traj_k = simulate_trajectory(wp_k, p, false);
    all_ranges(k+1) = traj_k.x_final; all_times(k+1) = t_k;
    all_exits(k+1)  = exit_k;         all_feas(k+1)  = traj_k.feasible;
    print_p3_result(traj_k, exit_k, t_k);
    if exit_k > 0 && traj_k.feasible && J_k < best_J
        best_J = J_k; best_wp = wp_k;
    end
end

wp_opt    = max(lb, min(ub, best_wp));
range_opt = -best_J;

info = pack_info(traj_ga, wp_ga, fracs_ga, traj_polish, wp_polish, ...
                 wp_opt, all_ranges, all_times, all_exits, all_feas, t_p1, t_p2, true);
info.wp_ranges_opt = p.wp_ranges;   % CRITICAL: return for main to use
info.total_time    = t_p1 + t_p2 + sum(all_times);

fprintf('\n  Phase 3 summary: best = %.1f km  |  total = %.0f s (%.1f min)\n', ...
        range_opt/1e3, info.total_time, info.total_time/60);
end


%% ==========================================================================
%% LOCAL FUNCTIONS
%% ==========================================================================

function [state, options, optchanged] = stall_on_best(options, state, flag, sg, st)
% Stop when BEST fitness has not improved by st [m] in the last sg generations.
% state.Best(k) = best objective value in generation k (monotonically non-increasing).
optchanged = false;
if strcmp(flag,'init') || strcmp(flag,'done'), return; end
n = state.Generation;
if n < sg, return; end
improvement = state.Best(n - sg + 1) - state.Best(n);   % >= 0
if improvement < st
    state.StopFlag = sprintf('best stalled: +%.0f m in last %d gen', improvement, sg);
end
end


function pop = build_joint_seed_population(p)
N = p.N_free_wp; a = p.alpha_max_deg; fracs = p.wp_fracs;
% Alpha seeds
s1 = zeros(1, N);
s2 = a * [0.05, 0.08, 0.12, 0.20, 0.30, 0.35, 0.30, 0.22, 0.15, 0.10, 0.05];
s3 = a * [-0.10, -0.10, -0.05, 0.10, 0.25, 0.32, 0.28, 0.20, 0.12, 0.08, 0.03];
s4 = a * [0.00, 0.05, 0.12, 0.28, 0.55, 0.65, 0.55, 0.38, 0.20, 0.12, 0.05];
s5 = (5.0/a) * a * ones(1, N);
alpha_seeds = max(-a, min(a, [s1; s2; s3; s4; s5]));
n_s = size(alpha_seeds, 1);
% Fraction seeds: nominal + tiny noise
frac_seeds = repmat(fracs, n_s, 1) + randn(n_s, N) * 0.02;
frac_seeds = max(repmat(p.frac_lb,n_s,1), min(repmat(p.frac_ub,n_s,1), frac_seeds));
seeds = [alpha_seeds, frac_seeds];
n_rand = max(0, p.ga_pop_size - n_s);
if n_rand > 0
    ra = -a + 2*a * rand(n_rand, N);
    rf = repmat(p.frac_lb,n_rand,1) + rand(n_rand,N) .* repmat(p.frac_ub-p.frac_lb,n_rand,1);
    pop = [seeds; ra, rf];
else
    pop = seeds(1:p.ga_pop_size, :);
end
end


function J = ga_objective_joint(x, p)
% Phase 1 objective: joint alpha + fracs. Fracs sorted internally.
N = p.N_free_wp;
alpha_x     = x(1:N);
fracs       = sort(x(N+1:end));
p.wp_ranges = fracs * p.ball.x_max;
traj = simulate_trajectory(alpha_x, p, false);
pen_alt = max(0, traj.c_viol(1) / p.h_max_ga);
pen_g   = max(0, traj.c_viol(2) / p.n_max);
pen_vel = max(0, traj.c_viol(3) / p.V_min_ga);
J = -traj.x_final + p.ball.x_max * ...
    (p.ga_w_alt * pen_alt^2 + p.ga_w_g * pen_g^2 + p.ga_w_vel * pen_vel^2);
end


function J = ga_objective_alpha_only(alpha, p2)
% Phase 2 objective: alpha only, waypoints fixed in p2.
% Same form as Phase 1 but p2 has tighter margins and higher weights.
% Contains the range term (-x_final) so the degenerate "crash early"
% solution (zero constraint violation at short range) is always dominated
% by any long-range trajectory — the range term is large and negative.
traj    = simulate_trajectory(alpha, p2, false);
pen_alt = max(0, traj.c_viol(1) / p2.h_max_ga);
pen_g   = max(0, traj.c_viol(2) / p2.n_max);
pen_vel = max(0, traj.c_viol(3) / p2.V_min_ga);
J = -traj.x_final + p2.ball.x_max * ...
    (p2.ga_w_alt * pen_alt^2 + p2.ga_w_g * pen_g^2 + p2.ga_w_vel * pen_vel^2);
end


function v = violation_scalar(traj, p)
v = max(0, traj.c_viol(1)/p.h_max)^2 + max(0, traj.c_viol(3)/p.V_min_impact)^2;
end


function wp = pick_by_range(wp_a, traj_a, wp_b, traj_b)
if traj_b.x_final > traj_a.x_final; wp = wp_b; else; wp = wp_a; end
end


function J = objective_only(x, p)
traj = simulate_trajectory(x, p, false);
J    = -traj.x_final;
end


function [c, ceq] = constraints_normalised(x, p)
traj = simulate_trajectory(x, p, false);
c    = [traj.c_viol(1)/p.h_max; traj.c_viol(3)/p.V_min_impact];
ceq  = [];
end


function [wp_out, J_out, exit_out] = run_fmincon(wp0, obj_fn, con_fn, lb, ub, opts)
try
    [wp_out, J_out, exit_out] = fmincon(obj_fn, wp0, [], [], [], [], lb, ub, con_fn, opts);
catch ME
    warning('run_fmincon: %s', ME.message);
    wp_out = wp0; J_out = obj_fn(wp0); exit_out = -99;
end
end


function print_p3_result(traj, exit_k, t_k)
fprintf('    range=%.1f km  feas=%d  exit=%d  %.0f s  |  alt=%+.0f m  V=%+.1f m/s\n', ...
        traj.x_final/1e3, traj.feasible, exit_k, t_k, traj.c_viol(1), traj.c_viol(3));
end


function info = pack_info(traj_ga, wp_ga, fracs_ga, traj_polish, wp_polish, ...
                           wp_best, all_ranges, all_times, all_exits, all_feas, ...
                           t_p1, t_p2, phase3_ran)
info.traj_ga      = traj_ga;
info.wp_ga        = wp_ga;
info.fracs_ga     = fracs_ga;
info.traj_polish  = traj_polish;
info.wp_polish    = wp_polish;
info.wp_best      = wp_best;
info.all_ranges   = all_ranges;
info.all_times    = all_times;
info.all_exits    = all_exits;
info.all_feas     = all_feas;
info.t_p1         = t_p1;
info.t_p2         = t_p2;
info.phase3_ran   = phase3_ran;
info.total_time   = t_p1 + t_p2 + sum(all_times);
end
