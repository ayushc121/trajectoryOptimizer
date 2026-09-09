function generate_aero_tables(filename)
% GENERATE_AERO_TABLES  Placeholder aerodynamic coefficient tables
%
% Generates CL and CD as functions of Mach and angle of attack for a
% generic hypersonic cone-cylinder-fin projectile.
%
% ======================================================================
% !! PLACEHOLDER DATA — MUST BE REPLACED WITH REAL CBAero OUTPUT !!
% ======================================================================
%
% HOW TO REPLACE:
%   Run CBAero across the Mach and AoA grids defined below and export
%   CL_table and CD_table in the same [nMach x nAlpha] format with the
%   same Mach_vec and alpha_vec breakpoints. Save to aero_tables.mat
%   using the same struct field names. No other file needs changing.
%
% ======================================================================
% ASSUMPTIONS (document which ones change when real data arrives):
%
% [A1] REFERENCE AREA & LENGTH
%   S_ref = pi*(0.075)^2 = 0.01767 m^2   (body base area, D_body = 0.15 m)
%   L_ref = 0.15 m                         (body diameter)
%   This convention MUST match whatever CBAero uses. If CBAero uses a
%   different reference area (e.g. cross-section of full fin span), all
%   CL and CD values will be wrong by a constant factor. Verify with
%   the aerodynamics team before running the optimizer.
%
% [A2] TRIMMED COEFFICIENTS
%   Tables represent trimmed (pitch-equilibrium) coefficients. At each
%   (M, alpha) the control surface deflection required to zero pitching
%   moment has already been applied, and the resulting trim drag is
%   included in CD. No CM table is needed for the point-mass trajectory
%   simulation. This matches standard CBAero trimmed-output mode.
%
% [A3] REYNOLDS NUMBER
%   CBAero incorporates skin-friction corrections that implicitly embed
%   Reynolds number effects. The tables here are NOT parameterised by Re.
%   At extreme conditions (sea-level Mach 8 or 30 km Mach 3) this may
%   introduce 5–15% error in CD. This is acceptable for trajectory
%   planning; validate with ANSYS Fluent at the frozen design point.
%
% [A4] AERODYNAMIC MODEL FORM
%   CL(M, alpha) = CLa(M) * alpha          [linear; valid |alpha| < ~12 deg]
%   CD(M, alpha) = CD0(M) + k(M)*CL^2      [standard polar]
%
%   CLa(M): lift-curve slope [1/deg]
%     Inspired by linearised supersonic theory (4/sqrt(M^2-1) per radian)
%     blended with a Newtonian floor at high Mach. Values are calibrated
%     to give L/D_max ~ 2–3 at Mach 5–8, consistent with published data
%     for fin-stabilised slender-body projectiles of similar fineness ratio.
%
%   CD0(M): zero-lift drag coefficient
%     Wave drag (Newtonian impact theory, 5-deg half-angle cone) plus
%     skin-friction contribution. Decreases with Mach beyond Mach 3.
%
%   k(M): induced drag factor
%     Typical values for slender hypersonic bodies. Decreases slightly
%     with Mach as Newtonian pressure distributions become more efficient.
%
% [A5] GEOMETRY ASSUMED
%   • 5-deg half-angle conical nose, ~0.8 m length
%   • Cylinder, D = 0.15 m, ~1.2 m length
%   • 4 planar fins, semi-span 0.125 m (body to tip), root chord ~0.3 m
%   • Total diameter including fins: 0.4 m (meets competition limit)
%
% [A6] MACH RANGE
%   Table spans Mach 2–9 for robustness near the Mach 3 terminal floor.
%   Extrapolation below Mach 1.5 is invalid (transonic, separate model
%   needed). The controller clamps queries to the table range.
%
% [A7] AoA RANGE
%   ±15 deg, matching the competition control authority limit.
%   Linear CL model degrades above ~12 deg; stall is not modelled.
%   The optimizer bounds AoA to ±15 deg, so stall should not be reached.
%
% ======================================================================

fprintf('  Building placeholder aero tables (replace with CBAero data)\n');

%% ---- Reference values [A1] ----
S_ref = pi * (0.075)^2;   % [m^2]  body base area
L_ref = 0.15;              % [m]    body diameter

%% ---- Breakpoint grids ----
Mach_vec  = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];  % [-]
alpha_vec = (-15:1:15);                                    % [deg]

nM = numel(Mach_vec);
nA = numel(alpha_vec);

%% ---- Mach-dependent polar parameters [A4] ----
%
%       M:   2.0   3.0   4.0   5.0   6.0   7.0   8.0   9.0
CLa_M = [0.180 0.140 0.120 0.110 0.100 0.100 0.090 0.090]; % [1/deg]
CD0_M = [0.250 0.200 0.150 0.100 0.080 0.070 0.060 0.060]; % [-]
k_M   = [1.000 0.800 0.750 0.700 0.680 0.670 0.670 0.670]; % [-]

% Interpolate smoothly onto Mach_vec (already defined at those points, but
% kept explicit for readability and future in-filling of additional Mach points)
CLa_interp = interp1(Mach_vec, CLa_M, Mach_vec, 'pchip');
CD0_interp = interp1(Mach_vec, CD0_M, Mach_vec, 'pchip');
k_interp   = interp1(Mach_vec, k_M,   Mach_vec, 'pchip');

%% ---- Build 2-D tables ----
CL_table = zeros(nM, nA);
CD_table = zeros(nM, nA);

for iM = 1:nM
    CL_row = CLa_interp(iM) .* alpha_vec;           % [1 x nA]
    CD_row = CD0_interp(iM) + k_interp(iM) .* CL_row.^2;
    CL_table(iM,:) = CL_row;
    CD_table(iM,:) = CD_row;
end

%% ---- Print derived L/D_max summary ----
LD_max    = 1 ./ (2 * sqrt(k_M .* CD0_M));
CL_opt    = sqrt(CD0_M ./ k_M);
alpha_opt = CL_opt ./ CLa_M;

fprintf('  L/D summary (PLACEHOLDER values):\n');
fprintf('  %6s %8s %8s %8s\n','Mach','L/D_max','CL_opt','a_opt(deg)');
for iM = 1:nM
    fprintf('  %6.1f %8.2f %8.3f %8.1f\n', ...
            Mach_vec(iM), LD_max(iM), CL_opt(iM), alpha_opt(iM));
end
fprintf('\n');

%% ---- Package and save ----
aero.Mach_vec  = Mach_vec;
aero.alpha_vec = alpha_vec;
aero.CL_table  = CL_table;   % [nMach x nAlpha]
aero.CD_table  = CD_table;   % [nMach x nAlpha]
aero.S_ref     = S_ref;
aero.L_ref     = L_ref;
aero.CLa_M     = CLa_M;
aero.CD0_M     = CD0_M;
aero.k_M       = k_M;
aero.LD_max    = LD_max;
aero.notes     = [...
    'PLACEHOLDER — replace with CBAero trimmed output. '...
    'See generate_aero_tables.m for full assumption list. '...
    'Assumptions: [A1] S_ref=body base, [A2] trimmed, '...
    '[A3] Re implicit, [A4] linear-polar model. '...
    'Generated: ' datestr(now)];

save(filename, '-struct', 'aero');
fprintf('  Saved to %s\n', filename);
end
