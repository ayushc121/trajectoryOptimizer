function generate_aero_tables(filename, csv_filename, S_ref, L_ref)
% GENERATE_AERO_TABLES  Import AoA sweep from FOSTRAD into trajectory format.
%
% CSV columns (no header): [AoA_deg, CL, CD, L_over_D, CMy, Qmax_Wm2].
% The supplied FOSTRAD sweep has no Mach coordinate. CL, CD, and CMy are
% repeated across Mach 2:9 as an explicit modelling assumption, not as
% independent validation at those Mach numbers. Qmax is retained only as
% source metadata: heat flux is NOT Mach-invariant or a trajectory history.
%
% S_ref must equal the opt.SREF used to nondimensionalise FOSTRAD CL/CD.
% L_ref is stored for future moment analysis; CMy is not used by the 2D
% point-mass equations. Both must match the source run before using CMy.

if nargin < 1 || isempty(filename), filename = 'aero_tables.mat'; end
if nargin < 2 || isempty(csv_filename)
    csv_filename = fullfile(fileparts(mfilename('fullpath')), ...
        'owen_tricon_v1_data.csv');
end
if nargin < 3 || isempty(S_ref), S_ref = 0.8; end
if nargin < 4 || isempty(L_ref), L_ref = 2.0; end
if ~isscalar(S_ref) || ~isfinite(S_ref) || S_ref <= 0 || ...
        ~isscalar(L_ref) || ~isfinite(L_ref) || L_ref <= 0
    error('generate_aero_tables:reference', ...
        'Positive finite FOSTRAD reference area and length are required.');
end
if ~exist(csv_filename, 'file')
    error('generate_aero_tables:missingCSV', 'FOSTRAD CSV not found: %s', csv_filename);
end

raw = readmatrix(csv_filename);
if size(raw,2) ~= 6 || size(raw,1) < 2 || any(~isfinite(raw(:)))
    error('generate_aero_tables:format', ...
        'CSV must contain at least two finite rows and exactly six numeric columns.');
end
alpha_vec = raw(:,1).';
CL_alpha  = raw(:,2).';
CD_alpha  = raw(:,3).';
LD_csv    = raw(:,4).';
CMy_alpha = raw(:,5).';
Qmax_Wm2  = raw(:,6).';
if any(diff(alpha_vec) <= 0) || any(CD_alpha <= 0) || any(Qmax_Wm2 < 0)
    error('generate_aero_tables:values', ...
        'AoA must increase strictly, CD must be positive, and Qmax nonnegative.');
end
LD_calc = CL_alpha ./ CD_alpha;
if any(abs(LD_csv - LD_calc) > 1e-8 .* max(1, abs(LD_calc)))
    error('generate_aero_tables:LD', 'CSV L/D disagrees with CL/CD.');
end

% Find an interpolated zero-lift AoA for the reference trajectory. A cambered
% waverider can generate significant lift at zero geometric AoA.
iz = find(CL_alpha(1:end-1) .* CL_alpha(2:end) <= 0, 1, 'first');
if isempty(iz)
    error('generate_aero_tables:zeroLift', ...
        'AoA sweep does not bracket a zero-lift reference condition.');
end
alpha_zero_lift_deg = alpha_vec(iz) - CL_alpha(iz) * ...
    (alpha_vec(iz+1) - alpha_vec(iz)) / ...
    (CL_alpha(iz+1) - CL_alpha(iz));

Mach_vec = 2:9;
nM = numel(Mach_vec);
aero.Mach_vec  = Mach_vec;
aero.alpha_vec = alpha_vec;
aero.CL_table  = repmat(CL_alpha, nM, 1);
aero.CD_table  = repmat(CD_alpha, nM, 1);
aero.CMy_table = repmat(CMy_alpha, nM, 1);  % diagnostic only
aero.LD_source = LD_csv;
aero.Qmax_source_Wm2 = Qmax_Wm2;             % source-run values only
aero.alpha_zero_lift_deg = alpha_zero_lift_deg;
aero.S_ref     = S_ref;
aero.L_ref     = L_ref;
aero.source    = 'FOSTRAD CSV AoA sweep';
aero.source_csv = csv_filename;
aero.mach_invariant_assumed = true;
aero.notes = ['CL/CD/CMy repeated at Mach 2:9; source has no Mach grid. ' ...
    'Qmax is not used in trajectory integration. Trim and moment reference ' ...
    'point have not been validated.'];

save(filename, '-struct', 'aero');
fprintf('Imported %d FOSTRAD AoA rows from %s\n', numel(alpha_vec), csv_filename);
fprintf('AoA: %.2f to %.2f deg | zero-lift AoA: %.3f deg | Sref: %.6g m^2\n', ...
    alpha_vec(1), alpha_vec(end), alpha_zero_lift_deg, S_ref);
fprintf('Replicated CL/CD/CMy over Mach %.1f to %.1f; Qmax retained as source data.\n', ...
    Mach_vec(1), Mach_vec(end));
end
