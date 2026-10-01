function [CL, CD, CMy] = aero_lookup(M, alpha_deg, aero)
% AERO_LOOKUP  Interpolate CL and CD from the aerodynamic coefficient tables
%
% Inputs:
%   M         - Mach number [-], scalar
%   alpha_deg - angle of attack [deg], scalar
%   aero      - struct from the FOSTRAD CSV import:
%               .Mach_vec   [1 x nM]   Mach breakpoints
%               .alpha_vec  [1 x nA]   AoA breakpoints [deg]
%               .CL_table   [nM x nA]  lift coefficient
%               .CD_table   [nM x nA]  drag coefficient
%
% Outputs:
%   CL - lift coefficient [-]  (positive = upward lift)
%   CD - drag coefficient [-]  (positive = drag, always > 0)
%   CMy - optional pitching-moment coefficient [-]; diagnostic only.
%
% Interpolation: bilinear, so values stay within their neighboring table
% entries. Queries outside the table range are clamped for ODE robustness;
% simulate_trajectory marks any out-of-range trajectory infeasible.
%
% NOTE: Tables indexed as CL_table(iMach, iAlpha).
%       interp2(X,Y,Z,Xq,Yq) convention: X = alpha (columns), Y = Mach (rows).

% ---- Clamp to table bounds (no extrapolation) ----
M_q     = max(aero.Mach_vec(1),  min(aero.Mach_vec(end),  double(M)));
alp_q   = max(aero.alpha_vec(1), min(aero.alpha_vec(end), double(alpha_deg)));

% The supplied FOSTRAD sweep has one AoA axis; its Mach rows are identical.
% Avoid a two-dimensional interpolation at every ODE evaluation in that case.
if isfield(aero, 'mach_invariant_assumed') && aero.mach_invariant_assumed
    CL = interp1(aero.alpha_vec, aero.CL_table(1,:), alp_q, 'linear');
    CD = interp1(aero.alpha_vec, aero.CD_table(1,:), alp_q, 'linear');
else
    CL = interp2(aero.alpha_vec, aero.Mach_vec, aero.CL_table, alp_q, M_q, 'linear');
    CD = interp2(aero.alpha_vec, aero.Mach_vec, aero.CD_table, alp_q, M_q, 'linear');
end

if ~isfinite(CL) || ~isfinite(CD) || CD <= 0
    error('aero_lookup:invalidCoefficients', ...
        'Aerodynamic tables produced nonfinite lift or nonpositive drag.');
end
if nargout > 2
    if ~isfield(aero, 'CMy_table')
        error('aero_lookup:missingMoment', 'CMy_table is not present.');
    end
    if isfield(aero, 'mach_invariant_assumed') && aero.mach_invariant_assumed
        CMy = interp1(aero.alpha_vec, aero.CMy_table(1,:), alp_q, 'linear');
    else
        CMy = interp2(aero.alpha_vec, aero.Mach_vec, aero.CMy_table, ...
            alp_q, M_q, 'linear');
    end
end
end
