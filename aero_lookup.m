function [CL, CD] = aero_lookup(M, alpha_deg, aero)
% AERO_LOOKUP  Interpolate CL and CD from the aerodynamic coefficient tables
%
% Inputs:
%   M         - Mach number [-], scalar
%   alpha_deg - angle of attack [deg], scalar
%   aero      - struct from generate_aero_tables / CBAero import:
%               .Mach_vec   [1 x nM]   Mach breakpoints
%               .alpha_vec  [1 x nA]   AoA breakpoints [deg]
%               .CL_table   [nM x nA]  lift coefficient
%               .CD_table   [nM x nA]  drag coefficient
%
% Outputs:
%   CL - lift coefficient [-]  (positive = upward lift)
%   CD - drag coefficient [-]  (positive = drag, always > 0)
%
% Interpolation: bicubic spline (interp2 'spline') for smooth gradients.
% Queries outside the table range are clamped (no extrapolation).
%
% NOTE: Tables indexed as CL_table(iMach, iAlpha).
%       interp2(X,Y,Z,Xq,Yq) convention: X = alpha (columns), Y = Mach (rows).

% ---- Clamp to table bounds (no extrapolation) ----
M_q     = max(aero.Mach_vec(1),  min(aero.Mach_vec(end),  double(M)));
alp_q   = max(aero.alpha_vec(1), min(aero.alpha_vec(end), double(alpha_deg)));

% ---- Bicubic spline interpolation ----
CL = interp2(aero.alpha_vec, aero.Mach_vec, aero.CL_table, alp_q, M_q, 'spline');
CD = interp2(aero.alpha_vec, aero.Mach_vec, aero.CD_table, alp_q, M_q, 'spline');

% ---- Physical floor on CD ----
CD = max(CD, 1e-6);
end
