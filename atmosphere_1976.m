function [rho, T, p_atm, a, mu] = atmosphere_1976(h)
% ATMOSPHERE_1976  US Standard Atmosphere 1976
%
% Inputs:
%   h     - geometric altitude [m], scalar or vector, clamped to [0, 86000] m
%
% Outputs:
%   rho   - density          [kg/m^3]
%   T     - temperature      [K]
%   p_atm - static pressure  [Pa]
%   a     - speed of sound   [m/s]
%   mu    - dynamic viscosity [Pa·s]  (Sutherland's law)
%
% Valid to 86 km. All trajectory altitudes are <= 30 km so this model
% is more than adequate.
%
% Reference: US Standard Atmosphere 1976, NOAA/NASA/USAF

% ---- Constants ----
R_air  = 287.058;    % specific gas constant for air  [J/(kg·K)]
gam    = 1.4;        % ratio of specific heats [-]
g0     = 9.80665;    % standard gravity [m/s^2]
C1     = 1.458e-6;   % Sutherland coefficient [kg/(m·s·K^0.5)]
S_suth = 110.4;      % Sutherland temperature [K]

% ---- Layer table: [h_base(m), T_base(K), lapse_rate(K/m)] ----
layers = [
       0,  288.15, -6.5e-3;   % Troposphere
   11000,  216.65,  0.0;      % Tropopause (isothermal)
   20000,  216.65,  1.0e-3;   % Lower stratosphere
   32000,  228.65,  2.8e-3;   % Upper stratosphere
   47000,  270.65,  0.0;      % Stratopause
   51000,  270.65, -2.8e-3;   % Lower mesosphere
   71000,  214.65, -2.0e-3;   % Upper mesosphere
];

% ---- Base pressures at each layer boundary (Pa) ----
p_base = [101325.0; 22632.1; 5474.89; 868.019; 110.906; 66.9389; 3.95642];

% ---- Convert geometric altitude to geopotential height for layer equations ----
% The standard lapse rates and pressure boundaries are in geopotential metres.
R_earth = 6356766.0;   % nominal Earth radius [m]
h = max(0, min(double(h(:)), 86000));
H = R_earth .* h ./ (R_earth + h);
n = numel(h);
T     = zeros(n,1);
p_atm = zeros(n,1);

for i = 1:n
    idx = find(layers(:,1) <= H(i), 1, 'last');
    dh  = H(i) - layers(idx,1);
    T0  = layers(idx,2);
    L   = layers(idx,3);
    p0  = p_base(idx);

    if abs(L) < 1e-12          % isothermal layer
        T(i)     = T0;
        p_atm(i) = p0 * exp(-g0 * dh / (R_air * T0));
    else                        % gradient layer
        T(i)     = T0 + L * dh;
        p_atm(i) = p0 * (T(i)/T0)^(-g0/(L*R_air));
    end
end

rho = p_atm ./ (R_air .* T);
a   = sqrt(gam .* R_air .* T);
mu  = C1 .* T.^1.5 ./ (T + S_suth);

% Return scalars if input was scalar
if numel(h) == 1
    rho = rho(1); T = T(1); p_atm = p_atm(1); a = a(1); mu = mu(1);
end
end
