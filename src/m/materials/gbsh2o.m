function rigidity = gbsh2o(temperature, phi)

if nargin < 2
    phi = 0.01; % default value for phi
end 

if any(temperature < 0)
    error('input temperature should be in Kelvin (positive)');
  end

% Beyond-melting-point case
if any(temperature >= 273.15)
    warning('gbsh2o.m: H2O ICE - GUARANTEED MELTING. Some temperatures are beyond 273.15K.\nLook at indices: %s', mat2str(find(temperature>=273.15))');
  end

% Coefficients
Rg = 8.3144598;    % J mol^-1 K^-1
n = 1.8;           % Glen's exponent

rd = 1.5e-6;
f = 0.25;
A0 = 3.9e-3;       % s^-1 MPa
p = -1.4;

d = (13.4e-6 / sqrt(f)) * (rd / 1.5e-6) .* ((0.05 ./ phi)^0.5);
A1 = A0 * 0.5 * (3^((n + 1) / 2));
A2 = A1 .* exp(-2 * n * phi) .* (d^p);

Q = 49000.;        % J mol^-1


% Arrhenius Law
A = A2 .* exp(-Q ./ (temperature * Rg));  % s^-1 MPa
rigidity = 1e6 .* (A^(-1/n));                   % s^(1/n) Pa

end
