function satP = saturationPressure(T, material) %#codegen
%       Calculates the saturation vapor pressure of water above ice at a
%   given temperature.

% For water
if material == "H2O" % IAPWS Formulation
    %       Data retrieved from NIST.gov:
    %   Bielska et al. (2013), High‐accuracy measurements of the vapor
    % pressure of ice referenced to the triple point,
    % doi:10.1002/2013GL058474.
    
    % Constants
        a1 = -0.212144006e2;    a2 = 0.273203819e2;    a3 = -0.610598130e1;
        b1 = 0.333333333e-2;    b2 = 0.120666667e1;    b3 = 0.170333333e1;

    % Reference Temperature and Pressure to the triple point
    refT = 273.16; % [K]
    refP = 611.657; % [Pa]

    % Ratio of temperature profile to the reference temperature
    ratio = T./refT;

    % Saturation Vapor Pressure
    P = ratio.^(-1).*(a1*ratio.^(b1) + a2*ratio.^(b2) + a3*ratio.^(b3));
    satP = refP*exp(P); % [Pa]

% For Carbon Dioxide
elseif material == "CO2" % IAPWS Formulation

    % Constants
    a1 = -14.740846;    a2 = 2.4327015;     a3 = -5.3061778;
    b1 = 1;             b2 = 1.9;           b3 = 2.9;
    
    % Reference Temperature and Pressure to the triple point
    refT = 216.592; % [K]
    refP = 0.51795e6; % [Pa]

    % Ratio of temperature profile to the reference temperature
    ratio = T./refT;

    P = ratio.^(-1).*(a1*(1-ratio).^b1 + a2*(1-ratio).^b2 + a3*(1-ratio).^b3);
    satP = refP*exp(P);

% For Carbon Monoxide
elseif material == "CO" % Empirical Vapor Pressure Relations 

    % Temperature validity limits
    if any(T < 14)
        error('Temperature is outside the lower valid bound.')
    elseif any(T > 68.127)
        error('Temperature is outside the upper valid bounds.')
    end

    % Polynomial extrapolation coefficients
    a0 =  1.80741183e1;  % [~]
    a1 = -7.69842078e2;  % [K]
    a2 = -1.21487759e4;  % [K^2]
    a3 =  2.73500950e5;  % [K^3]
    a4 = -2.90874670e6;  % [K^4]
    a5 =  1.20319418e7;  % [K^5]

    c = zeros(1,length(T));

    Tless = find(T < 61.544);
    Tmore = find(T >= 61.544);

    temp = T;
    T_ = temp(Tless);
    c(Tless) = a0 + a1./T_ + a2./T_.^2 + a3./T_.^3 + a4./T_.^4 + a5./T_.^5;
    
    % Polynomial extrapolation coefficients
    a0 =  1.68655152e1;  % [~]
    a1 = -7.48151471e2;  % [K]
    a2 = -5.84330795e3;  % [K^2]
    a3 =  3.93853859e4;  % [K^3]

    T_ = temp(Tmore);
    c(Tmore) = a0 + a1./T_ + a2./T_.^2 + a3./T_.^3;
    satP = 133.322*exp(c);

% For Methane
elseif material == "CH4"
    % Polynomial extrapolation coefficients
    a0 =  1.051e1;  % [~]
    a1 = -1.110e3;  % [K]
    a2 = -4.341e3;  % [K^2]
    a3 =  1.035e5;  % [K^3]
    a4 = -7.910e5;  % [K^4]

    c = a0 + a1./T + a2./T.^2 + a3./T.^3 + a4./T.^4;
    satP = exp(c);

% For Nitrogen
elseif material == "N2"
    % Temperature validity limits
    if any(T < 10)
        error('Temperature is outside the lower valid bound.')
    end
    if any(T > 63.14)
        error('Temperature is outside the upper valid bounds.')
    end

    % Polynomial extrapolation coefficients
    a0 =  1.240e1;  % [~]
    a1 = -8.074e2;  % [K]
    a2 = -3.926e3;  % [K^2]
    a3 =  6.297e4;  % [K^3]
    a4 = -4.633e5;  % [K^4]
    a5 =  1.325e6;  % [K^5]

    c = zeros(1,length(T));

    Tless = find(T < 35.61);
    Tmore = find(T >= 35.61);

    temp = T;
    T = temp(Tless);
    c(Tless) = a0 + a1./T + a2./T.^2 + a3./T.^3 + a4./T.^4 + a5./T.^5;
    
    % Polynomial extrapolation coefficients
    a0 =  8.514;  % [~]
    a1 = -4.584e2;  % [K]
    a2 = -1.987e4;  % [K^2]
    a3 =  4.800e5;  % [K^3]
    a4 =  -4.524e6;  % [K^4]

    T = temp(Tmore);
    c(Tmore) = a0 + a1./T + a2./T.^2 + a3./T.^3 + a4./T.^4;
    satP = exp(c);

% For Ammonia
elseif material == "NH3"
    % Temperature validity limits
    if any(T < 15)
        error('Temperature is outside the lower valid bound.')
    end
    if any(T > 195.41)
        error('Temperature is outside the upper valid bounds.')
    end

    % Polynomial extrapolation coefficients
    a0 =  1.596e1;  % [~]
    a1 = -3.537e3;  % [K]
    a2 = -3.310e4;  % [K^2]
    a3 =  1.742e6;  % [K^3]
    a4 = -2.995e7;  % [K^4]

    c = a0 + a1./T + a2./T.^2 + a3./T.^3 + a4./T.^4;
    satP = exp(c);

else
    error('Species selected is not installed.')
end
end