function mFlux = sublimationRate(tempProfile,SMC) %# codegen
%       Calculates the Number Flux of water molecules based on the
%   saturation vapor pressure of water from NIST and the temperature.

    % Boltzman Constant
    kB = 1.380649e-23; % [J/K]

    % Molecular mass of the species
    Mw = molecularMass(SMC.species); % [kg/kmol]

    % Saturation vapor pressure
    satP = saturationPressure(tempProfile, SMC.species); % [Pa]

    % Number density
    numberDensity = satP./(kB.*tempProfile); % [m^-3]

    % Mean molecular speed
    cBar = ((8*kB.*tempProfile)/(pi*Mw)).^(0.5); % [m/s]

    % Number flux
    nFlux = numberDensity.*cBar/4; % [molecules/m^2-s]
    
    % Mass flux
    mFlux = nFlux * Mw; % [kg/m^2-s]
    fprintf('  - Sublimation rate: %.2e kg/m^2-s\n',abs(sum(mFlux)));

end