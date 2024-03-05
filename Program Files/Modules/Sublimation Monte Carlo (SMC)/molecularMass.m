function Mw = molecularMass(species) %#codegen
%       Calculates the molecular mass of the species by dividing its molar
% mass by Avogadro's Number.

    % Avogadro's Number
    Na = 6.0221409e26; % [molecule/kmol]

    % Molar mass of the species
    % Water/snow
    if species == "H2O"
        molarMass = 18.01528; % [kg/kmol]

    % Carbon dioxide
    elseif species == "CO2"
        molarMass = 44.01; % [kg/kmol]

    % Methane
    elseif species == "CH4"
        molarMass = 16.04; % [kg/kmol]

    % Ammonia
    elseif species == "NH3"
        molarMass = 17.031; % [kg/kmol]

    % Carbon monoxide
    elseif species == "CO"
        molarMass = 28.01; % [kg/kmol]

    % Nitrogen
    elseif species == "N2"
        molarMass = 28.02; % [kg/kmol]

    else
        error('Species not installed. Select: H2O, CO2, CH4, NH3, CO, N2.')
    end

    % Molecular mass
    Mw = molarMass/Na; % [kg/molecule]

end