function mass = planetaryBodyMass(input) %#codegen
%       Retrieves the planetary body mass input by the user or selects a
% predetermined mass input from the planet selected. Data extracted from
% NASA Planetary Fact Sheet at:
% <https://nssdc.gsfc.nasa.gov/planetary/factsheet/

    % Planetary mass
    if isstring(input) % planets and moons
        switch input % predetermined inputs
            case "Mercury"
                mass = 0.33010e24; % [kg] mass of Mercury
                airlessWarning
            case "Venus"
                mass = 4.8673e24; % [kg] mass of Venus
                airlessWarning
            case "Earth"
                mass = 5.9722e24; % [kg] mass of Earth
                airlessWarning
            case "Mars"
                mass = 0.64169e24; % [kg] mass of Mars
                airlessWarning
            case "Jupiter"
                mass = 1898.13e24; % [kg] mass of Jupiter
                airlessWarning
            case "Saturn"
                mass = 568.32e24; % [kg] mass of Saturn
                airlessWarning
            case "Uranus"
                mass = 86.811e24; % [kg] mass of Uranus
                airlessWarning
            case "Neptune"
                mass = 102.409e24; % [kg] mass of Neptune
                airlessWarning
            case "Pluto"
                mass = 0.01303e24; % [kg] mass of Pluto
            case "Io"
                mass = 893.2e20; % [kg] mass of Io
            case "Europa"
                mass = 480.0e20; % [kg] mass of Europa
            case "Ganymede"
                mass = 1481.9e20; % [kg] mass of Ganymede
            case "Callisto"
                mass = 1075.9e20; % [kg] mass of Callisto
            case "Metis"
                mass = 0.001e20; % [kg] mass of Metis
            case "Adrastea"
                mass = 0.0002e20; % [kg] mass of Adrastea
            case "Amalthea"
                mass = 0.075e20; % [kg] mass of Amalthea
            case "Thebe"
                mass = 0.008e20; % [kg] mass of Thebe
            case "Mimas"
                mass = 0.379e20; % [kg] mass of Mimas
            case "Enceladus"
                mass = 1.08e20; % [kg] mass of Enceladus
            case "Tethys"
                mass = 6.18e20; % [kg] mass of Tethys
            case "Dione"
                mass = 11.0e20; % [kg] mass of Dione
            case "Rhea"
                mass = 23.1e20; % [kg] mass of Rhea
            case "Titan"
                mass = 1345.5e20; % [kg] mass of Titan
                airlessWarning
            case "Hyperion"
                mass = 0.056e20; % [kg] mass of Hyperion
            case "Iapetus"
                mass = 18.1e20; % [kg] mass of Iapetus

            otherwise
                error('Planetary body %s not found.', input)
        end

    elseif isa(input,'double') % planetary body mass input by user
        mass = input; % [kg]

    else
        error('Input type %s is not supported.', class(input))
    end

    function airlessWarning
        warning('Selected planetary body is not an airless world.')
    end
end