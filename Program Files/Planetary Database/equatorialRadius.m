function radius = planetaryBodyRadius(input) %#codegen
%       Retrieves the planet radius input by the user or selects a
% predetermined mass input from the planet selected. Data extracted from
% NASA Planetary Fact Sheet at:
% <https://nssdc.gsfc.nasa.gov/planetary/factsheet/

    % Planetary radius
    if isstring(input) % planets and moons
        switch input % predetermined inputs
            case "Mercury"
                radius = 2440.5e3; % [m] radius of Mercury
                airlessWarning
            case "Venus"
                radius = 6051.8e3; % [m] radius of Venus
                airlessWarning
            case "Earth"
                radius = 6378.1e3; % [m] radius of Earth
                airlessWarning
            case "Mars"
                radius = 3396.2e3; % [m] radius of Mars
                airlessWarning
            case "Jupiter"
                radius = 71492e3; % [m] radius of Jupiter
                airlessWarning
            case "Saturn"
                radius = 60268e3; % [m] radius of Saturn
                airlessWarning
            case "Uranus"
                radius = 25559e3; % [m] radius of Uranus
                airlessWarning
            case "Neptune"
                radius = 24764e3; % [m] radius of Neptune
                airlessWarning
            case "Pluto"
                radius = 1188e3; % [m] radius of Pluto
            case "Io"
                radius = 1821.5e3; % [m] radius of Io
            case "Europa"
                radius = 1560.8e3; % [m] radius of Europa
            case "Ganymede"
                radius = 2631.2e3; % [m] radius of Ganymede
            case "Callisto"
                radius = 2410.3e3; % [m] radius of Callisto
            case "Metis"
                radius = mean([30, 20, 17])*1e3; % [m] radius of Metis
            case "Adrastea"
                radius = mean([10, 8, 7])*1e3; % [m] radius of Adrastea
            case "Amalthea"
                radius = mean([125, 73, 64])*1e3; % [m] radius of Amalthea
            case "Thebe"
                radius = mean([58, 49, 42])*1e3; % [m] radius of Thebe
            case "Mimas"
                radius = mean([208, 197, 191])*1e3; % [m] radius of Mimas
            case "Enceladus"
                radius = mean([257, 251, 248])*1e3; % [m] radius of Enceladus
            case "Tethys"
                radius = mean([538, 528, 526])*1e3; % [m] radius of Tethys
            case "Dione"
                radius = mean([563, 561, 560])*1e3; % [m] radius of Dione
            case "Rhea"
                radius = mean([765, 763, 762])*1e3; % [m] radius of Rhea
            case "Titan"
                radius = 2575e3; % [m] radius of Titan
                airlessWarning
            case "Hyperion"
                radius = mean([180, 133, 103])*1e3; % [m] radius of Hyperion
            case "Iapetus"
                radius = mean([746, 746, 712])*1e3; % [m] radius of Iapetus

            otherwise
                error('Planetary body %s not found.', input)
        end

    elseif isa(input,'double') % planetary body radius input by user
        radius = input; % [m]

    else
        error('Input type %s is not supported.', class(input))
    end

    function airlessWarning
        warning('Selected planetary body is not an airless world.')
    end
end