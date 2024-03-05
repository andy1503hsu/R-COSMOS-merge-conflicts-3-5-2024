function ecc = eccentricity(input)
%       Retrieves the eccentricity of the orbit from the user input or
% selects from the predetermined eccentricities of the selected planet.
%  Data extracted from NASA Planetary Fact Sheet at:
% <https://nssdc.gsfc.nasa.gov/planetary/factsheet/

    % Orbit eccentricity
    if isstring(input) % if input is class string
        switch input % planets predetermined inputs
            case "Mercury"
                ecc = 0.2056; % eccentricity of Mercury
            case "Venus"
                ecc = 0.0068; % eccentricity of Venus
            case "Earth"
                ecc = 0.0167; % eccentricity of Earth
            case "Mars"
                ecc = 0.0935; % eccentricity of Mars
            case "Jupiter"
                ecc = 0.0487; % eccentricity of Jupiter
            case "Saturn"
                ecc = 0.0520; % eccentricity of Saturn
            case "Uranus"
                ecc = 0.0469; % eccentricity of Uranus
            case "Neptune"
                ecc = 0.0097; % eccentricity of Neptune
            case "Pluto"
                ecc = 0.2444; % eccentricity of Pluto

            otherwise
                error('Only planet name is available as a string input.')
        end

    elseif isa(input,'double') % if input is class double
        ecc = input; % eccentricity input by user

    else
        error('Input type %s is not supported.', class(input))
    end
end