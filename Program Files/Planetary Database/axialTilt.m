function tilt = axialTilt(input)
%       Retrieves the axial tilt input by the user or selects a
% predetermined axial tilt from the planet selected. Data extracted from
% NASA Planetary Fact Sheet at:
% <https://nssdc.gsfc.nasa.gov/planetary/factsheet/

    % Axial tilt
    if isstring(input) % planets
        switch input % predetermined inputs
            case "Mercury"
                tilt = 0.034; % [deg] axial tilt of Mercury
            case "Venus"
                tilt = 2.64; % [deg] axial tilt of Venus
            case "Earth"
                tilt = 23.44; % [deg] axial tilt of Earth
            case "Mars"
                tilt = 25.19; % [deg] axial tilt of Mars
            case "Jupiter"
                tilt = 3.13; % [deg] axial tilt of Jupiter
            case "Saturn"
                tilt = 26.73; % [deg] axial tilt of Saturn
            case "Uranus"
                tilt = 82.23; % [deg] axial tilt of Uranus
            case "Neptune"
                tilt = 28.32; % [deg] axial tilt of Neptune
            case "Pluto"
                tilt = 122.53; % [deg] axial tilt of Pluto

            otherwise
                error('Only planet name is available as a string input.')
        end

    elseif isa(input,'double') % axial tilt input by user
        tilt = input; % [deg]

    else
        error('Input type %s is not supported.', class(input))
    end
end