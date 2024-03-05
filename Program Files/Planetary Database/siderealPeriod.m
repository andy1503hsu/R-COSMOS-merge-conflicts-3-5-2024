function T = siderealPeriod(input)
%       Retrieves the sidereal period of the orbit given by the user input
% or by selecting from the predetermined sidereal period of the selected
% planet.  Data extracted from NASA Planetary Fact Sheet at:
% <https://nssdc.gsfc.nasa.gov/planetary/factsheet/

    % Sidereal period
    if isstring(input) % planets
        if length(input) == 1
            switch input % predetermined inputs
                case "Mercury"
                    sP = 87.969; % [days]
                case "Venus"
                    sP = 224.701; % [days]
                case "Earth"
                    sP = 365.256; % [days]
                case "Mars"
                    sP = 686.980; % [days]
                case "Jupiter"
                    sP = 4332.589; % [days]
                case "Saturn"
                    sP = 10759.22; % [days]
                case "Uranus"
                    sP = 30685.4; % [days]
                case "Neptune"
                    sP = 60189.0; % [days]
                case "Pluto"
                    sP = 90560; % [days]

                otherwise
                    error('Only planet name is available as a string input.')
            end
            % Calculate the sidereal period
            T = timePeriod([sP, "days"]); % [s]

        elseif length(input) == 2 % generic time type
            T = timePeriod(input); % [s]

        else
            error('Select either: "PlanetName" or [value, "timeUnit"].')
        end

    elseif isa(input,'double') % sidereal period input by user
        T = input; % [s]

    else
        error('Input type %s is not supported.', class(input))
    end
end