function T = rotationPeriod(input)
%       Retrieves the rotation period of the planetary body given by the 
% user input or by selecting from the predetermined rotation period of the 
% selected planetary body.  Data extracted from NASA Planetary Fact Sheet 
% at: <https://nssdc.gsfc.nasa.gov/planetary/factsheet/

    % Rotation period
    if isstring(input) % planetary bodies
        if length(input) == 1
            switch input % predetermined inputs
                case "Mercury"
                    sP = 1407.6/23.9345; % [days] rotation period of Mercury
                case "Venus"
                    sP = 5832.6/23.9345; % [days] rotation period of Venus
                case "Earth"
                    sP = 1; % [days] rotation period of Earth
                case "Mars"
                    sP = 24.6229/23.9345; % [days] rotation period of Mars
                case "Jupiter"
                    sP = 9.9250/23.9345; % [days] rotation period of Jupiter
                case "Saturn"
                    sP = 10.656/23.9345; % [days] rotation period of Saturn
                case "Uranus"
                    sP = 17.24/23.9345; % [days] rotation period of Uranus
                case "Neptune"
                    sP = 16.11/23.9345; % [days] rotation period of Neptune
                case "Pluto"
                    sP = 153.2928/23.9345; % [days] rotation period of Pluto
                case "Io"
                    sP = 1.769138; % [days] rotation period of Io
                case "Europa"
                    sP = 3.551181; % [days] rotation period of Europa
                case "Ganymede"
                    sP = 7.154553; % [days] rotation period of Ganymede
                case "Callisto"
                    sP = 16.689017; % [days] rotation period of Callisto
                case "Metis"
                    sP = 0.294779; % [days] rotation period of Metis
                case "Adrastea"
                    sP = 0.298260; % [days] rotation period of Adrastea
                case "Amalthea"
                    sP = 0.498179; % [days] rotation period of Amalthea
                case "Thebe"
                    sP = 0.6745; % [days] rotation period of Thebe
                case "Mimas"
                    sP = 0.9424218; % [days] rotation period of Mimas
                case "Enceladus"
                    sP = 1.370218; % [days] rotation period of Enceladus
                case "Tethys"
                    sP = 1.887802; % [days] rotation period of Tethys
                case "Dione"
                    sP = 2.736915; % [days] rotation period of Dione
                case "Rhea"
                    sP = 4.517500; % [days] rotation period of Rhea
                case "Titan"
                    sP = 15.945421; % [days] rotation period of Titan
                case "Hyperion"
                    sP = 21.276609; % [days] rotation period of Hyperion
                case "Iapetus"
                    sP = 79.330183; % [days] rotation period of Iapetus


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