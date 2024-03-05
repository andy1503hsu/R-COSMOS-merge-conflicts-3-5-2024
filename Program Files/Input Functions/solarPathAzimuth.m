function azi = solarPathAzimuth(tt, dayL)
%       Take in the time and return the solar azimuth angle of the sun in
%   the sky.
    
    perc = tt/dayL - floor(tt/dayL); % percentage completion of day
    
    if perc <= 0.5
        azi = pi/2*ones(size(tt));
        %azi = zeros(size(tt));
    else
        azi = 3*pi/2*ones(size(tt));
        %azi = pi*ones(size(tt));
    end
    
    % 0, pi azimuth means Sun goes from east to west of penitente
    % pi/2, 3pi/2 means Sun moved parallel to penitente's ridges (no
    % shadowing)
end