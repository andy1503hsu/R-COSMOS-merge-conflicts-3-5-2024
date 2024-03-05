function zen = solarPathZenith(tt, dayL)
%       Take in the time and return the solar azimuth angle of the sun in
%   the sky (in radians).
    
    perc = tt/dayL - floor(tt/dayL); % percentage completion of day
    
    zen = pi*abs(1 - 2*perc);
end