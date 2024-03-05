function tempProfile = surfaceTemperatureProfile(surface, temperature)
%       Creates a temperature profile for the surface.

    switch temperature.type
        case "sinusoidal"
            % Temperature profile shape
            shape = -sin(surface.x/surface.depth*pi*2-pi/2);
        
            % Amplitude
            amplitude = (temperature.trough-temperature.peak)/2;
            mid = (temperature.trough+temperature.peak)/2;
            profile = amplitude.*shape + mid;
        
            % Averaged temperature
            tempProfile = zeros(1,length(surface.x)-1);
            for ii = 1:length(surface.x)-1
                tempProfile(1,ii) = mean(profile(1,ii:ii+1));
            end
    end
end