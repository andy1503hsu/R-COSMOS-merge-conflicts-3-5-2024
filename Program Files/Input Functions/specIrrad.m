function irr = specIrrad(lightSrc, dist)
%       Return the spectral irradiance at given points. First column of irr
%   are the wavelengths. Second column  are the corresponding spectral
%   irradiances. Dist is distance from the Sun in AU.

    if lightSrc == "Sun"
        
        % Old code
        % filename = ['irradSpectrum_',lightSrc,'.mat'];
        % load(filename,'irrad');

        % New code
        load("Smoothed ASTM-E490-00.mat", "irrad")
        
        irrad(:, 2) = irrad(:, 2) / dist ^ 2; % Scale accordingly
        irr = irrad;
    
    elseif lightSrc == "LEDBar"
        filename = ['irradSpectrum_',lightSrc,'.mat'];
        load(filename,'irrad');
        
        irr = irrad;
        amp = 20.05; % increase irrad so that surf. flux is 102 W/m^2
        irr(:,2) = amp*irr(:,2); % final irradiance

    elseif lightSrc == "BlackBodySpectrum"

        % Capture full spectrum here, mainUI wavelength bounds will
        % truncate the distribution using sampleDistribution.m later on

        % 0.1 to 300 microns capture 99.999993% of black body irradiance
        wL = 0.1e-6; wU = 300e-6;
        irr = sunIrradianceEuropa(wL,wU, dist);
    else 
        disp(lightSrc + " is not a valid irradiance source.")
        disp("See specIrradSource.m for valid sources.")
    end
end