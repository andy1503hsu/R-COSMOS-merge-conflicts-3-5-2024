function ledPhotons = getLEDPhotons(temp, PMC)
    
    % Physical constants
    c = 2.99792458e8; % [m/s] speed of light
    h = 6.62607015e-34; % [J-s] Plank's constant

    % Ark constants
    arkWidth = PMC.arkDimensions(1);
    arkHeight = PMC.arkDimensions(2);
    ledWidth = PMC.LED_Dimensions(1);
    led_deltaz = PMC.LED_Dimensions(2);

    ww = PMC.solarIrrad(:,1); % [m] wavelengths
    photonIrrad = PMC.solarIrrad(:,2).*ww/h/c; % [phot/s/m^2/m] photon spectral irradiance
    photonsPerWw = photonIrrad*diff(temp.tInt)*PMC.LED_Dimensions(1); 
    
    if trapz(ww, photonsPerWw) == 0
        ledPhotons = zeros(0, 16);
        return;
    end
    
    % Sample and get wavelengths of each bundle, along with how many photons each bundle represents
    [bundleWavelengths, photonsPerBundle] = sampleDistribution(ww, photonsPerWw, PMC.wavelengthBounds*1e-6, PMC.wallPhotons); 
    
    ledPhotons = zeros(PMC.wallPhotons, 16);
    
    % Wavelength and Energy
    ledPhotons(:, 8) = bundleWavelengths;
    ledPhotons(:, 9) = h*c./bundleWavelengths*photonsPerBundle;
    
    % Assume that half of photons are emitted from NIR LEDS, while other
    % half are emitted from visible LEDs
    zen_oneHalf = sampleZenithAngles(floor(PMC.wallPhotons/2), "NIR");
    zen_secondHalf = sampleZenithAngles(PMC.wallPhotons - length(zen_oneHalf), "Visible");
    zen = [zen_oneHalf; zen_secondHalf] * pi/180; % Convert to radians
    azi = 2*pi*rand(PMC.wallPhotons,1); % [rad] azimuth angle
    
    % Trajectory-related values
    ledPhotons(:, 6) = -sin(zen).*cos(azi); % [-] trajectory vector x comp
    ledPhotons(:, 7) = -sin(zen).*sin(azi); % [-] trajectory vector y comp
    ledPhotons(:, 8) = -cos(zen); % [-] trajectory vector z comp

    % Location-related values
    ledPhotons(:, 1) = rand(PMC.wallPhotons, 1)*ledWidth + (arkWidth-ledWidth)/2; % [m] Uniformly sampled along x-length of bar
    ledPhotons(:, 2) = rand(PMC.wallPhotons, 1); % [m] Uniformly sampled along unit depth
    ledPhotons(:, 3) = (arkHeight - led_deltaz)*ones(PMC.wallPhotons,1); % [m] starting z coordinates
    ledPhotons(:, 4) = ones(PMC.wallPhotons,1); % All bundles start in ray-tracing region
    ledPhotons(:, 15) = zeros(PMC.wallPhotons,1); % Photons do not start on a surface

    % Times when bundles are launched (determines zenith/azimuth angle)
    ledPhotons = addRemainingProperties(ledPhotons, PMC);
end