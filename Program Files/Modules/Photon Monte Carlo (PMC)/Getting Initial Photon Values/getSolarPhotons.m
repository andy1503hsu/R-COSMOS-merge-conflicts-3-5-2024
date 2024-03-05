% Creates photon bundles that simulate radiation from the Sun within a specified time interval

% This function accounts for the movement of the sun across the sky during this time interval (use of zenith
% azimuth functions), as well as the changing irradiance (proportional to
% the cosine of the zenith angle). In addition, only wavelengths within a
% certain interval are sampled from (instead of 0 to infinity). After sampling 
% to get the wavelength and launch time of a photon bundle, its initial values are found and set accordingly.
% These are:
    % Wavelength, initial energy, kill energy
    % Initial trajectory (in all three dimensions)
    % Initial location (in all three dimensions)
        % x-coordinate is uniformly sampled along width of domain
        % y-coordinate is uniformly sampled along a unit depth [0, 1)
        % z-coordinate is set to height of the domain
    % Initial region (set to ray-tracing region)

function solarPhotons = getSolarPhotons(temp, PMC) %#codegen
    
    if all(cosd(temp.zen) <= 1e-10)  % Sun is below horizon
        solarPhotons = zeros(0, 16);  % No solar photons during this interval
        return
    end

    % Physical constants
    c = 2.99792458e8; % [m/s] speed of light
    h = 6.62607015e-34; % [J-s] Plank's constant
    
    ww = temp.specIrrad(:,1); % [m] wavelengths
    photonSpecIrrad = temp.specIrrad(:,2).*ww/h/c; % [phot/s/m^2/m] photon spectral irradiance
    photonIrrad = trapz(ww,photonSpecIrrad); % [phot/s/m^2], photon irradiance
    photonPower = photonIrrad*max(temp.xS); % [phot/s], photon "power"

    tSections = linspace(0, temp.time, 1e4); % Time sections of current timestep
    zens = linspace(temp.zen(1), temp.zen(2), 1e4); % Corresponding zenith angles
    
    irradAmp = cosd(zens);  % Solar intensity over time
    irradAmp(irradAmp < 1e-10) = 0; % Ignore negative values
    
    photonsPerSec = photonPower*irradAmp; % [phot/s] photon time distribution
    
    if trapz(tSections,photonsPerSec) == 0 % No solar photons during this time interval
        solarPhotons = zeros(0, 16); % empty matrix
        return;
    end
    
    photonsperWw = trapz(tSections,irradAmp)*max(temp.xS)*photonSpecIrrad; % [phot/m] photon spectral distribution
    [bundleWavelengths, photonsPerBundle] = sampleDistribution(ww, photonsperWw, PMC.wavelengthBounds*1e-6, PMC.solarPhotons); 
    
    % Times when bundles are launched (determines zenith/azimuth angle)
    [launchTime, ~] = sampleDistribution(tSections, photonsPerSec, [0 temp.time], PMC.solarPhotons); 
    
    solarPhotons = zeros(PMC.solarPhotons, 16);
    
    % Wavelength and Energy
    solarPhotons(:, 8) = bundleWavelengths;
    solarPhotons(:, 9) = h*c./bundleWavelengths*photonsPerBundle;
    
    % Get solar angles at specified times
    timeNorm = launchTime / temp.time; % [0, 1]

    %{
    % Using zenith + azimuth to get solar location
    zen = timeNorm * diff(temp.zen) + temp.zen(1);
    azi = timeNorm * diff(temp.azi) + temp.azi(1);
    
    % Trajectory-related values
    solarPhotons(:, 5) = -sind(zen).*cosd(azi); % [-] trajectory vector x comp
    solarPhotons(:, 6) = -sind(zen).*sind(azi); % [-] trajectory vector y comp
    solarPhotons(:, 7) = -cosd(zen); % [-] trajectory vector z comp
    %}

    % From solar path
    x_comp = -(temp.lightPath(1,1) + timeNorm*diff(temp.lightPath(1, :)));
    y_comp = -(temp.lightPath(2,1) + timeNorm*diff(temp.lightPath(2, :)));
    z_comp = -(temp.lightPath(3,1) + timeNorm*diff(temp.lightPath(3, :)));
    len_traj = sqrt(x_comp.^2 + y_comp.^ 2 + z_comp.^2);

    solarPhotons(:, 5) = x_comp ./ len_traj; % [-] trajectory vector x comp
    solarPhotons(:, 6) = y_comp ./ len_traj; % [-] trajectory vector y comp
    solarPhotons(:, 7) = z_comp ./ len_traj; % [-] trajectory vector z comp
    
    % Located-related values
    solarPhotons(:, 1) = rand(PMC.solarPhotons, 1) * max(temp.xS); % Uniformly random along width of domain
    solarPhotons(:, 2) = rand(PMC.solarPhotons, 1); % Uniformly random along unit depth
    solarPhotons(:, 3) = max(temp.zS) * ones(PMC.solarPhotons, 1); % Ceiling of domain

    if any(abs(diff(temp.zS)) > 1e-10)
        solarPhotons(:, 4) = ones(PMC.solarPhotons, 1);
        solarPhotons(:, 15) = zeros(PMC.solarPhotons, 1); % Typically, solar bundles won't start on a surface
    else
        solarPhotons(:, 4) = 2*ones(PMC.solarPhotons, 1); % Flat snow field -- solar bundles start in absorptive region
        solarPhotons(:, 15) = ones(PMC.solarPhotons, 1); % Flat snow field -- solar bundles start on surface
    end

    % Single-scattering coefficients and energy-related values
    solarPhotons = addRemainingProperties(solarPhotons, PMC);
end