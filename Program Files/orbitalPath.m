function orbit = orbitalPath(model)
%       Calculates the orbital path and returns the radii in AU, the time
% between two radii, and the latitude adjustment of the light source based
% on the axial tilt.

    % Heliocentric gravitational constant
    GMsun = 1.32712440041279419e20; % [m^3/s^2]

    % True anomaly
    iTa = model.orbit.initialTrueAnomaly; % [deg]
    theta = wrapTo360(linspace(iTa, 360+iTa, model.resolution.orbitSteps+1)); % [deg]

    % Sidereal period
    T = siderealPeriod(model.orbit.period); % [s] orbit period
    orbit.period = T; % [s] orbit period

    % Semimajor axis
    sma = ((T/(2*pi))^2*GMsun)^(1/3); % [m]

    % Eccentricity
    ecc = eccentricity(model.orbit.eccentricity);

    % Heliocentric distance along the orbit
    radii = sma*(1-ecc^2)./(1+ecc*cosd(theta)); % [m] orbit radii
    orbit.radii = radii/149597870700; % [AU] orbit radii

    % Elapsed time between radii
    orbit.times = elapsedTime(radii,theta,T,sma,ecc,GMsun); % [s]

    % Latitude adjustment from axial tilt
    orbit.lFactor = axialTilt(model.orbit.axialTilt)*cosd(theta); % [deg]

end