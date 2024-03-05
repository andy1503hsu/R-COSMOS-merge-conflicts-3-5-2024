                    %% Main User Interface %%
%% Workspace
    clc % command windows text
    clear % workspace variables
    close all % windows
    addpath(genpath("Program Files")); % workspace folders
    addpath(genpath("Visualization")); % visualization folders

%% General
    model.filename = "testing"; % name of output file
    model.type = "Planetary"; % simulation type
    model.gravity = false; % mass conservation from gravity
    model.orbitalMorphology = false; % surface changes along the orbit

%% Geologic Time
    model.time.total = [10, "years"]; % total simulation time

%% Planetary body
    model.body.mass = "Europa"; % [kg] planetary body mass
    model.body.radius = "Europa"; % [m] equatorial radius
    model.body.rotation = "Europa"; % [s] rotation period

%% Heliocentric orbit
    model.orbit.period = "Jupiter"; % [s] sidereal orbit period
    model.orbit.axialTilt = 0; % [deg] obliquity to orbit
    model.orbit.eccentricity = 0; % orbital eccentricity
    model.orbit.initialTrueAnomaly = 0; % [deg] initial true anomaly

%% Light
    model.light.source = "Sun"; % light source
    model.light.wavelengthBounds = [0.2 300]; % [microns] wavelength limits

%% Surface
    % Shape
    model.surface.type = "sinusoidal"; % surface shape
    model.surface.width = 1.00; % [m] shape width
    model.surface.height = 0.00; % [m] shape height
    model.surface.depth = 5.00; % [m] depth below trough
    model.surface.number = 1; % number of surface replicas
    model.surface.facets = 200; % number of facets per shape
    model.surface.divisions = 0; % number of divisions per facet
    model.surface.roughness = 0.00; % [%] surface roughness percentage
    model.surface.modifier = 0; % modifying factor

    % Latitude and orientation
    model.surface.latitude = 0; % [deg] latitude
    model.surface.orientation = "E-W"; % field orientation

%% Material properties
    model.material.species = "H2O"; % species
    model.material.density = 400; % [kg/m^3] density
    model.material.emissivity = 1; % emissivity
    model.material.grainRadii = 1000; % [microns] grain radii
    model.material.thermalInertia = 70; % [J/m^2-K-s^0.5] thermal inertia
    model.material.initialTemperature = 110; % [K] initial temperature

%% Numerical resolution
    % Time steps
    model.resolution.geologicSteps = 1; % number of geologic periods
    model.resolution.orbitSteps = 1; % number of orbital periods
    model.resolution.daySteps = 40; % [steps/day] number of day periods

    % Steady-state
    model.resolution.maxDays = 40; % [day] max days to reach steady-state
    model.resolution.tolerance = 0.25; % [K] temperature stopping tolerance

    % Thermal radiative mesh
    model.resolution.Hmax = 5e-1; % [m] max element side length
    model.resolution.Hgrad = 1.1; % mesh growth rate [1-2]

    % Number of particle bundles
    model.resolution.sourcePhotons = 2e5; % light source photons
    model.resolution.thermalPhotons = 1e3; % thermal photons per facet
    model.resolution.molecules = 1e4; % molecules per facet

%% Exitance
    model.exitance.save = false; % save thermal/solar exitance results
    if model.exitance.save
  
    % Exitance resolution (typically > model.resolution values)
    model.exitance.sourcePhotons = model.resolution.sourcePhotons; % 2e6; % Light source photons
    model.exitance.thermalPhotons = model.resolution.thermalPhotons; % Thermal photons per facet

    % Hemispherical angular grid for trajectory binning
    model.exitance.numZenithBins = 6; % Encompassing [0, pi/2)
    model.exitance.numAzimuthBins = 24; % Encompassing [0, 2pi)

    % Wavelength bounds for exitance detection
    model.exitance.solarBounds = [0 3.5]; % [microns]
    model.exitance.thermalBounds = [0 300]; % [microns]

    % Number of uniformally spaced bins within wavelength bounds
    model.exitance.numSolarBins = 70;
    model.exitance.numThermalBins = 225;
    
    end

%% Simulation
    R_COSMOS(model);

        %% Questions, comments, or improvement suggestions? %%
%% Authors
% Antonio Macias | antonio.macias@jpl.nasa.gov | (832) 806-0343
% Andy Hsu | andy.hsu@jpl.nasa.gov | (408) 707-8852
% Anthony Carreon | anthony.carreon@jpl.nasa.gov | (469) 233-9428

%% The University of Texas at Austin, Aerospace Department
% Dr. David Goldstein | david@oden.utexas.edu | (512) 471-4187
% Dr. Philip Varghese | varghese@mail.utexas.edu | (512) 471-3110
% Dr. Laurence Trafton | lmt@astro.as.utexas.edu | (512) 471-1476
% <https://www.ae.utexas.edu>

%% NASA Jet Propulsion Laboratory, Ocean Worlds Laboratory
% Dr. Kevin Hand | kevin.p.hand@jpl.nasa.gov | (626) 487-5379
% Dr. Daniel Berisford | daniel.berisford@jpl.nasa.gov | (626) 318-3861
% <https://oceanworldslab.jpl.nasa.gov>

%% Copyright
% The University of Texas at Austin, NASA Jet Propulsion Laboratory
% © 2021 Ocean Worlds Laboratory. All rights reserved.