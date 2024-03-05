function temp = temporaryStructure(model)
%       Initializes a temporary structure used to transfer data between
%   modules and to ensure time continuity during the simulation.
    
    % Dimensions to allocate fields in the temp struct
    dS = model.resolution.daySteps; % number of day steps
    nF = model.surface.facets; % number of surface facets
    dayLength = rotationPeriod(model.body.rotation); % length of day

    % Time
    temp.geoStep = 0; % geologic period
    temp.iter = 0; % steady-state local day
    temp.dayStep = 0; % local day step
    temp.orbitStep = 0; % orbital step
    temp.time = dayLength/dS; % [s] module simulation time

    % Orbit
    temp.orbitTime = 0; % [s] elapsed time between radii
    temp.lFactor = 0; % [deg] latitude adjustment

    % Light
    temp.lightDistance = 0; % [AU] light source distance
    temp.zen = [0; 0];  % [deg] current zenith of light location
    temp.azi = [0; 0];  % [deg] current azimuth of light location
    temp.numPartitions = 1; % Number of partitions for PMC simulations

    % Surface
    temp.xS = zeros(1,nF+1); % current x-coordinates
    temp.zS = zeros(1,nF+1); % current z-coordinates
    temp.netFlux = zeros(1,nF); % [kg/m^2-s] net mass flux
    temp.dispFacet = zeros(1,nF); % [m] displacement of each facet
    temp.latentHeat = zeros(1,nF); % [W/m^2] latent heat of sublimation

    % Temperature of each facet [K]
    temp.tempProfile = model.material.initialTemperature.*ones(1,nF);

    % Steady-state temperatures [K]
    temp.actualSurfTemp = model.material.initialTemperature.*ones(dS, nF);
    temp.netDisplacement = zeros(1,nF); % [m] net facet displacement

    % Cached absorption information, used to optimize PMC
    temp.prevSolarAbsorption = zeros(0, dS);
    temp.thermalDistribution = zeros(0, nF);

    % Spectral irradiance values for the light source
    temp.specIrrad = [];

    % Bond albedo (calculated regardless of model.exitance.save)
    temp.bondAlbedo = 0;
    
    % Exitance information
    if model.exitance.save
        nZen = model.exitance.numZenithBins;
        nAzi = model.exitance.numAzimuthBins;
        nSolar = model.exitance.numSolarBins;
        nThermal = model.exitance.numThermalBins;
        temp.solarExitance = zeros(nSolar, nZen, nAzi);
        temp.thermalExitance = zeros(nThermal, nZen, nAzi);
    else
        temp.solarExitance = "No Exitance Saved";
        temp.thermalExitance = "No Exitance Saved";
    end
    
    % Mesh
    temp.meshNodes = [];
    temp.meshElements = [];
    temp.meshTemperatures = model.material.initialTemperature;
    temp.meshLightAbsorption = [];
    temp.meshThermalAbsorption = [];
    temp.thermalEmission = zeros(1,nF);

    % Thermal model used in HT solver
    temp.tModel = [];
    
    % Module execution markers
    temp.lightRadiationExecuted = false;
    temp.blackBodyRadiationExecuted = false;

end