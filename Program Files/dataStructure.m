function data = dataStructure(model)
%       Preallocates the data structure to store the simulation data.

    % Dimensions to allocate fields in the data struct
    gS = model.resolution.geologicSteps; % number of global steps
    oS = model.resolution.orbitSteps; % number of orbit steps
    dS = model.resolution.daySteps; % number of day steps
    nF = model.surface.facets; % number of surface facets
    
    % Surface data
    data.surface.x = zeros(oS,nF+1,gS); % x-coordinates
    data.surface.z = zeros(oS,nF+1,gS); % z-coordinates
    data.surface.dispFacet = zeros(dS,nF,oS,gS); % [m] facet displacement
    data.surface.netFlux = zeros(dS,nF,oS,gS); % [kg/m^2-s] net mass flux
    data.surface.temperature = zeros(dS,nF,oS,gS); % [K] facet temperature

    % Steady state temperature
    data.surface.steadyStateT = zeros(dS,nF,model.resolution.maxDays,oS,gS); % [K]

    % Bond albedo
    data.exitance.bondAlbedo = zeros(dS, oS, gS);

    % Exitance data -- only steady state exitance is saved
    data.exitance.thermal = cell(dS, oS, gS);
    data.exitance.solar = cell(dS, oS, gS);

    % Mesh data
    data.mesh.nodes = cell(oS, gS);
    data.mesh.elements = cell(oS, gS);

    % Steady state mesh temperatures and internal heat sources
    data.mesh.temperatures = cell(oS, gS);
    data.mesh.lightAbsorption = cell(oS, gS);
    data.mesh.thermalAbsorption = cell(oS, gS);

end