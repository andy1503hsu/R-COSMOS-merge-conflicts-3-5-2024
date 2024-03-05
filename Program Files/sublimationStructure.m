function in = sublimationStructure(temp)
%       Transforms the temporary structure from R-COSMOS into the essential
%   inputs for the Sublimation Monte Carlo (SMC) model. Also preallocates
%   the outputs.

    % Surface coordinates
    in.xS = temp.xS; % x-coordinates
    in.zS = temp.zS; % z-coordinates

    % Surface temperature profile
    in.tempProfile = temp.tempProfile; % [K]

    % Simulation time
    in.time = temp.time; % [s]

    % Number of facets
    nF = length(in.xS)-1;

    % Net molecular outflux
    in.netFlux = zeros(1,nF); % [kg]

    % Surface displacement
    in.dispFacet = zeros(1,nF); % [m]

    % Latent heat of sublimation
    in.latentHeat = zeros(1,nF); % [W/m^2]

end