function temp = sublimationSolver(temp,SMC)
%       Uses the Sublimation Monte Carlo (SMC) solver to simulate
%   sublimation and ensures compatibility with the R-COSMOS architecture.

    % Sublimation Monte Carlo (SMC) inputs
    in = sublimationStructure(temp);

    % Sublimation solver
    in = sublimationMonteCarlo_mex(in,SMC); % SMC

    % Outputs
    temp.netFlux = in.netFlux; % [kg] net molecular outflux
    temp.dispFacet = in.dispFacet; % [m] surface displacement
    temp.latentHeat = in.latentHeat; % [W/m^2] latent heat of sublimation

end