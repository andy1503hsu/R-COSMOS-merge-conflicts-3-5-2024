function in = sublimationMonteCarlo(in,SMC) %#codegen
% Sublimation Monte Carlo (SMC) solver version 3.4
%       Simulates sublimation of a species at hard vacuum conditions given
%   an initial surface shape distribution, temperature profile, and
%   material properties.

fprintf('Sublimation Monte Carlo (SMC)\n'); tic

    % Surface discretization
    sData = surfaceData(in);

    % Sublimation rate
    mFlux = sublimationRate(in.tempProfile,SMC);

    % Molecular emission distributions
    lData = monteCarloData(in,sData,SMC);

    % Sublimation Monte Carlo simulation
    loc = molecularPath(sData,lData,SMC);

    % Net molecular flux and latent heat of sublimation
    [in.netFlux, in.latentHeat] = netFlux(in,sData,mFlux,lData,loc,SMC);

    % Surface displacement
    in.dispFacet = surfaceDisplacement(sData,in.netFlux,SMC);

fprintf('Simulation completed in %.2f seconds\n\n',toc);

end