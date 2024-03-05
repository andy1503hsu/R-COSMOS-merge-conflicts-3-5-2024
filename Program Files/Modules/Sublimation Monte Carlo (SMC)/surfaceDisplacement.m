function displacement = surfaceDisplacement(sData,netFlux,SMC)
%       Calculates the surface displacement based on the net molecular flux
%   and the material properties.

    % Surface displacement
    displacement = netFlux ./sData.areaFacet ./ SMC.rho; % [m]
    fprintf('  - Average displacement: %.2e m\n',mean(displacement));
    
end