function [nF,sH] = netFlux(in,sData,mFlux,lData,location,SMC) %#codegen
%       Calculates the net mass flux of water molecules based on their
% current location and the total mass lost from the system.
    
    % Molecular velocity magnitude (speed)
    vL = lData.vL; % [m/s]

    % Launch facet index
    iL = lData.iL;

    % Total mass ejected from each facet
    mass = mFlux .* sData.areaFacet * in.time; % [kg]
    
    % Net mass flux (defined negative for outflux)
    nF = -mass; % [kg] initial mass desorbed

    % Sublimation heat (latent heat of sublimation)
    sH = zeros(1,length(mass)); % [J]

    % Net molecular flux
    for ii = 1:length(location) % for all molecular bundles launched
        
        % Facet from which the molecular bundle was originally launched
        ind = iL(ii);
        
        % Mass of the molecular bundle
        massBundle = mass(ind)/SMC.nParticles; % [kg/bundle]
        
        % Kinetic energy of the molecular bundle
        KE = 1/2*massBundle*vL(ii).^2; % [J]
        
        % Energy lost from launch facet
        sH(ind) = sH(ind) - KE; % [J]

        % Final location of the molecular bundle
        loc = location(ii);
        
        % Net molecular flux
        if loc > 0 % if molecular bundle impacted a surface region
            nF(loc) = nF(loc) + massBundle; % [kg] mass adsorbed
            sH(loc) = sH(loc) + KE; % [J] energy deposited
        end
    end

    % Latent heat of sublimation
    sH = sH./sData.areaFacet/in.time; % [W/m^2]
end