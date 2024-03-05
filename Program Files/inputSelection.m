function [PMC, HT, SMC, ASM] = inputSelection(model)
%       Selects the essential inputs needed by a specific model and stores 
%   them in shorthand names within the 'in' struct. 'PMC' stands for Photon
%   Monte Carlo, 'HT' for Heat Transfer and 'SMC' for Sublimation Monte
%   Carlo.

    % Photon Monte Carlo (PMC)
    PMC.type = model.type; % simulation type
    PMC.wavelengthBounds = model.light.wavelengthBounds; % [microns]
    PMC.lightSource = model.light.source; % light source
    PMC.density = model.material.density; % [kg/m^3] density
    PMC.grainRadii = model.material.grainRadii; % [microns] grain radii
    PMC.solarPhotons = model.resolution.sourcePhotons; % light photons
    PMC.thermalPhotons = model.resolution.thermalPhotons; % thermal photons
    PMC.emissivity = model.material.emissivity; % emissivity
    PMC.exitance = model.exitance; % exitance parameters
    PMC.exitance.on = false; % 0 during PMC-HT loop, 1 during exitance runs

    % Heat Transfer (HT)
    HT.type = model.type; % simulation type
    HT.species = model.material.species; % species
    HT.rho = model.material.density; % [kg/m^3] density
    HT.I = model.material.thermalInertia; % [J/m^2-K-s^0.5] thermal inertia
    HT.Hmax = model.resolution.Hmax; % [m] max element side length
    HT.Hgrad = model.resolution.Hgrad; % mesh growth rate

    % Sublimation Monte Carlo (SMC)
    SMC.type = model.type; % simulation type
    SMC.gravity = model.gravity; % gravity effect
    SMC.species = model.material.species; % species
    SMC.rho = model.material.density; % [kg/m^3] density
    SMC.planetaryMass = model.body.mass; % [kg] planetary body mass
    SMC.equatorialRadius = model.body.radius; % [m] planetary body radius
    SMC.nParticles = model.resolution.molecules; % molecules per facet
    
    % Analytic Surface Morphology (ASM)
    ASM.type = model.type; % simulation type
    ASM.declare.adjustByDistance = true; % full morphology method
    
end