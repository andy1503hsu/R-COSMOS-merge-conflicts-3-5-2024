% Add the remaining photon properties that can be universally applied to
% all photons without any conditionals for photon "type" (solar, snow,
% or wall photon).

% Currently, these are:
%   The kill energy threshold
%   The 3 single-scattering properties
%   Maximum number of absorption events
%   Reflectivity when striking annodized aluminium (only Ark Chamber case)

function photons = addRemainingProperties(photons, inPMC)
    
    % If matrix is just empty
    if size(photons, 1) == 0
        return
    end
    
    % Kill energy, or the energy at which a photon's tracing is considered
    % completed and the remaining energy is drained at the same location
    killEnergy = 0.01;  % Default is 1%
    photons(:, 10) = photons(:, 9) * killEnergy;
    
    % Single scattering properties, from the delta-Eddington approximation
    % mu --> attenuation coefficent
    % omega --> albedo
    % g --> asymmetry parameter / anisotropy factor
    photons(:, 11:13) = interp1(inPMC.mieValues(:, 1), ...
                                inPMC.mieValues(:, 2:4), ...
                                photons(:, 8), "linear");
    
    % Maximum number of absorption events
    photons(:, 14) = ceil(log(killEnergy) ./ log(photons(:, 12))) + 1;

    % Reflectivity of photons striking annodized aluminium
    if inPMC.type == "Ark Chamber"
        photons(:, 16) = inPMC.arkWalls.reflectivityFunction(photons(:, 8));
    end
end