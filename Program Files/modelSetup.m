function [model, temp, data, PMC, HT, SMC, ASM] = modelSetup(model)
%       Setup the R-COSMOS simulation environment: time, surface, orbit,
%   light source path. Display a detailed description of the user inputs.

    % Model description
    modelDescription(model)

    % Generate planetary body
    model = planetaryBody(model);

    % Initialize the internal temporary structure
    temp = temporaryStructure(model);

    % Preallocate the output data structure
    data = dataStructure(model);

    % Generate orbital path
    data.orbit = orbitalPath(model);

    % Generate light source path
    data.light = lightPath(data.orbit,model);

    % Retrieve essential inputs for each module
    [PMC, HT, SMC, ASM] = inputSelection(model);

    % Generate surface
    model.surface = surfaceShape(model.surface);
    temp.xS = model.surface.x; % x-coordinates
    temp.zS = model.surface.z; % z-coordinates

    % Create thermal model
    temp = thermalModel(temp,HT);

    % Compile the Monte Carlo modules
    moduleCompiler(temp, PMC, SMC);

    % Retrieve single-scattering coefficients for photon scattering
    PMC.mieValues = mieCalculator(PMC);

end