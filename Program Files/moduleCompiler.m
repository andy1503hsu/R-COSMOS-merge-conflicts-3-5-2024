function moduleCompiler(temp, PMC, SMC)
%       Compiles the individual modules into C and organizes them into the
% corresponding temporary subfolder.

    % Photon Monte Carlo (PMC)
    fprintf('Compiling PMC model (solar photons case)... ')
    sData = surfaceData(temp);
    photons.values = zeros(PMC.solarPhotons/temp.numPartitions, 16);
    photons.thermalStart = 1;
    codegen -d Temporary/PMC planetarySimulateSolar.m -args {sData, photons}
    movefile planetarySimulateSolar_mex.* Temporary/PMC  % Format could be either .mexa64 (pwyll) or .mexw64 (?)

    fprintf('Compiling PMC model (thermal photons case)... ')
    photons.values = zeros((length(temp.xS)-1)*PMC.thermalPhotons/temp.numPartitions, 16);
    photons.thermalStart = 1;
    codegen -d Temporary/BMC planetarySimulateThermal.m -args {sData, photons}
    movefile planetarySimulateThermal_mex.* Temporary/BMC
    
    if PMC.exitance.save
        fprintf('Compiling PMC model (solar exitance case)... ')
        photons.values = zeros(PMC.exitance.sourcePhotons, 16);
        photons.thermalStart = 1;
        codegen -d Temporary/PMC_Exitance planetarySimulateSolarExitance.m -args {sData, photons}
        movefile planetarySimulateSolarExitance_mex.* Temporary/PMC_Exitance

        fprintf('Compiling PMC model (thermal exitance case)... ')
        photons.values = zeros((length(temp.xS)-1)*PMC.exitance.thermalPhotons, 16);
        photons.thermalStart = 1;
        codegen -d Temporary/BMC_Exitance planetarySimulateThermalExitance.m -args {sData, photons}
        movefile planetarySimulateThermalExitance_mex.* Temporary/BMC_Exitance
    end

    % Sublimation Monte Carlo (SMC)
    in = sublimationStructure(temp); %#ok<NASGU> 'in' used by codegen 
    fprintf('Compiling SMC model... ')
    codegen -d Temporary/SMC sublimationMonteCarlo -args {in,SMC}
    movefile sublimationMonteCarlo_mex.* Temporary/SMC
end