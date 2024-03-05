function refreshCompilation(temp, PMC)
%   Refreshes PMC compiled files to account for new partitioning
%   More partitioning --> smaller matrices as inputs

    % Delete old-partitioned compiled folders and code
    rmpath('Temporary/PMC')  % Avoids deleting-on-path warning
    rmpath('Temporary/BMC')  % Avoids deleting-on-path warning
    [~, ~] = rmdir('Temporary/PMC', 's');
    [~, ~] = rmdir('Temporary/BMC', 's');

    % Re-create folders for PMC compiled code
    mkdir Temporary PMC % Photon Monte Carlo
    mkdir Temporary BMC % Blackbody Monte Carlo

    % Photon Monte Carlo (PMC)
    fprintf('Recompiling PMC model (solar photons case) with %d partitions... ', temp.numPartitions)
    sData = surfaceData(temp);
    photons.values = zeros(PMC.solarPhotons/temp.numPartitions, 16);
    photons.thermalStart = 1;
    codegen -d Temporary/PMC planetarySimulateSolar.m -args {sData, photons}
    movefile planetarySimulateSolar_mex.* Temporary/PMC  % Format could be either .mexa64 (pwyll) or .mexw64 (?)

    fprintf('Recompiling PMC model (thermal photons case) with %d partitions... ', temp.numPartitions)
    photons.values = zeros((length(temp.xS)-1)*PMC.thermalPhotons/temp.numPartitions, 16);
    photons.thermalStart = 1;
    codegen -d Temporary/BMC planetarySimulateThermal.m -args {sData, photons}
    movefile planetarySimulateThermal_mex.* Temporary/BMC

    addpath(genpath("Temporary"))
end
