function temp = getThermalExitance(temp, sData, PMC)

    if ~PMC.exitance.save
        fprintf("\nSkipping Steady-State Thermal Exitance\n\n");
        return
    end
    
    fprintf("\nRunning Photon Monte Carlo: Thermal Exitance\n");
    startTime = tic();
    PMC.exitance.on = true;

    % Get initial photon values
    PMC.thermalPhotons = PMC.exitance.thermalPhotons;
    PMC.emiFunc = @snowEmissivity; % emissivity of snow
    [thermalPhotons, ~] = getThermalPhotons(temp, sData, PMC, "Snow");

    photons.values = thermalPhotons;
    photons.thermalStart = 1;

    % Photon tracing
    traceStart = tic();
    [~, eData] = simulatePhotons(sData, photons, PMC);
    fprintf("\n%30s %f s\n","Photon tracing time:", toc(traceStart));

    % Binning exitance events into exitance angular bins
    binStart = tic();
    spectralRadiance = binExitData(temp, eData.thermal, PMC, "thermal");
    fprintf("\n%30s %f s\n","Binning time:", toc(binStart));
    
    temp.thermalExitance = spectralRadiance;
    fprintf("%30s %f s\n\n","Completion time:", toc(startTime));
end