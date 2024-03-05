function temp = getSolarExitance(temp, sData, PMC)
    
    if ~PMC.exitance.save
        fprintf("\nSkipping Steady-State Solar Exitance\n\n");
        return
    end

    fprintf("Running Photon Monte Carlo: Solar Exitance \n");
    startTime = tic();

    % Get initial photon values
    PMC.solarPhotons = PMC.exitance.sourcePhotons;  % Hmm...

    % To simplify phase angles, only emit photons from 1 incidence angle
    temp.zen(2) = temp.zen(1);
    temp.azi(2) = temp.azi(1);

    solarPhotons = getSolarPhotons(temp, PMC);

    photons.values = solarPhotons;
    photons.thermalStart = size(solarPhotons, 1) + 1;

    % Photon tracing
    traceStart = tic();
    [~, eData] = simulatePhotons(sData, photons, PMC);
    fprintf("\n%30s %f s\n","Photon tracing time:", toc(traceStart));

    % Binning exitance events into exitance angular bins
    binStart = tic();
    spectralRadiance = binExitData(temp, eData.solar, PMC, "solar");
    fprintf("\n%30s %f s\n","Binning time:", toc(binStart));
    
    temp.solarExitance = spectralRadiance;
    fprintf("%30s %f s\n\n","Completion time:", toc(startTime));

    %temp.solarEmission = sum(photons.values(:, 9));
end