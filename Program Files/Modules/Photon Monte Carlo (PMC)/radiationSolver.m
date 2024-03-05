function temp = radiationSolver(temp, sData, PMC) %#codegen

    fprintf('\nRunning Photon Monte Carlo: %s\n', PMC.type);
    startTime = tic();

    % Get mesh info and preallocate to-be-reused matrices
    temp = initializePMC(temp, sData, PMC);

    %% TODO: Streamline this better
    if PMC.type == "Ark Chamber"
        PMC = CONSTANTS_PMC(PMC);
    end

    % Get values of to-be-simulated photons (if any), as well as
    % the energy emitted by each facet
    [photons, energyPerFacet] = getPhotonValues(temp, sData, PMC);

    % Photon tracing
    traceStart = tic();
    [aData, ~] = simulatePhotons(sData, photons, PMC);
    fprintf("\n%30s %f s\n","Photon tracing time:", toc(traceStart));

    % Binning absorption events into heat transfer mesh
    binStart = tic();
    temp = binAbsorptionEvents(temp, aData, energyPerFacet);
    fprintf("\n%30s %f s\n","Binning time:", toc(binStart));
    
    % Save exitance data -- see "TO FIX" for explanation
    %temp = saveExitanceData(temp, eData, inPMC);
    
    fprintf("%30s %f s\n\n","Completion time:", toc(startTime));
end

%% TO FIX:
    % Ultimately, we want to save both thermal and solar exitance/spectra
    % data.
    % See getThermalExitance() for thermal photon simulation after
    % temperatures have reached a steady state. This has been completed.

    % See getSolarExitance() for solar photon simulation. However, work
    % has NOT been done yet to determine the optimal # of wavelength bins
    % as well as wavelength bounds of the solar spectrum. The # of angular
    % bins will likely be the same as thermal simulations for the sake of
    % consistency, but # of wavelength bins and bounds is independent. More
    % work has to be done here.

    % PROBLEM: When to simulate solar exitance data? During the simulations
    % solar exitance *is* simulated, but given experience with thermal the
    % # of photons needed for resolved exitance >>> # of photons needed for
    % resolved internal heat sources.

    % For thermal, this issue is superseded since thermal exitance is only
    % meaningful for steady-state diurnal temperatures, so a second run of
    % thermal photons is required even without this resolution dilemma
    % (and the second run can simply have a much higher # of thermal
    % photons). 

    % For solar, it may be more prudent to run a very very fine simulation
    % from the get-go and just bite the (potentially signficantly) longer
    % computational time, rather than a low-resolution simulation followed
    % by a high-resolution simulation. Will talk to Antonio about this.