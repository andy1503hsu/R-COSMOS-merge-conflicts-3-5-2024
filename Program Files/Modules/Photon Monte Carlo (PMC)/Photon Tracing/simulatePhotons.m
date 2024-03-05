function [aData, eData] = simulatePhotons(sData, photons, PMC) %#codegen
    
    photonValues = photons.values;
    thermalStart = photons.thermalStart;
    
    numPhotons = size(photonValues, 1);
    energyEmitted = sum(photonValues(:, 9));
     
    % No photons will be traced. This occurs from the 2nd day onwards.
    if energyEmitted == 0
       fprintf("\n%30s", "No photons simulated.")
       % Dummy empty variables
       aData.solar = zeros(0, 3); aData.thermal = cell(0, 1);
       eData.solar = zeros(0, 4); eData.thermal = zeros(0, 3);
       return
    else % Consequently, all code past here assumes a 1st day timestep.
        snowEmitted = sum(photonValues(thermalStart:end, 9));
        sunEmitted = sum(photonValues(1:thermalStart-1, 9));
    end

    % Display total energy to be simulated
    fprintf("\n%30s %11.4f J\n", "Total energy emitted:", energyEmitted)
    
    % Display amount of solar energy to be simulated
    fprintf("%30s %11.4f J\n", "Energy from Sun:", sunEmitted)

    % Display amount of thermal energy to be simulated
    if numPhotons - thermalStart + 1 == 0
       fprintf("%30s %11.4f J (reused)\n", "Energy from Snow:", 0)
    else
       fprintf("%30s %11.4f J\n", "Energy from Snow:", snowEmitted)
    end
     
    if PMC.type == "Planetary"
        
        solarPhotons.values = photons.values(1:thermalStart-1, :);
        solarPhotons.thermalStart = photons.thermalStart;
        
        thermalPhotons.values = photons.values(thermalStart:end, :);
        thermalPhotons.thermalStart = 1;

        if ~isempty(solarPhotons.values) % Solar radiation is present
            if PMC.exitance.save && PMC.solarPhotons == PMC.exitance.sourcePhotons % Exitance case
                [aDataSolar, eDataSolar, energyDataSolar] = planetarySimulateSolarExitance_mex(sData, solarPhotons);
            else % Non-exitance case (PMC-HT loop)
                [aDataSolar, eDataSolar, energyDataSolar] = planetarySimulateSolar_mex(sData, solarPhotons);
            end
        else % Solar radiation wasn't simulated
            aDataSolar.solar = zeros(0, 3);
            eDataSolar.solar = zeros(0, 4);
            energyDataSolar.absorbed = 0;
            energyDataSolar.reflected = 0;
            energyDataSolar.exited = 0;
        end

        if ~isempty(thermalPhotons.values) % Thermal radiation is present
            if PMC.exitance.save && PMC.thermalPhotons == PMC.exitance.thermalPhotons % Exitance case
                [aDataThermal, eDataThermal, energyDataThermal] = planetarySimulateThermalExitance_mex(sData, thermalPhotons);
            else % Non-exitance case (PMC-HT loop)
                [aDataThermal, eDataThermal, energyDataThermal] = planetarySimulateThermal_mex(sData, thermalPhotons);
            end
        else % Thermal radiation wasn't simulated
            aDataThermal.thermal = cell(0, 1);
            eDataThermal.thermal = zeros(0, 4);
            energyDataThermal.absorbed = 0;
            energyDataThermal.reflected = 0;
            energyDataThermal.exited = 0;
        end

        aData.solar = aDataSolar.solar; aData.thermal = aDataThermal.thermal;
        eData.solar = eDataSolar.solar; eData.thermal = eDataThermal.thermal;

        energyData.absorbed = energyDataSolar.absorbed + energyDataThermal.absorbed;
        energyData.reflected = energyDataSolar.reflected + energyDataThermal.reflected;
        energyData.exited = energyDataSolar.exited + energyDataThermal.exited;

    elseif PMC.type == "Ark Chamber"
        [aData, energyData] = arkSimulate(sData, photons);
         % Unfinished section
    else
        disp(PMC.type + " is not a valid PMC simulation type.")
    end

    absorbed = energyData.absorbed;
    reflected = energyData.reflected;
    exited = energyData.exited;

    fprintf("\n%30s %11.4f J\n", "Energy absorbed by snow:", absorbed)
    fprintf("%30s %11.4f J\n", "Energy absorbed by surfaces:", reflected)
    fprintf("%30s %11.4f J\n", "Energy that exited domain:", exited)
    
    % Check for energy conservation
    if abs(energyEmitted / ...
           (absorbed+reflected+exited)) -1 < 1e-10
        fprintf("%30s\n", "Energy was conserved.")
    else
        fprintf("%30s\n", "Energy was NOT conserved.")
    end
end

