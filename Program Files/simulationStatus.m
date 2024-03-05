function simulationStatus(temp, model, marker) %#codegen
%       Function displays the status of the simulation.

% Header
fprintf("\n================================\n")
% fprintf("\n--------------------------------\n")

% Time period
fprintf("|| Time Period: %3d of %-3d\n", temp.geoStep,model.resolution.geologicSteps);

% Radiative heat transfer
if marker == "RHT"
    fprintf("||   Iteration: %3d of %-3d (max)\n",temp.iter,model.resolution.maxDays);
end

% Day steps
fprintf("||        Step: %3d of %-3d\n",temp.dayStep,model.resolution.daySteps);

end