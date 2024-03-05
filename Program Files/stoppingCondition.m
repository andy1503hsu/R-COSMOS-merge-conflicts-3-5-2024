function condition = stoppingCondition(temp, data, model, marker) %#codegen
%       Function determines whether the diurnal cycle of the radiative heat
% transfer process has reached steady-state.

    % Diurnal cycle
    if marker == "steadyState"
        % First condition: steady-state
        % Actual surface temperature
        actual = temp.actualSurfTemp;
    
        % New iteration
        new = data.surface.temperature(:,:,temp.orbitStep,temp.geoStep);
    
        % Relative error in surface temperature
        relativeError = abs(actual-new);

        fprintf("%35s %.3f K\n", "Average relative error:", mean(relativeError, "all"))        
        fprintf("%35s %.3f K\n", "Maximum relative error:", max(relativeError, [], "all"))
        fprintf("%35s %.3f K\n", "Tolerance for max relative error:", model.resolution.tolerance)

        % Stopping condition
        steadyState = all(all(relativeError < model.resolution.tolerance));

        % Second condition: maximum number of iterations
        iterations = temp.iter == model.resolution.maxDays;

        % Stopping condition
        condition = steadyState || iterations;
    
    % Time period
    elseif marker == "timePeriod"
        condition = temp.timePeriod == model.resolution.geologicSteps || model.resolution.geologicSteps == 1;

    % Error condition
    else
        error('Marker not found: select 1) steadyState or 2) timePeriod.')
    end
end
