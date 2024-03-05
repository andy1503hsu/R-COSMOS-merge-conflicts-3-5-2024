function temp = projectedDisplacement(temp,model,data,marker) %#codegen
%       Function calculates the projected net surface displacement of each 
% surface facet over the global step.
    
    gS = temp.geoStep; % geologic step

    % Length of a world day
    dayLength = rotationPeriod(model.body.rotation); % [s/worldDay]

    % Orbital morphology
    if marker == "orbital"
        oS = temp.orbitStep; % orbital step

        % Calculated facet displacement over a local day
        disp = sum(data.surface.dispFacet(:,:,oS,gS)); % [m/worldDay]

        % Facet displacement rate
        dispFacet = disp/dayLength; % [m/s]

        % Projection time
        projTime = data.orbit.times(oS); % [m/orbitSegment]

    % Geological projection
    elseif marker == "geological"

        % Surface displacement along the orbit
        orbitDisp = zeros(1,model.surface.facets); % [m/orbitSegment]
        for oS = 1:model.resolution.orbitSteps % for each orbit step
            % Daily facet displacement in the current orbit step
            disp = sum(data.surface.dispFacet(:,:,oS,gS)); % [m/worldDay]
    
            % Length of a world day
            dayLength = rotationPeriod(model.body.rotation); % [s/worldDay]

            % Facet displacement rate
            dispRate = disp/dayLength; % [m/s]

            % Displacement of the orbital segment
            segDisp = dispRate*data.orbit.times(oS); % [m/orbitSegment]

            % Surface (added for all orbital segments)
            orbitDisp = orbitDisp + segDisp; % [m/orbit]
        end

        % Net facet displacement over the orbit period
        dispFacet = orbitDisp/data.orbit.period; % [m/s]

        % Total simulation time
        simulationTime = timePeriod(model.time.total); % [s]
    
        % Projection time
        projTime = simulationTime/model.resolution.geologicSteps; % [s]

    end

    % Net surface displacement
    temp.netDisplacement = dispFacet * projTime; % [m]

end