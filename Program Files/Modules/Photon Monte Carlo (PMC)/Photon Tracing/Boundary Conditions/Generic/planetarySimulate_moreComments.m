function [aData, eData, energyData] = planetarySimulate3(sData, photons) %#codegen
    
    photonValues = photons.values;
    thermalStart = photons.thermalStart;
    numPhotons = size(photonValues, 1);
    
    aEvents = cell(numPhotons, 1);
    parfor index = 1:numPhotons  % Cell array elements must be assigned
        aEvents{index} = index;
    end

    numAEvents = zeros(numPhotons, 1);
    energyThroughBottom = 0;
    eDataAll = zeros(numPhotons, 4); % wavelength, energy, zenith, azimuth
    
    % Cached info about domain/snow surfaces
    [xi,xf,zi,zf, mF, ~,zM, uC] = commonVariables(sData);
    [x34, z34, determinant34] = cachedInfo(xi, xf, zi, zf);
    maxX = max([xi xf(end)]);
    minZ = min([zi zf(end)]);
    zExit = max([zi zf(end)])*1.00;
    maxSlope = max(mF);
    minSlope = min(mF);
    
    % Extract and group variables for efficiency
    gs = photonValues(:, 13); % scattering asymmetry parameter
    g1s = (1+gs.^2)*0.5./gs;
    g2s = 1-gs.^2;
    g3s = 0.5./gs; % grouping of terms for computational efficiency

    pie2 = 6.28318530718;

    %disp("For loop used instead of parfor in order to profile") % Comment in/out and change for/parfor accordingly
    %disp("Simulating " + numPhotons + " photons")

    %{
    numWhileLoopIter = 0;
    numIntersectAlg = 0;
    numIntersectInf = 0;
    numTunneling = 0;
    numNoIntersectShortcuts = 0;
    %}

    parfor index = 1:numPhotons
        %if mod(index, 5000) == 0
        %    disp("Current photon: " + index)
        %end
        % 1 2 3 -- x, y, z locations
        % 4 -- region (1 is ray-tracing, 2 is absorptive)
        % 5 6 7 -- x, y, z trajectory
        % 8 -- wavelength
        % 9 10 -- energy and kill energy
        % 11 12 13 -- mu, omega, g (attenomegatuation, albedo, anistropy)
        % 14 -- maximum # of absorption events
        % 15 -- Whether photon is on a surface
        
        onePhoton = photonValues(index, :);

        g1 = g1s(index);
        g2 = g2s(index);
        g3 = g3s(index);
        
        % Things that are accessed more than once *and* are unchanging
        killEnergy = onePhoton(10);
        mu = onePhoton(11);
        omega = onePhoton(12);
        g = onePhoton(13);

        onSurface = onePhoton(15);
        % 0 -- Photon is NOT on a surface
        % 1 -- Photon *is* residing on a surface
    
        absorptionData_ = zeros(onePhoton(14), 3);
        numAEvents_ = 0;
        exitInfo = zeros(1, 4);

        while 1 % Keep looping until specific conditions are met
                % Energy < kill energy threshold, or photon exits domain
            %numWhileLoopIter = numWhileLoopIter + 1;
            %{
            %% CONFIRM THAT PHOTON IS IN RIGHT LOCATION
            % Ray tracing region: Photon is "above" a surface
            % Absorptive region: Photon is below a surface
            % Find surface that has "overlapping" x-coordinates
            currentSurfaceIndex = find(xf >= onePhoton(1), 1, "first");
            if onePhoton(4) == 1 % Photon is ~ supposedly ~ in ray-tracing region
                if onePhoton(3) - (mF(currentSurfaceIndex)*onePhoton(1) + nF(currentSurfaceIndex)) >= -1e-10 
                else
                    disp("PHOTON IS REGION 1 BUT SHOULD BE REGION 2")
                end
            else % Photon is ~ supposedly ~ in absorptive region
                if onePhoton(3) - (mF(currentSurfaceIndex)*onePhoton(1) + nF(currentSurfaceIndex)) <= 1e-10 
                else
                    disp("PHOTON IS REGION 2 BUT SHOULD BE REGION 1")
                end
            end
            if onePhoton(1) < 0 || onePhoton(1) > xRightBoundary
                disp("PHOTON IS NOT WITHIN X_BOUNDS")
            end
            if onePhoton(3) < 0 || onePhoton(1) > zExitHeight
                disp("PHOTON IS NOT WITHIN Z_BOUNDS")
            end
            %}
                        % Particle's trajectory vector
            x1 = onePhoton(1); % initial X coordinate
            z1 = onePhoton(3); % initial Z coordinate
            x2 = x1 + onePhoton(5); % final X coordinate
            z2 = z1 + onePhoton(7); % final Z coordinate
            m = (z1-z2)/(x1-x2); % slopes

            photon_slope = onePhoton(7)/onePhoton(5);
            if onePhoton(4) == 1 % Ray-tracing region
                if onePhoton(7) > 0 && (photon_slope > maxSlope || photon_slope < minSlope) % Guaranteed to not hit a surface
                    intersectDistance = Inf;
                    %numNoIntersectShortcuts = numNoIntersectShortcuts + 1;
                else % Run intersection algorithm
                    %numIntersectAlg = numIntersectAlg + 1;
                    % Get all the x-coordinates of intersections between the
                    % current trajectory and snow surfaces
                    [xSc,zSc] = intersectSurfacePMC(x1,z1,x2,z2,xi,zi,xf,zf,zM,uC,x34,z34,determinant34);
                    
                    % Get the distance to the nearest (and correct) intersection
                    distance_xz = correctIntersectionPMC(xSc,zSc,x1,z1,onSurface,onePhoton(5),onePhoton(7),m);
                                
                    mag_xz = (onePhoton(5)^2 + onePhoton(7)^2)^0.5;
                    distance_y = distance_xz/mag_xz*onePhoton(6);
                    intersectDistance = (distance_xz^2 + distance_y^2)^0.5;
                    %if intersectDistance == Inf 
                    %    numIntersectInf = numIntersectInf + 1;
                    %end
                end

                if intersectDistance < Inf % Photon does strike a snow surface
                    %onePhoton(1:3) = onePhoton(1:3) + intersectDistance*onePhoton(5:7); % Update location
                    onePhoton(1) = onePhoton(1) + intersectDistance*onePhoton(5);
                    onePhoton(2) = onePhoton(2) + intersectDistance*onePhoton(6);
                    onePhoton(3) = onePhoton(3) + intersectDistance*onePhoton(7);
                    onePhoton(4) = 2; % Move photon to absorptive region
                    onSurface = 1; % Photon is on surface
                    continue
                else % Photon does NOT strike a snow surface

                    % Line properties
                    m = (z1-z2)/(x1-x2); % slopes
                    n = z1 - m.*x1; % z-intersects

                    xExitHeight = (zExit - n)/m; % Get exitHeight
                    if onePhoton(7) > 0 && xExitHeight > 0 && xExitHeight < maxX  % Photon is exiting!
                        
                        % Get zenith and azimuth angle (IN RADIANS) of exitance
                        zenith = acos(onePhoton(7)); % zenith angle wrt local basis vectors
                        azimuth = atan2(onePhoton(6),onePhoton(5)); % azimuth angle wrt local basis vectors
                        if azimuth < 0 % Ensure azimuth is in between 0 and 2pi
                            azimuth = azimuth + pie2; %% Precomputer 2*pi before parfor loop
                        end
                        % wavelength, energy, zenimuth, azimuth
                        exitInfo(1) = onePhoton(8);
                        exitInfo(2) = onePhoton(9);
                        exitInfo(3) = zenith;
                        exitInfo(4) = azimuth;
                        %onePhoton(9) = 0; % Drain energy
                        break
                    end

                    % Photon moves left...
                    if onePhoton(5) < 0
                        distToLeft = -onePhoton(1) / onePhoton(5);
                        onePhoton(2) = onePhoton(2) + distToLeft*onePhoton(6);
                        onePhoton(3) = onePhoton(3) + distToLeft*onePhoton(7);
                        onePhoton(1) = maxX; % Tunneled to right side

                    elseif onePhoton(5) > 0 % Photon moves right...
                        distToRight = (maxX - onePhoton(1)) / onePhoton(5);
                        onePhoton(2) = onePhoton(2) + distToRight*onePhoton(6);
                        onePhoton(3) = onePhoton(3) + distToRight*onePhoton(7);
                        onePhoton(1) = 0; % Tunneled to left side

                    end

                     % Particle's trajectory vector
                    x1 = onePhoton(1); % initial X coordinate
                    z1 = onePhoton(3); % initial Z coordinate
                    x2 = x1 + onePhoton(5); % final X coordinate
                    z2 = z1 + onePhoton(7); % final Z coordinate

                    %numTunneling = numTunneling + 1;

                    % After tunneling, do intersection alg again
                    % Get all the x-coordinates of intersections between the
                    % current trajectory and snow surfaces
                    [xSc,zSc] = intersectSurfacePMC(x1,z1,x2,z2,xi,zi,xf,zf,zM,uC,x34,z34,determinant34);
                    if length(xSc) == 1 && xSc(1) == 0  && zSc(1) == 0 % No intersections even after tunneling -- photon exited from ceiling
                         % Get zenith and azimuth angle (IN RADIANS) of exitance
                        zenith = acos(onePhoton(7)); % zenith angle wrt local basis vectors
                        azimuth = atan2(onePhoton(6),onePhoton(5)); % azimuth angle wrt local basis vectors
                        if azimuth < 0 % Ensure azimuth is in between 0 and 2pi
                            azimuth = azimuth + pie2; %% Precomputer 2*pi before parfor loop
                        end
                        % wavelength, energy, zenimuth, azimuth
                        exitInfo(1) = onePhoton(8);
                        exitInfo(2) = onePhoton(9);
                        exitInfo(3) = zenith;
                        exitInfo(4) = azimuth;
                        %onePhoton(9) = 0; % Drain energy
                        break
                    else  % After tunneling, photon strikes snow surface
                        % Correct Intersection is one with smallest x_distance away from current location
                        [x_dist, minDistIndex] = min(abs(xSc-onePhoton(1)));
                        z_dist = zSc(minDistIndex) - onePhoton(3);  % Positive or negative won't matter
                        distance_xz = (x_dist^2 + z_dist^2)^0.5;
                        
                        mag_xz = (onePhoton(5)^2 + onePhoton(7)^2)^0.5;
                        distance_y = distance_xz/mag_xz*onePhoton(6);
                        intersectDistance = (distance_xz^2 + distance_y^2)^0.5;

                        onePhoton(1) = onePhoton(1) + intersectDistance*onePhoton(5);
                        onePhoton(2) = onePhoton(2) + intersectDistance*onePhoton(6);
                        onePhoton(3) = onePhoton(3) + intersectDistance*onePhoton(7);
                        onePhoton(4) = 2; % Move photon to absorptive region
                        onSurface = 1;
                        continue
                    end
                end
            else % Photon is in absorptive region
                % Generate an advancement distance, dependent on attenuation coefficient
                advanceDistance = -log(rand)/mu;

                % Get intersection distance
                if onePhoton(7) < 0 && (onePhoton(3) < minZ  || photon_slope < minSlope || photon_slope > maxSlope)
                    intersectDistance = Inf;
                    %numNoIntersectShortcuts = numNoIntersectShortcuts + 1;
                else % Run intersection algorithm
                    %numIntersectAlg = numIntersectAlg + 1;
                    % Get all the x-coordinates of intersections between the
                    % current trajectory and snow surfaces
                    [xSc,zSc] = intersectSurfacePMC(x1,z1,x2,z2,xi,zi,xf,zf,zM,uC,x34,z34,determinant34);
                    
                    % Get the distance to the nearest (and correct) intersection
                    distance_xz = correctIntersectionPMC(xSc,zSc,x1,z1,onSurface,onePhoton(5),onePhoton(7),m);
                                
                    mag_xz = (onePhoton(5)^2 + onePhoton(7)^2)^0.5;
                    distance_y = distance_xz/mag_xz*onePhoton(6);
                    intersectDistance = (distance_xz^2 + distance_y^2)^0.5;
                    %if intersectDistance == Inf 
                    %    numIntersectInf = numIntersectInf + 1;
                    %end
                end
                
                if advanceDistance >= intersectDistance % Intersection will occur first
                    onePhoton(1) = onePhoton(1) + intersectDistance*onePhoton(5);
                    onePhoton(2) = onePhoton(2) + intersectDistance*onePhoton(6);
                    onePhoton(3) = onePhoton(3) + intersectDistance*onePhoton(7); % Update location
                    onePhoton(4) = 1; % Move photon to ray-tracing region
                    onSurface = 1; % Photon is on surface
                    continue
                else % Advancement and absorptive event occur, unless photon strikes sides of domain (check first)
                    
                    % Photon moves left...
                    if onePhoton(5) < 0
                        distToLeft = -onePhoton(1) / onePhoton(5);
                        if distToLeft < advanceDistance % photon will strike left side of domain
                            onePhoton(1) = maxX; % Tunneled to right side
                            onePhoton(2) = onePhoton(2) + distToLeft*onePhoton(6);
                            onePhoton(3) = onePhoton(3) + distToLeft*onePhoton(7);
                            continue
                        end
                    elseif onePhoton(5) > 0 % Photon moves right...
                        distToRight = (maxX - onePhoton(1)) / onePhoton(5);
                        if distToRight < advanceDistance % photon will strike right side of domain
                            onePhoton(2) = onePhoton(2) + distToRight*onePhoton(6);
                            onePhoton(3) = onePhoton(3) + distToRight*onePhoton(7);
                            onePhoton(1) = 0; % Tunneled to left side
                            continue
                        end
                    end
                
                    %onePhoton(1:3) = onePhoton(1:3) + advanceDistance*onePhoton(5:7); % Update location
                    onePhoton(1) = onePhoton(1) + advanceDistance*onePhoton(5);
                    onePhoton(2) = onePhoton(2) + advanceDistance*onePhoton(6);
                    onePhoton(3) = onePhoton(3) + advanceDistance*onePhoton(7);
                    
                    onSurface = 0; % Photon is not on surface
                    if onePhoton(3) < 0 % If photon leaves domain (exits bottom)
                        %% Photon has left domain from the bottom
                        energyThroughBottom = energyThroughBottom + onePhoton(9);
                        %onePhoton(9) = 0; % Drain energy
                        break
                    end
                    
                    % Get energy absorbed at location (1 - albedo)
                    energyAbsorbed = onePhoton(9)*(1-omega);
                    % Record absorption event
                    numAEvents_ = numAEvents_ + 1;
                    absorptionData_(numAEvents_, 1) = onePhoton(1);
                    absorptionData_(numAEvents_, 2) = onePhoton(3);
                    absorptionData_(numAEvents_, 3) = energyAbsorbed;
                    
                    onePhoton(9) = onePhoton(9)*omega; % Reduce energy
                    if onePhoton(9) < killEnergy % If energy is under kill energy threshold
                        numAEvents_ = numAEvents_ + 1;
                        absorptionData_(numAEvents_, 1) = onePhoton(1);
                        absorptionData_(numAEvents_, 2) = onePhoton(3);
                        absorptionData_(numAEvents_, 3) = onePhoton(9);
                        %onePhoton(9) = 0;
                        break
                    end
                    
                    % Generate a deflection zenith and azimuth angle
                    zen = acos(g1-((g2/(1+g-2*g*rand))^2)*g3); % [rad] zenith
                    azi = rand*pie2; % [rad] azimuth
                    
                    % Compute new trajectory unit vector
                    u0 = onePhoton(5); v0 = onePhoton(6); w0 = onePhoton(7); % old unit vectors
                    sinzen = sin(zen);
                    if abs(w0) > 0.9999
                        onePhoton(5) = sinzen*cos(azi);
                        onePhoton(6) = sinzen*sin(azi);
                        onePhoton(7) = w0*cos(zen)/abs(w0);
                    else
                        sinazi = sin(azi);  % Precomputing reused terms
                        cosazi = cos(azi);
                        w0cosazi = w0*cosazi;
                        coszen = cos(zen);
                        sqrt_term = sqrt(1-w0^2);
                        onePhoton(5) = sinzen/sqrt_term*(u0*w0cosazi-v0*sinazi) + u0*coszen; % new x comp
                        onePhoton(6) = sinzen/sqrt_term*(v0*w0cosazi+u0*sinazi) + v0*coszen; % new y comp
                        onePhoton(7) = -sinzen*cosazi*sqrt_term +w0*coszen; % new z comp
                    end
                    
                end % End of advance-reduce-deflect algorithm
            end % End of if-else for regions
        end % End of while loop (energy > kill-energy)

        aEvents{index} = absorptionData_(1:numAEvents_, :);
        numAEvents(index) = numAEvents_;
        eDataAll(index, :) = exitInfo;
    end % End of parfor loop (looping through each and every photon)

    if thermalStart == 1  % Solar photons weren't simulated
        aEventsSolar = zeros(0, 3);
        eData.solar = zeros(0, 4);
    else
        numAbsorptionEventsSolar = sum(numAEvents(1:thermalStart-1));
        aEventsSolar = zeros(numAbsorptionEventsSolar, 3);
        currIndex = 1;
        for index = 1:thermalStart-1
            newIndex = currIndex + numAEvents(index);
            aEventsSolar(currIndex:newIndex-1, :) = aEvents{index};
            currIndex = newIndex;
        end
        eData.solar = eDataAll(eDataAll(1:thermalStart-1, 2) ~= 0, :);
    end
    aData.solar = aEventsSolar;
    
    thermalAbsorbed = 0;
    numFacets = length(sData.xi);
    aEventsThermal = cell(numFacets, 1);
    if numPhotons - thermalStart + 1 == 0 % Thermal photons weren't sim-ed
        eData.thermal = zeros(0, 4);
    else
        aEventsThermal = cell(numFacets, 1);
        photonsPerFacet = (numPhotons - thermalStart + 1)/numFacets;

        for i = 1:numFacets
            indices = thermalStart + photonsPerFacet*(i-1) + (0:photonsPerFacet-1);
            aEvents_ = zeros(sum(numAEvents(indices)), 3);
            currIndex = 1; % Index of next empty row in aEvents_
            for j = indices(1):indices(end)  % Every photon simulated from this facet
                newIndex = currIndex + numAEvents(j);
                oneFacetEvents = aEvents{j};
                aEvents_(currIndex:newIndex-1, :) = oneFacetEvents;
                thermalAbsorbed = thermalAbsorbed + sum(oneFacetEvents(:, 3));
                currIndex = newIndex;
            end
            aEventsThermal{i} = aEvents_;
        end
        eDataThermal = eDataAll(thermalStart:end, :);
        eData.thermal = eDataThermal(eDataThermal(:, 2) ~= 0, :);
    end
    aData.thermal = aEventsThermal;

    % "Summary" energy numbers -- Used for energy conservation test
    energyData.absorbed = sum(aData.solar(:, 3)) + thermalAbsorbed;
    energyData.reflected = energyThroughBottom;
    energyData.exited = sum([eData.solar(:, 2); eData.thermal(:, 2)]);

    %toc;
end