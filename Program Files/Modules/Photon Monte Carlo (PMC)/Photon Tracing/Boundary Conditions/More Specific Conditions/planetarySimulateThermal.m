% Identical to planetarySimulate, but assumes that *only* thermal photons
% are simulated (no solar photons)

% Allows for less code after the parfor loop

function [aData, eData, energyData] = planetarySimulateThermal(sData, photons) %#codegen
    
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
    [xi,xf,zi,zf, ~, ~,zM, uC] = commonVariables(sData);
    [x34, z34, determinant34] = cachedInfo(xi, xf, zi, zf);
    maxX = max([xi xf(end)]);
    minZ = min([zi zf(end)]);
    zExit = max([zi zf(end)])*1.00;
    
    % Extract and group variables for efficiency
    gs = photonValues(:, 13); % scattering asymmetry parameter
    g1s = (1+gs.^2)*0.5./gs;
    g2s = 1-gs.^2;
    g3s = 0.5./gs;

    pie2 = 6.28318530718;

    parfor index = 1:numPhotons
        % 1 2 3 -- x, y, z locations
        % 4 -- region (1 is ray-tracing, 2 is absorptive)
        % 5 6 7 -- x, y, z trajectory
        % 8 -- wavelength
        % 9 10 -- energy and kill energy
        % 11 12 13 -- mu, omega, g (attenutation, albedo, anistropy)
        % 14 -- maximum # of absorption events
        % 15 -- Whether photon is on a surface (1 is on surface, 0 isn't)
        
        onePhoton = photonValues(index, :);

        % Retrieve all values that remain constant throughout tracing *and*
        % are accessed more than once
        killEnergy = onePhoton(10);
        mu = onePhoton(11);
        omega = onePhoton(12);
        g = onePhoton(13); 
        g1 = g1s(index); g2 = g2s(index); g3 = g3s(index);
        onSurface = onePhoton(15);
    
        absorptionData_ = zeros(onePhoton(14), 3);
        numAEvents_ = 0;
        exitInfo = zeros(1, 4);

        while 1 % Keep looping until specific conditions are met
                % Energy < kill energy threshold, or photon exits domain

            % Particle's trajectory vector
            x1 = onePhoton(1); % initial X coordinate
            z1 = onePhoton(3); % initial Z coordinate
            x2 = x1 + onePhoton(5); % final X coordinate
            z2 = z1 + onePhoton(7); % final Z coordinate
            m = (z1-z2)/(x1-x2); % slopes

            if onePhoton(4) == 1 % Ray-tracing region

                % Get all the x-coordinates of intersections between the
                % current trajectory and snow surfaces
                [xSc,zSc] = intersectSurfacePMC(x1,z1,x2,z2,xi,zi,xf,zf,zM,uC,x34,z34,determinant34);
                
                % Get projected distance to the nearest/next intersection
                distance_xz = correctIntersectionPMC(xSc,zSc,x1,z1,onSurface,onePhoton(5),onePhoton(7),m);
                
                % Get 3D distance to next intersection
                mag_xz = (onePhoton(5)^2 + onePhoton(7)^2)^0.5;
                distance_y = distance_xz/mag_xz*onePhoton(6);
                intersectDist = (distance_xz^2 + distance_y^2)^0.5;

                if intersectDist < Inf % Photon does strike a snow surface

                    % Move photon to intersection coordinate
                    onePhoton(1) = onePhoton(1) + intersectDist*onePhoton(5);
                    onePhoton(2) = onePhoton(2) + intersectDist*onePhoton(6);
                    onePhoton(3) = onePhoton(3) + intersectDist*onePhoton(7);
                    onePhoton(4) = 2; % Move photon to absorptive region
                    onSurface = 1; % Photon is on surface
                    continue

                else % Photon does NOT strike a snow surface

                    % Line properties
                    m = (z1-z2)/(x1-x2); % slopes
                    n = z1 - m.*x1; % z-intersects

                    xExitHeight = (zExit - n)/m;
                    if onePhoton(7) > 0 && xExitHeight > 0 && xExitHeight < maxX  % Photon is exiting
                        
                        % Get zenith and azimuth angle (radians) of exitance
                        zenith = acos(onePhoton(7)); % zenith angle wrt local basis vectors
                        azimuth = atan2(onePhoton(6),onePhoton(5)); % azimuth angle wrt local basis vectors
                        if azimuth < 0 % Ensure azimuth is in between 0 and 2pi
                            azimuth = azimuth + pie2;
                        end
                        % wavelength, energy, zenimuth, azimuth
                        exitInfo(1) = onePhoton(8);
                        exitInfo(2) = onePhoton(9);
                        exitInfo(3) = zenith;
                        exitInfo(4) = azimuth;
                        break
                    end

                    % Photon moves left...
                    if onePhoton(5) < 0
                        tunnelDistance = -onePhoton(1) / onePhoton(5);
                        onePhoton(2) = onePhoton(2) + tunnelDistance*onePhoton(6);
                        onePhoton(3) = onePhoton(3) + tunnelDistance*onePhoton(7);
                        onePhoton(1) = maxX; % Tunneled to right side

                    elseif onePhoton(5) > 0 % Photon moves right...
                        tunnelDistance = (maxX - onePhoton(1)) / onePhoton(5);
                        onePhoton(2) = onePhoton(2) + tunnelDistance*onePhoton(6);
                        onePhoton(3) = onePhoton(3) + tunnelDistance*onePhoton(7);
                        onePhoton(1) = 0; % Tunneled to left side

                    end

                    % Particle's trajectory vector
                    x1 = onePhoton(1); % initial X coordinate
                    z1 = onePhoton(3); % initial Z coordinate
                    x2 = x1 + onePhoton(5); % final X coordinate
                    z2 = z1 + onePhoton(7); % final Z coordinate

                    % Get all the x-coordinates of intersections between the
                    % current trajectory and snow surfaces
                    [xSc,zSc] = intersectSurfacePMC(x1,z1,x2,z2,xi,zi,xf,zf,zM,uC,x34,z34,determinant34);

                    if length(xSc) == 1 && xSc(1) == 0  && zSc(1) == 0 % No intersections even after tunneling -- photon exited from ceiling
                         % Get zenith and azimuth angle (radians) of exitance
                        zenith = acos(onePhoton(7)); % zenith angle wrt local basis vectors
                        azimuth = atan2(onePhoton(6),onePhoton(5)); % azimuth angle wrt local basis vectors
                        if azimuth < 0 % Ensure azimuth is in between 0 and 2pi
                            azimuth = azimuth + pie2;
                        end
                        % wavelength, energy, zenimuth, azimuth
                        exitInfo(1) = onePhoton(8);
                        exitInfo(2) = onePhoton(9);
                        exitInfo(3) = zenith;
                        exitInfo(4) = azimuth;
                        break
                    else  % After tunneling, photon strikes snow surface

                        % Correct Intersection is one with smallest x_distance away from current location
                        [x_dist, minDistIndex] = min(abs(xSc-onePhoton(1)));
                        z_dist = zSc(minDistIndex) - onePhoton(3);
                        distance_xz = (x_dist^2 + z_dist^2)^0.5;
                        
                        mag_xz = (onePhoton(5)^2 + onePhoton(7)^2)^0.5;
                        distance_y = distance_xz/mag_xz*onePhoton(6);
                        intersectDist = (distance_xz^2 + distance_y^2)^0.5;

                        % Move photon to intersection coordinate
                        onePhoton(1) = onePhoton(1) + intersectDist*onePhoton(5);
                        onePhoton(2) = onePhoton(2) + intersectDist*onePhoton(6);
                        onePhoton(3) = onePhoton(3) + intersectDist*onePhoton(7);
                        onePhoton(4) = 2; % Move photon to absorptive region
                        onSurface = 1;
                        continue
                    end
                end
            else % Photon is in absorptive region
                % Generate an advancement distance
                advanceDistance = -log(rand)/mu;

                % Under these conditions, photon is guaranteed to not hit a
                % surface (no need to run expensive intersection algorithm)

                % 1st shortcut: Even a directly upward path won't hit the surface
                % 2nd shortcut: Photon is heading downwards and is already below minZ
                if onePhoton(3) + advanceDistance < minZ || (onePhoton(7) < 0 && onePhoton(3) < minZ)
                    intersectDist = Inf;
                else % Run intersection algorithm

                    % Get all the x-coordinates of intersections between the
                    % current trajectory and snow surfaces
                    [xSc,zSc] = intersectSurfacePMC(x1,z1,x2,z2,xi,zi,xf,zf,zM,uC,x34,z34,determinant34);
                    
                    % Get the projected distance to the nearest/next intersection
                    distance_xz = correctIntersectionPMC(xSc,zSc,x1,z1,onSurface,onePhoton(5),onePhoton(7),m);
                    
                    % Get 3D distance to next intersection
                    mag_xz = (onePhoton(5)^2 + onePhoton(7)^2)^0.5;
                    distance_y = distance_xz/mag_xz*onePhoton(6);
                    intersectDist = (distance_xz^2 + distance_y^2)^0.5;

                end
                
                if advanceDistance >= intersectDist % Intersection will occur first
                     % Update location to intersection
                    onePhoton(1) = onePhoton(1) + intersectDist*onePhoton(5);
                    onePhoton(2) = onePhoton(2) + intersectDist*onePhoton(6);
                    onePhoton(3) = onePhoton(3) + intersectDist*onePhoton(7);
                    onePhoton(4) = 1; % Move photon to ray-tracing region
                    onSurface = 1; % Photon is on surface
                    continue

                else % Advancement and absorptive event occur, unless photon strikes sides of domain
                    
                    tunneled = 0;

                    % Photon moves left...
                    if onePhoton(5) < 0
                        tunnelDistance = -onePhoton(1) / onePhoton(5);
                        if tunnelDistance < advanceDistance % photon will strike left side of domain
                            onePhoton(1) = maxX; % Tunneled to right side
                            onePhoton(2) = onePhoton(2) + tunnelDistance*onePhoton(6);
                            onePhoton(3) = onePhoton(3) + tunnelDistance*onePhoton(7);
                            advanceDistance = advanceDistance - tunnelDistance;
                            tunneled = 1;
                        end
                    elseif onePhoton(5) > 0 % Photon moves right...
                        tunnelDistance = (maxX - onePhoton(1)) / onePhoton(5);
                        if tunnelDistance < advanceDistance % photon will strike right side of domain
                            onePhoton(1) = 0; % Tunneled to left side
                            onePhoton(2) = onePhoton(2) + tunnelDistance*onePhoton(6);
                            onePhoton(3) = onePhoton(3) + tunnelDistance*onePhoton(7);
                            advanceDistance = advanceDistance - tunnelDistance;
                            tunneled = 1;
                        end
                    end

                    if tunneled

                        % Particle's trajectory vector
                        x1 = onePhoton(1); % initial X coordinate
                        z1 = onePhoton(3); % initial Z coordinate
                        x2 = x1 + onePhoton(5); % final X coordinate
                        z2 = z1 + onePhoton(7); % final Z coordinate
    
                        % Under these conditions, photon is guaranteed to not hit a
                        % surface (no need to run expensive intersection algorithm)
        
                        % 1st shortcut: Even a directly upward path won't hit the surface
                        % 2nd shortcut: Photon is heading downwards and is already below minZ
                        if onePhoton(3) + advanceDistance < minZ || (onePhoton(7) < 0 && onePhoton(3) < minZ)
                            intersectDist = Inf;
                        else
                            % Rerun intersection algorithm with tunneled location
                            [xSc,zSc] = intersectSurfacePMC(x1,z1,x2,z2,xi,zi,xf,zf,zM,uC,x34,z34,determinant34);
    
                            if length(xSc) > 1 || (xSc(1) ~= 0 || zSc(1) ~= 0) % Intersection after tunneling
                                % Correct Intersection is one with smallest x_distance away from current location
                                [x_dist, minDistIndex] = min(abs(xSc-onePhoton(1)));
                                z_dist = zSc(minDistIndex) - onePhoton(3);
                                distance_xz = (x_dist^2 + z_dist^2)^0.5;
                                
                                mag_xz = (onePhoton(5)^2 + onePhoton(7)^2)^0.5;
                                distance_y = distance_xz/mag_xz*onePhoton(6);
                                intersectDist = (distance_xz^2 + distance_y^2)^0.5;
                            else % No intersection post-tunneling
                                intersectDist = Inf;
                            end
                        end

                        if advanceDistance >= intersectDist % Intersection will occur after tunneling
                             % Update location to intersection
                            onePhoton(1) = onePhoton(1) + intersectDist*onePhoton(5);
                            onePhoton(2) = onePhoton(2) + intersectDist*onePhoton(6);
                            onePhoton(3) = onePhoton(3) + intersectDist*onePhoton(7);
                            onePhoton(4) = 1; % Move photon to ray-tracing region
                            onSurface = 1; % Photon is on surface
                            continue
                        end % Otherwise, absorption event will occur
                    end
                
                    % Update location
                    onePhoton(1) = onePhoton(1) + advanceDistance*onePhoton(5);
                    onePhoton(2) = onePhoton(2) + advanceDistance*onePhoton(6);
                    onePhoton(3) = onePhoton(3) + advanceDistance*onePhoton(7);
                    
                    onSurface = 0; % Photon is not on surface
                    if onePhoton(3) < 0 % Photon exits domain through floor
                        energyThroughBottom = energyThroughBottom + onePhoton(9);
                        break
                    end
                    
                    % Get energy absorbed at location
                    energyAbsorbed = onePhoton(9)*(1-omega); % [Joules]

                    % Record absorption event
                    % x and z coordinates, amount of energy absorbed
                    numAEvents_ = numAEvents_ + 1;
                    absorptionData_(numAEvents_, 1) = onePhoton(1);
                    absorptionData_(numAEvents_, 2) = onePhoton(3);
                    absorptionData_(numAEvents_, 3) = energyAbsorbed;
                    
                    onePhoton(9) = onePhoton(9)*omega; % Reduce energy
                    if onePhoton(9) < killEnergy

                        % Drain photon at current location
                        numAEvents_ = numAEvents_ + 1;
                        absorptionData_(numAEvents_, 1) = onePhoton(1);
                        absorptionData_(numAEvents_, 2) = onePhoton(3);
                        absorptionData_(numAEvents_, 3) = onePhoton(9);

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

    aData.solar = zeros(0, 3);
    eData.solar = zeros(0, 4);
    
    thermalAbsorbed = 0;
    numFacets = length(sData.xi);
    aEventsThermal = cell(numFacets, 1);
    parfor index = 1:numFacets  % Cell array elements must be assigned
        aEventsThermal{index} = index;
    end
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
    aData.thermal = aEventsThermal;

    % "Summary" energy numbers -- Used for energy conservation test
    energyData.absorbed = sum(aData.solar(:, 3)) + thermalAbsorbed;
    energyData.reflected = energyThroughBottom;
    energyData.exited = sum([eData.solar(:, 2); eData.thermal(:, 2)]);

end