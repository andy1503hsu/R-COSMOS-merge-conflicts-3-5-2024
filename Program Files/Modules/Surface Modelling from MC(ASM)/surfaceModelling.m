function [xNew,yNew] = surfaceModelling(surface,dispFacet,type) %#codegen
%% DESCRIPTION:
%       Function adjusts the surface according to the displacement of the
%   facets in space. The inputs should be an array containing the X and Y
%   coordinates of the surface, a defined rectangular domain, and the
%   displacement of all the facets.
%
%       Input variables:
%           surface:
%               xS --> Y coordinates of the original surface [m]
%               yS --> Y coordinates of the original surface [m]
%           domain:
%               h --> height of the domain [m]
%               w --> width of the domain [m]
%
%           dispFacet --> displacement of all the facets [m]
%               if "dispFacet(facet)" > 0 --> facet moved outwards 
%               || "dispFacet(facet)" < 0 --> facet moved inwards
%               || "dispFacet(facet)" = 0 --> facet did not moved

                                %% Steps %%
%   1. Extract the surface, domain information, and calculate the area of
% the facets and their angle with respect to the horizontal. The area is
% measured from one surface node to another, taking in consideration a unit
% depth into the screen. Calculate the angle of each facet with respect to
% the horizontal and adjust it to be normal to the facet pointing outwards
% from the surface.
%
%   variable "angleFacet" --> angle w.r.t horizontal in degrees
%
%   Compare the Y components between surface nodes:
%       if yS(2) > yS(1) -->  angleFacet = 90 +  angleFacet
%       || yS(2) < yS(1) -->  angleFacet = 90 -  angleFacet
%       || yS(2) = yS(1) -->  angleFacet = 90
%
%   Note: Units are in degrees
% The displacement of each facet is in a direction 
% normal to the facet.

% -------------------------------------------------------------------------
%   2. Find the location and number of critical points in the original
% surface. Establish a dynamic tolerance based on the surface.
%       2.1 Find the location and number of flat maxima and minima.
%               Compare the Y components between 2 adjacent surface nodes:
%                   if yS(1) - yS(2) < tolerance then it is a flat surface
%
%               Note: The critical point defined by a flat surface will be
%               only yS(1).
%               
%               Note: The concept of a flat surface here is defined in this
%               function as a flat maximun or a flat minimun.
%
%       2.2 Find the location and number of maxima and minima.
%               Compare the Y components between 3 surface nodes:
%                   -yS(1) is the first surface node in the loop
%                   -yS(3) is the last surface node in the loop
%                   -yS(2) is a surface node adjacent to yS(1) and yS(3),
%                   therefore is in between them.
%
%                   if yS(2) > (yS(1) && yS(3)), then yS(2) is a maximun
%                   || yS(2) < (yS(1) && yS(3)), then yS(2) is a minimun
%
%               Note: The loop runs from 2 to the length(yS)-1, meaning 
%               that it does not check for these conditions on the 
%               endpoints. This is designed to maintain consistency through 
%               the model. The endpoints are adjusted in step 7.

% -------------------------------------------------------------------------
%   3. Adjust the position in space of each facet. Move each facet 
% according to its respective displacement and direction. Each facet is
% stored in a cell that contains the X and Y coordinates of the endpoints
% of the facet.

% -------------------------------------------------------------------------
%   4. Adjust the surface 
%       4.1 Each critical location.
%           1. Join the critical locations together in the same array.
%           2. Sort the critical locations in ascending order.
%           3. Adjust the surface at each critical location by identifying
%           whether it is a flat surface, a minimun, or a maximun. See
%           additional functions for this step that corresponds to each
%           scenario.
%
%   Note: A variable called "I" is created in this step and it is in charge
%   of adjusting the location of the critical point based on how many cells
%   are deleted in the adjusting process.
%
%       4.2 Begining and end of the domain
%           1. Check how many facets have completely crossed the START of 
%           the domain and delete them.
%           2. Find the intersection between the START of the domain and
%           the FIRST facet and adjust it to intersect the Y axis.
%
%       Note: Process is SIMILAR for the end of the domain
%           3. Check how many facets have completely crossed the END of the
%           domain and delete them.
%           4. Find the intersection between the END of the domain and
%           the LAST facet and adjust it to intersect the Y axis.

% -------------------------------------------------------------------------
%   5. Calculate the average of each facet. This accounts for the rotation 
% of the facets and their increase/decrease in size with time. The
% average is calculated by taking the average between the X and Y
% coordinates of the adjacent points that correspond to each facet.

% -------------------------------------------------------------------------
%   6. Check the final surface for imperfection such as crossing of
% different facets and adjust the surface accordingly. This step helps
% maintain the surface free from imperfections and things that could
% affect future calculations. Crossings of the facets are checked 

% -------------------------------------------------------------------------
%   7. For consistency throughout the model, the endpoints are adjusted
% according to the process described in step 6. To maintain a periodic
% boundary condition (in case that it is selected) the Y coordinate of
% the final node will be equal to the Y coordinate of the initial node.

% -------------------------------------------------------------------------
%   8. To maintain compatibility with any other model or function, the
% number of surface nodes deleted during these processes are restored in a
% random location along the new surface.
%       8.1 Extract the new surface with established dimensions

%% 1. Surface data and domain
    % Surface coordinates
    xS = surface.xS;  yS = surface.zS;

    % Domain definition
    h = max(yS); % height of the domain [m]
    w = xS(end); % width of the domain [m]    

    % Area of the facets
    areaFacet = zeros(1,length(xS)-1); % preallocate the area of the facets
    for ii = 1:length(xS)-1
        areaFacet(ii) = ((xS(ii)-xS(ii+1))^2 + (yS(ii)-yS(ii+1))^2)^0.5;
    end

    % Angles w.r.t horizontal
    angleFacet = acosd((xS(1,2:end)-xS(1,1:end-1))./areaFacet);
    %%% Check and adjust the angle to point outwards from the surface
    for ind = 1:length(xS)-1
        if yS(ind+1) > yS(ind)
            angleFacet(ind) = 90 + angleFacet(ind);            
        elseif yS(ind+1) < yS(ind)
            angleFacet(ind) = 90 - angleFacet(ind);            
        else
            angleFacet(ind) = 90 + 0;
        end
    end

%% 2. Location of critical points in the original surface
%% 2.1 Flat surfaces
    tol = 0; % dynamic tolerance
    % Preallocate temporary array to store the location of flat surfaces
    tempFlat = zeros(1,length(yS)-2);
    %%% Find and count the number of flat surfaces
    countFlat = 0; % initialize the flat surfaces count
    for node = 2:length(yS)-2
        if abs(yS(node) - yS(node+1)) < tol
            % Add 1 to the current number of flat surfaces
            countFlat = countFlat + 1;
            % Store them in a temporary array
            if yS(node) > yS(node-1) && yS(node+1) > yS(node+2)
                tempFlat(1,countFlat) = node;
            elseif yS(node) < yS(node-1) && yS(node+1) < yS(node+2)
                tempFlat(1,countFlat) = node;
            elseif yS(node) > yS(node-1) && yS(node+1) < yS(node+2)
                tempFlat(1,countFlat) = node;
            elseif yS(node) < yS(node-1) && yS(node+1) > yS(node+2)
                tempFlat(1,countFlat) = node;
            else 
                countFlat = countFlat - 1;
            end
        end
    end

    %%% Store the location and the number of flat surfaces
    if countFlat == 0 % if there are not flat surfaces
        flat = []; % return an empty array
    else % if there are flat surfaces
        flat = tempFlat(1:countFlat);
    end

%% 2.2 Local maxima and minima
    % Preallocate temporary arrays to store the location of maxima and
    % minima
    tempMax = zeros(1,length(yS)-2);
    tempMin = zeros(1,length(yS)-2);

    %%% Find and count the number of local maxima and minima
    countMax = 0; % initialize the maxima count
    countMin = 0; % initialize the minima count
    for loc = 2:length(yS)-1
        if yS(loc) > yS(loc-1) && yS(loc) > yS(loc+1)
            countMax = countMax + 1; % update the maxima count
            tempMax(1,countMax) = loc; % store location of the maxima
            
        elseif yS(loc) < yS(loc-1) && yS(loc) < yS(loc+1)
            countMin = countMin + 1; % update the minima count
            tempMin(1,countMin) = loc; % store location of the minima
        end
    end

    %%% Store the location and the number of maxima
    if countMax == 0 % if there are no maxima
        maxima = []; % return an empty array
    else
        maxima = tempMax(1:countMax); % store the location of the maxima
    end

    %%% Store the location and the number of minima
    if countMin == 0 % if there are no minima
        minima = []; % return an empty array
    else
        minima = tempMin(1:countMin); % store the location of the minima
    end

%% 3. Position of each facet in space
    % Preallocate the cell array containing the facets    %units here
    facets = ones(4,length(xS)-1);
    coder.varsize('facets')
    % Adjust the facets based on their displacement and direction and store
    % their X and Y coordinates in the "facets" cell array.
    for i = 1:length(xS)-1
        % X coordinates of the 2 nodes in the facet
        facets(1,i) = xS(i)   + dispFacet(1,i) * cosd(angleFacet(1,i));
        facets(3,i) = xS(i+1) + dispFacet(1,i) * cosd(angleFacet(1,i));

        % Y coordinates of the 2 nodes in the facet
        facets(2,i) = yS(i)   + dispFacet(1,i) * sind(angleFacet(1,i));
        facets(4,i) = yS(i+1) + dispFacet(1,i) * sind(angleFacet(1,i));
    end
% plot([facets(1,:) facets(3,:)],[facets(2,:) facets(4,:)],'g-');
%% 4. Adjust the surface
%% 4.1 Critical points
    I = 0; % Initialize the global index variable in charge of adjusting 
           % the index of the critical locations throught the entire code.

    % Concatenate the critical locations
    criticalLocations = [flat maxima minima];

    % Sort the critical locations in ascending order
    sortedLocations = sort(criticalLocations);

    % Run through each critical point
    for jj = 1:length(criticalLocations)   
        % Extract the current location
        currentLoc = sortedLocations(jj);

        %%% Check if it is a flat surface and adjust it
        for locFlat = 1:length(flat)
            if currentLoc == flat(locFlat)
                [facets,I] = adjustFlat(facets,currentLoc,dispFacet,I);
            end
        end

        %%% Check if it is a maximun and adjust it
        for locMax = 1:length(maxima)
            if currentLoc == maxima(locMax)
                [facets,I] = adjustMaximun(facets,currentLoc,dispFacet,I);
            end
        end

        %%% Check if it is a minimun and adjust it
        for locMin = 1:length(minima)
            if currentLoc == minima(locMin)
                [facets,I] = adjustMinimun(facets,currentLoc,dispFacet,I);
            end
        end
    end

%% 4.2 Begining and end of the domain
    %%% Count how many facets have crossed the start of the domain.
    countStart = 0; % initialize the count
    for facet = 1:length(facets)
        if facets(3,facet) < 0 % if the coordinate of node 2 is < 0
            countStart = countStart + 1; % add 1 to the count
        end
    end
    % Delete the facets 
    if countStart > 0
        facets(:,1:countStart) = [];
        I = I + countStart;
    end

    %%% Adjust the begining point to intersect the Y axis for consistency.
    x0 = facets(3,1);
    y0 = facets(4,1);
    x1 = facets(1,1);
    y1 = facets(2,1);

    % Calculate the intersection between the first facet and the Y axis.
    yLc = (y1*x0 - y0*x1)/(x0-x1);
    if xS(1) < yLc && yLc < h && x1 < 0 % if there are intersections
        facets(1,1) = 0; % adjust X coordinate
        facets(2,1) = yLc; % adjust Y coordinate
    end

           %%% Process is similar for the end of the domain %%%
    %%% Count how many facets have crossed the end of the domain
    countEnd = 0; % initialize the count
    for facet = 1:length(facets)
        if facets(1,facet) > w % if the coordinate of node 1 is > w
            countEnd = countEnd + 1; % add 1 to the count
        end
    end
    % Delete the facets
    if countEnd > 0
        facets(:,end:-1:end-countEnd+1) = []; % delete those facets
        I = I + countEnd;
    end

    %%% Adjust the end point to intersect the Y axis for consistency.
    x0 = facets(1,end);
    y0 = facets(2,end);
    x1 = facets(3,end);
    y1 = facets(4,end);

    % Calculate the intersection between the last facet and the Y axis.
    yRc = (y1*(x0-w) - y0*(x1-w))/(x0-x1);
    if xS(1) < yRc && yRc < h && x1 > w % if there are intersections
        facets(3,end) = w; % adjust X coordinate
        facets(4,end) = yRc; % adjust Y coordinate
    end   

%% 5. Middle Points
    % Calculate the average between the adjacent points of two facets to
    % account for rotation and size of the facets.
    midX = zeros(1,length(facets)-1); % preallocate the X coordinates
    midY = zeros(1,length(facets)-1); % preallocate the Y coordinates
    for point = 1:length(facets)-1
        midX(point) = (facets(3,point) + facets(1,point+1))/2;
        midY(point) = (facets(4,point) + facets(2,point+1))/2;
    end

%% 6. Real Surface Imperfections
%% 6.1 Crossings
    %%% Check the surface nodes for intersections between facets
    tempStore = zeros(length(midX)-3,4);
    counting = 0; % count the imperfections
    for firstFacet = 1:length(midX) - 3 % for non-adjacent facets
        % X and Y coordinates of the first facet
        xF = [midX(firstFacet) midX(firstFacet+1)];
        yF = [midY(firstFacet) midY(firstFacet+1)];

        for secondaryFacet = firstFacet+3:length(midX)
            % X and Y coordinate of the third and successive facets
            xSF = [midX(secondaryFacet-1) midX(secondaryFacet)];
            ySF = [midY(secondaryFacet-1) midY(secondaryFacet)];

            % Check for an intersection between x,y and xS,yS
            [xSc,ySc] = intersectSurface2(xF,yF,xSF,ySF);

            % If there is an intersection between the first and third facet
            if ~(isempty(xSc)) && xSc(1,1) > xF(1) && xSc(1,1) < xF(2)
                counting = counting + 1; % add 1 to the count

                % Store the locations and the intersect
                tempStore(counting,1) = firstFacet;
                tempStore(counting,2) = secondaryFacet;
                tempStore(counting,3) = xSc(1,1);
                tempStore(counting,4) = ySc(1,1);

                % If the first facet has more than 1 intersection with a
                % secondary facet, keep the last intersection as true
                if counting > 1
                    if tempStore(counting,1) == tempStore(counting-1,1)
                        tempStore(counting-1,:) = [];
                        counting = counting - 1;
                    end
                end
            end
        end
    end

    %%% Check how many cells need to be deleted in order to remove the
    %%% imperfections on the surface
    tI = 0; % temporary index: this variable has the same function as "I"
    for node = 1:counting
        % Calculate how many facets are in between the intercepting facets
        diff = tempStore(node,2) - tempStore(node,1);

        % Adjust the new node to be on the intersection
        midX(tempStore(node,1)+1-tI) = tempStore(node,3);
        midY(tempStore(node,1)+1-tI) = tempStore(node,4);

        % If the difference is 3 there is only one facet
        if diff == 3
            % Delete next surface node
            midX(tempStore(node,1)+2-tI) = []; % delete X coord
            midY(tempStore(node,1)+2-tI) = []; % delete Y coord
            tI = tI + 1;

        % If the difference is bigger than 3 there is more than one facet
        elseif diff > 3
            % Delete all surface nodes in between
            midX(tempStore(node,1)+2-tI:tempStore(node,1)+diff-1-tI) = [];
            midY(tempStore(node,1)+2-tI:tempStore(node,1)+diff-1-tI) = [];
            tI = tI + diff - 2;
        end
    end
    I = I + tI;

%% 6.2 Non-Functions
    upH = zeros(1,length(midX)-2); % preallocate the uphill movement
    downH = zeros(1,length(midX)-2); % preallocate the downhill movement
    countUp = 0; % initialize the uphill non function count
    countDown = 0; % initialize the downhill non function count
    % Check the surface for non-function conditions 
    for i = 2:length(midX)-1
        % Slope of the current facet
        m = (midY(i-1)-midY(i))/(midX(i-1)-midX(i));
        % If there are 2 values of y for each x coordinate it means that a
        % non-function condition is present.
        if midX(i) > midX(i+1)
            % Check the slope sign to identify the direction of travel
            if m > 0
                countUp = countUp + 1; % update the uphill count
                upH(countUp) = i; % store the current location
            elseif  m < 0
                countDown = countDown + 1; % update the downhill count
                downH(countDown) = i; % store the current location
            end
        end
    end

    % Store the current locations of nodes that represent a non-function
    % condition, organize them, and sort them in an ascending order.
    if countUp == 0
        upHillLoc = []; % return empty array if there are no conditions
    else
        upHillLoc = upH(1:countUp); % store the nodes' location
    end

    if countDown == 0
        downHillLoc = []; % return empty array if there are no conditions
    else 
        downHillLoc = downH(1:countDown); % store the nodes' location
    end

    % Store the location of the nodes and sort them in ascending order 
    nodeLoc = sort([downHillLoc upHillLoc]);

    tIndex = 0; % initialize temporary index
    for a1 = 1:length(nodeLoc)
        % Retrieve the location of the node
        currentNode = nodeLoc(a1);

        % Check if the current node corresponds to the uphill case
        for u = 1:length(upHillLoc)
            if currentNode == upHillLoc(u)
                % Remove the current surface node
                midX(upHillLoc(u)-tIndex+1) = [];
                midY(upHillLoc(u)-tIndex+1) = [];
                tIndex = tIndex + 1; % update temporary index
            end
        end

        % Check if the current node corresponds to the downhill case
        for d = 1:length(downHillLoc)
            if currentNode == downHillLoc(d)
                % Remove the next surface node
                midX(downHillLoc(d)-tIndex) = [];
                midY(downHillLoc(d)-tIndex) = [];
                tIndex = tIndex + 1; % update temporary index
            end
        end
    end
    I = I + tIndex; % update global index

%% 7. Adjust endpoints
    % Check the type of domain used and adjust accordingly
    if type == "Planetary"
        % Initial surface node X and Y coordinates
        startX = 0;
        startY = facets(2,1);

        % Final Surface node X and Y coordinates
        endX = w;
        endY = facets(2,1);

    elseif type == "ArkChamber" || type == "ArkChamber2"
        % Initial surface node X and Y coordinates
        startX = facets(1,1);
        startY = 0;

        % Final Surface node X and Y coordinates
        endX = facets(3,end);
        endY = 0;
    else
        error('Select a correct Right boundary condition')
    end

    %%% Store the X and Y coordinates
    tX = [startX midX endX];
    tY = [startY midY endY];

%% 8. Restore compatibility
    % Calculate the number of points missing based on "I" and assign then 
    % a random position in the surface
    extraPoints = startX + (endX-startX)*rand(1,I);

    % Organize the extra points in ascending order
    sortedExtra = sort(extraPoints);

    % Preallocate the final X coordinates
    xN = zeros(1,length(tX)+I);
    xN(1:length(tX)) = tX; % existing values of X

    % Preallocate the final Y coordinates
    yN = zeros(1,length(tY)+I);
    yN(1:length(tY)) = tY; % existing values of Y

    %%% Add the extra points to the final coordinates array to restore
    %%% compatibility
    for eP = 1:length(extraPoints)
        % Check each interval to find the location of the extra point
        for iT = 1:length(xN)-1
            % If the extra point is in between the interval given by the X 
            % coordinates, update the final coordinates array
            if sortedExtra(eP) > xN(iT) && sortedExtra(eP) < xN(iT+1)

                % Calculate the average of the real points in the surface 
                avgX = (xN(iT) + xN(iT+1))/2;
                avgY = (yN(iT) + yN(iT+1))/2;

                % Make space for the extra point
                xN((iT+2):(length(tX)+eP)) = xN((iT+1):(length(tX)+eP-1));
                yN((iT+2):(length(tX)+eP)) = yN((iT+1):(length(tX)+eP-1));

                % Add the extra point to the real surface
                xN(iT+1) = avgX;
                yN(iT+1) = avgY;
                break
            end
        end
    end

%% 8.1 New surface
    %%% Extract the new surface with established dimensions
    xNew = zeros(1,length(xS)); % preallocate final surface
    yNew = zeros(1,length(xS)); % preallocate final surface
    for node = 1:length(xS)
        xNew(1,node) = xN(1,node); % equal to the xN array
        yNew(1,node) = yN(1,node); % equal to the yN array
    end
end