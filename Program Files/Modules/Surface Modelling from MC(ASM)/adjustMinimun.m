function [facets,I] = adjustMinimun(facets,minimun,dispFacet,I) %#codegen
%% DESCRIPTION
%       Function adjusts the surface at a minimun location. If the minimun
%   is expanding then it is adjusted, if it receeds, it is left intact.
%
%       Input variables:
%           facets --> array containing all the facets
%           minimun --> index for the location of the minimun
%           dispFacet --> displacement of the facet in space
%           I --> global index

                                %% Steps %%
%   1. Check if the surface is contracting at the location of the minimun
%   by taking the average of the displacement of the facets adjacent to the
%   minimun.

% -------------------------------------------------------------------------
%   2. Find the point at which the surface converges and calculate the
%   initial and final points of the adjacent facets.
%       2.1 Remove the cell that contains the facets that are no longer in
%       use to adjust the surface.

% -------------------------------------------------------------------------
%   3. Check the new peak for an intersection and adjust it accordingly.

%% Adjust the minimun
    %%% If the average of facets adjacent to the minimun is negative,
    %%% meaning that they minimun is receding:
    if (dispFacet(minimun-I-1)+dispFacet(minimun-I))/2 > 0
        % Find which pair of facets has contracted enough to create a new peak
        loc = 0; % initialize location index
        for i = 1:length(facets)/2
            loc = i; % store index in location variable
            % Final X coordinate of the facet on the left
            p2X = facets(3,minimun-I-loc);
            % Initial X coordinate of the facet on the right
            p3X = facets(1,minimun-I+loc-1);

            % Check if the first facet has been reached
            if p2X > p3X && p2X == facets(3,1)
                facets(:,minimun-I-loc) = [];
                I = I + 1;
                return
            end
            % Check if the last facet has been reached
            if p2X > p3X && p3X == facets(1,end)
                facets(:,minimun-I+loc-1) = [];
                I = I + 1;
                return
            end
            % If they are in the correct position, break the loop
            if p2X < p3X % if there is no intersection, break
                break
            end
        end

        % Remove the cells that contain the facets over
        tempIndex = 0;
        for j = 1:(loc-2)*2
            tempIndex = j;
            facets(:,minimun-I-loc+2) = [];
        end
        
        % Calculate the X and Y coordinates of the points in the adjacent
        % facets.
        p1X = facets(1,minimun-I-loc+1); p1Y = facets(2,minimun-I-loc+1);
        p2X = facets(3,minimun-I-loc+1); p2Y = facets(4,minimun-I-loc+1);

        p3X = facets(1,minimun-I-loc+2); p3Y = facets(2,minimun-I-loc+2);
        p4X = facets(3,minimun-I-loc+2); p4Y = facets(4,minimun-I-loc+2);

        % Check the new peak for an intersection
        x = [p1X p2X];
        y = [p1Y p2Y];
        xS = [p3X p4X];
        yS = [p3Y p4Y];
        [xSc,ySc] = intersectSurface2(x,y,xS,yS);

        % If there is an intersection, adjust the peak
        if ~(isempty(xSc))
%             facets(3,minimun-I-loc+1) = xSc(1,1);
%             facets(4,minimun-I-loc+1) = ySc(1,1);
% 
%             facets(1,minimun-I-loc+2) = xSc(1,1);
%             facets(2,minimun-I-loc+2) = ySc(1,1);
        else
%             internalIndex = 0; % preallocate the internal index
%             for k = 1:2 % 2 because 2 facets are over (the peak has 2 facets)
%                 internalIndex = k;
%                 facets(:,minimun-I-loc+1) = [];
%             end
%             tempIndex = internalIndex + tempIndex;
        end
    I = tempIndex + I; % update global index
    end
end