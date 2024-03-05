function Facets = intersectFacets(Facets) %#codegen
%       Determines if there is any intersection between adjacent facets and
%   adjusts the facets that are intersecting.

    % Create a temporary copy of the facets structure in order to modify
    % the original version while maintaining consistency.
    TempFacets = Facets; % local temporary structure

    % Run through each facet
    for ii = 1:length(Facets.xi)-1
        % Retrieve the initial and final coordinates of the current facet.
        %    Initial                    Final
        x1 = TempFacets.xi(ii);    x2 = TempFacets.xf(ii); % x-coordinate
        y1 = TempFacets.yi(ii);    y2 = TempFacets.yf(ii); % y-coordinate

        % Calculate the slope of the current facet and adjust it if it is
        % undefined.
        mCurrent = (y2-y1)./(x2-x1); % current facet slope
        mCurrent(mCurrent == -Inf) = Inf; % for vertical line segments

        % Retrieve the initial and final coordinates of the next facet.
        %    Initial                    Final
        x3 = TempFacets.xi(ii+1);  x4 = TempFacets.xf(ii+1); % x-coordinate
        y3 = TempFacets.yi(ii+1);  y4 = TempFacets.yf(ii+1); % y-coordinate

        % Calculate the slope of the current facet and adjust it if it is
        % undefined.
        mNext = (y4-y3)./(x4-x3); % next facet slope
        mNext(mNext == -Inf) = Inf; % for vertical line segments

        % Find the intersection point between both facets.
        [pX,pY] = intersectLines(x1,y1,x2,y2,x3,y3,x4,y4);

        % The following section accounts for floating point error and
        % undefined conditions.
        % Zero slopes (horizontal facets)
        if mCurrent == 0
            pY = TempFacets.yi(ii); % Y coordinates (for zero-slope)
        end

        % Zero slopes (horizontal facets)
        if mNext == 0
            pY = TempFacets.yi(ii+1); % Y coordinates (for zero-slope)
        end

        % Undefined conditions (vertical facets)
        if mCurrent == Inf
            pX = TempFacets.xi(ii); % X intersect coordinates
        end

        % Undefined conditions (vertical facets)
        if mNext == Inf
            pX = TempFacets.xi(ii+1); % X intersect coordinates
        end

        % Determine if the intersection point lies in both facets.
        if (x2-pX)*(x1-pX) <= 0 && (y2-pY)*(y1-pY) <= 0
            if (x4-pX)*(x3-pX) <= 0 && (y4-pY)*(y3-pY) <= 0
                % Adjust the final coordinates of the current facet to lie
                % on the intersection point.
                Facets.xf(ii) = pX; % final x-coordinate
                Facets.yf(ii) = pY; % final y-coordinate

                % Adjust the initial coordinates of the next facet to lie
                % on the intersection point.
                Facets.xi(ii+1) = pX; % initial x-coordinate
                Facets.yi(ii+1) = pY; % initial y-coordinate
            end
        end
    end
end