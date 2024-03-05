function minIndex = travelDist(xSI,ind,x0) %#codegen
%       Calculates the distance from the initial launch position (x0,y0) to
%   all the intersection points (xSI,ySI). Then, it swaps the location of 
%   the minimun distance and the first location in the array.

    % Minimun distance from x0 to xSI
    distX = abs(xSI-x0);
    minX = min(distX);

    % Index for the minimun distance and allocation
    minIndex = 0;
    for iD = 1:length(xSI)
        if distX(iD) == minX
            minIndex = ind(iD);
            break
        end
    end
end