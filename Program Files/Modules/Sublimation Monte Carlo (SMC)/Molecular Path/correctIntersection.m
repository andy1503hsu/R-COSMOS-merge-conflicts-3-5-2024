function ind = correctIntersection(xSc,ind,x0,cL,iL) %#codegen
%       Returns the correct intersection between the molecule and the
%   surface based on its initial launch position (LP) and the launch angle.

    if length(xSc) == 1 % for LP only
        ind = 0; % return zero index

    else
        % LP retrieval and storage
        distX = abs(xSc-x0);

        % Index for the minimun distance
        minX = min(distX);
        for iD = 1:length(xSc)
            if distX(iD) == minX
                xSc(iD) = xSc(1);
                ind(iD) = ind(1);
                break
            end
        end
        xSc(1) = x0; % remove floating point error in LP
        ind(1) = iL; % recall index of the LP

        if cL < 0 % for molecules moving to the left of LP
            iLeft = xSc < xSc(1); % intersections to the left of LP
            if any(iLeft) % if there are intersections to the left of LP
                % Store the intersections
                xSc = xSc(iLeft); ind = ind(iLeft);

                % Correct intersection contains the maximun X value
                ind = ind(xSc == max(xSc));

            else % if there are no intersections to the left of LP
                ind = 0; % return zero index
            end

        else % for molecules moving to the right of LP
            iRight = xSc > xSc(1); % intersections to the right of LP
            if any(iRight) % if there are intersections to the right
                % Store the intersections
                xSc = xSc(iRight); ind = ind(iRight);
                
                % Correct intersection contains the minimun X value
                ind = ind(xSc == min(xSc));
                
            else % if there are no intersections to the right of LP
                ind = 0; % return zero index 
            end
        end
    end
end