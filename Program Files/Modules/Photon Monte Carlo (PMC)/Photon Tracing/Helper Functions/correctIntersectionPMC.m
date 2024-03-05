% Returns the distance to the next intersection between a photon's current
% trajectory and a snow surface, using intersectSurface.m as a helper
% function

% Cannot use correctIntersection.m from SMC code because that only returns
% the index of the snow surface that a molecule strikes -- in PMC, we are
% concerned about the distance traveled rather than the index
function [distance] = correctIntersectionPMC(intersectX,intersectZ,currentX,currentZ,onSurface,x_traj,z_traj,slope)
    if length(intersectX) == 1
        if intersectX(1) == 0 && intersectZ(1) == 0 %% No intersection was found
            distance = Inf;
            return
        elseif onSurface % Only 1 intersection, and photon *is* on a surface ("Intersection" is really the current location of the photon)
            distance = Inf;
            return
        end
    else % There is more than 1 intersection
        if onSurface % Photon is on a surface -- one of the intersections is a fake one
            % Find fake intersection and delete it
            
            % Minimum distance is the fake intersection -- this will likely be < 1e-10
            [~, minDistanceIndex] = min(((intersectX-currentX).^2 + (intersectZ-currentZ).^2).^0.5);
            intersectX(minDistanceIndex) = [];
            intersectZ(minDistanceIndex) = [];
        end
    end
    
    % Typically, an "eligible" intersection is dependent on sign of the
    % x_trajectory (leftward or rightward). However, if the trajectory has
    % a near infinite slope (> 1e10), then the x-values of intersections
    % are error prone due to floating-point precision error. In these
    % cases, the z_trajectory sign will be used instead
    
    if (abs(slope) < 1e10 && x_traj > 0) || (abs(slope) > 1e10 && z_traj > 0)
        if abs(slope) < 1e10 && x_traj > 0 % Photon is moving rightward
            eligibleIndices = intersectX > currentX;
        else % Photon is moving upward (near-infinite slope case)
            eligibleIndices = intersectZ > currentZ;
        end
        if any(eligibleIndices) % There *are* intersections right of current intersection
            eligibleX = intersectX(eligibleIndices);
            eligibleZ = intersectZ(eligibleIndices);
            [correctX, correctIndex] = min(eligibleX);
        else % There are no intersections right of current intersection
            distance = Inf;
            return
        end
    else
        if abs(slope) < 1e10 && x_traj < 0 % Photon is moving leftward
            eligibleIndices = intersectX < currentX;
        else % Photon is moving downward (near-infinite slope case)
            eligibleIndices = intersectZ < currentZ;
        end
        if any(eligibleIndices) % There *are* intersections left of current intersection
            eligibleX = intersectX(eligibleIndices);
            eligibleZ = intersectZ(eligibleIndices);
            [correctX, correctIndex] = max(eligibleX);
        else % There are no intersections left of current intersection
            distance = Inf;
            return
        end
    end
    
    % We know the correct index corresponding to the minimum distance in
    % the right direction!
    correctZ = eligibleZ(correctIndex);
    
    % correctZ_ = intersectZ(intersectX == correctX); % slower

    % assert(isscalar(correctZ));  % for codegen purposes
    distance = ((correctX-currentX)^2 + (correctZ-currentZ)^2)^0.5;
end