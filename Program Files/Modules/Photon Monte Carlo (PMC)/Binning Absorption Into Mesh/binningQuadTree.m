function [energyInMesh] = binningQuadTree(absorptionEvents, nodes, elements, n)
    if isempty(absorptionEvents)
        energyInMesh = zeros(length(elements), 1);
        return
    end
    triangleIndices = helperFunc(absorptionEvents(:, 1:2), nodes, elements, n);
    energyInMesh = accumarray(triangleIndices,absorptionEvents(:,3),...
                   [length(elements), 1]);
end

function triangleIndices = helperFunc(coords, nodes, elements, n)
% pointLocationQuadTree  Triangle containing specified points
%   This function mimics the behavior of triangulation/pointLocation for 2d
%   triangulations using a quadtree search to accelerate the search for
%   many queryPoints in a large triangulation. If the triangulation is a
%   delaunayTriangulation, it directly uses
%   delaunayTriangulation/pointLocation, which is already very fast.
%
% See also triangulation/pointLocation,
% delaunayTriangulation/pointLocation.
%
% Copyright (C) Sebastian Ullmann 2018

  % We define the clusters according to the location of the query points.
  % This ensures that in the case of n==1, we end up with isolated
  % (multiple) query points, which is easy to detect. What is a good
  % termination criterion if we cluster according to the location of the
  % simplices?
    
    if size(elements,1) * size(coords,1) < n
        triangleIndices = pointLocationStandard(coords, nodes, elements);
        return
    end

  xmean = mean(coords(:, 1));
  ymean = mean(coords(:, 2));
  
  isEast  = @(x) x(:,1) >= xmean;
  isNorth = @(x) x(:,2) >= ymean;
  
  isQueryPointE = isEast(coords);
  isQueryPointN = isNorth(coords);
  
  isQueryPointInQuad = ...
    [( isQueryPointN &  isQueryPointE), ...
     (~isQueryPointN &  isQueryPointE), ...
     (~isQueryPointN & ~isQueryPointE), ...
     ( isQueryPointN & ~isQueryPointE)];
             
  if sum(any(isQueryPointInQuad))<2
    
    [triangleIndices] = pointLocationStandard(coords, nodes, elements);
        
    return
    
  end
  
  
  A = elements(:,1);
  B = elements(:,2);
  C = elements(:,3);
  
  isTriPointE = isEast(nodes);
  isTriPointW = ~isTriPointE;
  isTriPointN = isNorth(nodes);
  isTriPointS = ~isTriPointN;
  
  isTriangleE =  isTriPointE(A) | isTriPointE(B) | isTriPointE(C);
  isTriangleW =  isTriPointW(A) | isTriPointW(B) | isTriPointW(C);
  isTriangleN =  isTriPointN(A) | isTriPointN(B) | isTriPointN(C);
  isTriangleS =  isTriPointS(A) | isTriPointS(B) | isTriPointS(C);
    
  isTriangleInQuad = ...
    [(isTriangleN & isTriangleE), ...
     (isTriangleS & isTriangleE), ...
     (isTriangleS & isTriangleW), ...
     (isTriangleN & isTriangleW)];
   
  triangleIndices = zeros(size(coords,1), 1);
  
  indexes = find(any(isQueryPointInQuad) & any(isTriangleInQuad));
  for i_ = 1:length(indexes)
    i = indexes(i_);
    
    % We do not need to perform any computation if there is not any
    % query point in the i-th quad, because the result will be empty.
    %
    % We do not need to perform any computation if there is not any
    % triangle in the i-th quad, because the result will be NaNs.
    
    indTrianglesInQuad = find(isTriangleInQuad(:,i));
    
    % We create a local triangulation by picking selected triangles from the
    % global triangulation. Therefore we need a local point list and a local
    % triangle list pointing to the indices of the local point list.
    mat = elements(indTrianglesInQuad,:);
    indPointsInQuad = unique(mat(:));
    localPoints = nodes(indPointsInQuad,:);
    localPointMap = zeros(size(nodes,1),1);
    localPointMap(indPointsInQuad,:) = 1:length(indPointsInQuad);
    localTriangles = localPointMap(elements(indTrianglesInQuad,:));
    localTriangles = reshape(localTriangles,length(indTrianglesInQuad),[]);

    % We perform a tree search in the local triangulation.
    
    localQueryPoints = coords(isQueryPointInQuad(:,i),:);
    [localTriangleIndices] = helperFunc(localQueryPoints, localPoints, localTriangles, n);
    % We provide the output triangle indices in terms of indices into the
    % point list of the global triangulation.
    outputTriangleIndices = nan(size(localTriangleIndices));
    isNotNan = ~isnan(localTriangleIndices);
    outputTriangleIndices(isNotNan) = indTrianglesInQuad(localTriangleIndices(isNotNan));
    
    triangleIndices(isQueryPointInQuad(:,i),:) = outputTriangleIndices;
    
  end
end


function [indexes] = pointLocationStandard(absorptionEvents,nodes,elements)
% coder.varsize('absorptionEvents')
    A = elements(:,1);
    B = elements(:,2);
    C = elements(:,3);
    xA = nodes(A,1);
    yA = nodes(A,2);
    xB = nodes(B,1);
    yB = nodes(B,2);
    xC = nodes(C,1);
    yC = nodes(C,2);

    v0x = xB - xA;
    v0y = yB - yA;
    v1x = xC - xA;
    v1y = yC - yA;
    den = v0x .* v1y - v1x .* v0y;
    % They're... always positive??

    indexes = zeros(length(absorptionEvents(:, 1)), 1);

    % for loop
    for i = 1:length(absorptionEvents(:, 1))
        event = absorptionEvents(i, :);

        v2x = event(1) - xA;
        v2y = event(2) - yA;

        v = v2x .* v1y - v1x .* v2y;
        w = v0x .* v2y - v2x .* v0y;
        isInside = v >= 0 & w >= 0 & v + w <= den;

        index = find(isInside, 1);
        if isempty(index) % For whatever reason, this absorption event is NOT in the mesh
            % This is exceedingly rare (maybe floating point error?), but
            % if it happens, bin event into triangle with closest centroid

            % These two lines are NOT pre-calculated due to the exceedingly
            % rare circumstance of this edge case happening
            centroidX = (xA + xB + xC) / 3;
            centroidY = (yA + yB + yC) / 3;

            distFromCentroid = (event(1) - centroidX) .^2 + (event(2) - centroidY) .^ 2;
            [~, index] = min(distFromCentroid);
        end
        indexes(i) = index;

    end
end