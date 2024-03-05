function N = modelingPrecision(Surface,displacement) %#codegen
%       Measures the precision of the simulation according to the average 
%   area of the surface facets and the corresponding displacement.
    
    % Discretize the surface nodes into facets
    [xi,xf,yi,yf] = discretizeSurface(Surface.xS,Surface.yS);
    
    % Calculate the mean area of the facets
    area = mean(calculateFacetsArea(xi,xf,yi,yf));
    
    % Estimate the precision of the simulation using the dimensionless 
    % parameter N = 100(d^2/A)
    N = 100 * mean(displacement)^2 / area;
    
    % Print the results
    fprintf('      Estimated value of N: %.2f\n',N);
    
end