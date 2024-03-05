function area = calculateFacetsArea(xi,xf,yi,yf) %#codegen
%       Calculates the area of the facets given their initial and final x-
%   and y-coordinates.

    % Area of the facets (distance formula)
    area = ((xf-xi).^2 + (yf-yi).^2).^0.5;

end