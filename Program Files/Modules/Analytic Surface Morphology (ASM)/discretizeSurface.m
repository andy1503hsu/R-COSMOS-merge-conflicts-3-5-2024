function [xi,xf,yi,yf] = discretizeSurface(xS,yS) %#codegen
%       Uses the coordinates of the surface nodes (xS,yS) to create planar
%   facets denoted by their initial and final x- and y-coordinates.

    % Surface coordinates
    %    Initial                  Final
    xi = xS(1:end-1);        xf = xS(2:end); % x-coordinates [m]
    yi = yS(1:end-1);        yf = yS(2:end); % y-coordinates [m]

end