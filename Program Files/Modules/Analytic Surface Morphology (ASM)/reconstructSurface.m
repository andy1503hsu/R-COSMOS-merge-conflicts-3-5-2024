function [midX,midY] = reconstructSurface(Facets) %#codegen
%       Calculates the average between the adjacent points of two facets to
%   account for rotation and size of the facets.

    % Preallocate the x- and y-coordinates of the reconstructed surface
    midX = zeros(1,length(Facets.xi)-1); % x-coordinates
    midY = zeros(1,length(Facets.xi)-1); % y-coordinates

    % Average the endpoints of adjacent facets
    for ii = 1:length(Facets.xi)-1
        midX(ii) = (Facets.xf(ii) + Facets.xi(ii+1))/2; % x-coordinate
        midY(ii) = (Facets.yf(ii) + Facets.yi(ii+1))/2; % y-coordinate
    end

end