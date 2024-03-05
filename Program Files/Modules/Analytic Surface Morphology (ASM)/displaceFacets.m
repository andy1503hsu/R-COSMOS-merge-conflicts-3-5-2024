function [Facets,Nodes] = displaceFacets(sData,dS) %#codegen
%       Displaces the facets in the direction of their unit normal vector
%   and stores the information in the structure named 'Facets'.
%   Additionally, this function stores the path of the surface nodes by
%   averaging the displacement of adjacent facets.

    % Number of facets
    nFacets = length(sData.xi);

    % Preallocate the structure containing the displaced facets
    Facets.xi = zeros(1,nFacets); % initial x-coordinates
    Facets.xf = zeros(1,nFacets); % final x-coordinates
    Facets.yi = zeros(1,nFacets); % initial y-coordinates
    Facets.yf = zeros(1,nFacets); % final y-coordinates

    % Adjust the facets based on their displacement and direction
    for ii = 1:length(sData.xi)
        % Calculate the direction of travel
        sT = sind(sData.angleFacet(1,ii)); % sine of theta
        cT = cosd(sData.angleFacet(1,ii)); % cosine of theta

        % Compute and store the X coordinates of the 2 nodes in the facet
        Facets.xi(ii) = sData.xi(ii) + dS(1,ii) * cT; % initial x-coord
        Facets.xf(ii) = sData.xf(ii) + dS(1,ii) * cT; % final x-coord

        % Compute and store the Y coordinates of the 2 nodes in the facet
        Facets.yi(ii) = sData.yi(ii) + dS(1,ii) * sT; % initial y-coord
        Facets.yf(ii) = sData.yf(ii) + dS(1,ii) * sT; % final y-coord
    end

    % Preallocate the structure containing the path of the surface nodes
    nNodes = length(sData.xi)-1;
    Nodes.xi = zeros(1,nNodes); % initial x-coordinates
    Nodes.xf = zeros(1,nNodes); % final x-coordinates
    Nodes.yi = zeros(1,nNodes); % initial y-coordinates
    Nodes.yf = zeros(1,nNodes); % final y-coordinates

    for jj = 1:length(sData.xi)-1        
        % Initial node coordinates
        Nodes.xi(jj) = sData.xi(jj+1); % x-coordinates
        Nodes.yi(jj) = sData.yi(jj+1); % y-coordinates

        % Final node averaged coordinates
        Nodes.xf(jj) = (Facets.xf(jj)+Facets.xi(jj+1))/2; % x-coordinates
        Nodes.yf(jj) = (Facets.yf(jj)+Facets.yi(jj+1))/2; % y-coordinates

        % Travel distance
        Nodes.distance(jj) = sqrt((Nodes.xf(jj)-Nodes.xi(jj))^2 + ...
                                  (Nodes.yf(jj)-Nodes.yi(jj))^2);
    end
end