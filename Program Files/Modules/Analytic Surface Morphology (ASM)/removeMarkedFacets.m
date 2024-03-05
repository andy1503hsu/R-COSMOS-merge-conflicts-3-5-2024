function Facets = removeMarkedFacets(Facets,ind)
%       Removes the facets marked for elimination according to their index.

    % Remove the facets initial and final x- and y-coordinates.
    %      Initial                     Final
    Facets.xi(ind) = [];        Facets.xf(ind) = []; % x-coordinates
    Facets.yi(ind) = [];        Facets.yf(ind) = []; % y-coordinates

end