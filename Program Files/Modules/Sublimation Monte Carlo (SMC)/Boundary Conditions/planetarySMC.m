function location = planetarySMC(sData,launchData) %#codegen
%       Ray-traces particles leaving from the surface inside a domain with
%   periodic boundary conditions on the sides and open at the top.

    [xi,xf,zi,zf,mF,nF,zM,uC] = commonVariables(sData);

    % Launch coordinates and facet
    xL = launchData.xL;    zL = launchData.zL;     iL = launchData.iL;

    % Sine and cosine of the launch angle
    cP = cosd(launchData.Psi);    sP = sind(launchData.Psi);

    h = max(sData.zi); % height of the domain
    w = sData.xf(end); % width of the domain

    location = zeros(length(xL),1); % preallocate location
    parfor ii = 1:length(xL)
        % Particle's trajectory vector
        x1 = xL(ii); % initial X coordinate
        z1 = zL(ii); % initial Z coordinate
        x2 = x1 + cP(ii); % final X coordinate
        z2 = z1 + sP(ii); % final Z coordinate

        % Line properties
        m = (z1-z2)/(x1-x2); % slopes
        n = z1 - m.*x1; % z-intersects

        % Direct surface intersection
        [xSc,ind] = intersectSurface(m,n,xi,xf,zi,zf,mF,nF,zM,uC);
        % Correct intersection as a function of the launch angle
        ind = correctIntersection(xSc,ind,x1,cP(ii),iL(ii));
        if ind > 0
            location(ii) = ind; % arriving location
            continue
        end

        % Top Boundary Condition
        if sP(ii) > 0
            xTI = ((h-z1)*(x2-x1) + (z2-z1)*x1)/(z2-z1);
            if 0 < xTI && xTI < w
                location(ii) = 0; % lost particle
                continue
            end
        end

        % Periodic boundary conditions
        if cP(ii) < 0
            yLI = (z2*x1 - z1*x2)/(x1-x2); % Z left intersect
            x1 = w; % updated initial X coordinate
            z1 = yLI; % updated initial Z coordinate
        else
            yRI = (w-x1)*(z2-z1)/(x2-x1) + z1; % Z right intersect
            x1 = 0; % updated initial X coordinate
            z1 = yRI; % updated initial Z coordinate
        end
        n = z1 - m.*x1; % updated line's z-intersect

        % Direct surface intersection
        [xSc,ind] = intersectSurface(m,n,xi,xf,zi,zf,mF,nF,zM,uC);
        if ind > 0
            ind = travelDist(xSc,ind,x1);
            location(ii) = ind; % arriving location
            continue
        end

        % Particle's lost
        location(ii) = 0; % lost condition
    end
end