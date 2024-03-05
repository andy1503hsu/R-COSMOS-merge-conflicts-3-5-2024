function location = openSMC(sData,launchData) %#codegen
%       Ray-traces particles leaving from the surface without any boundary
%   conditions.

    [xi,xf,zi,zf,mF,nF,zM,uC] = commonVariables(sData);

    % Launch coordinates and facet
    xL = launchData.xL;    zL = launchData.zL;     iL = launchData.iL;

    % Sine and cosine of the launch angle
    sP = sind(launchData.Psi);    cP = cosd(launchData.Psi);

    location = zeros(length(xL),1); % preallocate location
    parfor ii = 1:length(xL)
        % Particle's trajectory vector
        x1 = xL(ii); % initial X coordinate
        z1 = zL(ii); % initial Z coordinate
        x2 = x1 + cP(ii); % final X coordinate
        z2 = z1 + sP(ii); % final Z coordinate

        % Line properties
        m = (z1-z2)/(x1-x2); % slopes
        n = z1 - m.*x1; % y-intersects

        % Direct surface intersection
        [xSc,ind] = intersectSurface(m,n,xi,xf,zi,zf,mF,nF,zM,uC);
        % Correct intersection as a function of the launch angle
        ind = correctIntersection(xSc,ind,x1,cP(ii),iL(ii));
        if ind > 0
            location(ii) = ind; % arriving location
            continue
        end

        % Particle's lost
        location(ii) = 0; % lost condition
    end
end