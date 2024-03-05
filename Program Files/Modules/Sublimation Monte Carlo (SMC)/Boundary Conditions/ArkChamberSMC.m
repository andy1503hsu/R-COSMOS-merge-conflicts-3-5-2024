function location = ArkChamberSMC(sData,launchData)
%       Ray-traces particles leaving from the surface inside a domain
%   according to the Ark experiment made at JPL:
%       - Diffusive top boundary condition
%       - 4 cold traps (one at each side of the snow and 2 in the walls)

    [xi,xf,zi,zf,mF,nF,zM,uC] = commonVariables(sData);

    % Launch coordinates and facet
    xL = launchData.xL;    zL = launchData.zL;     iL = launchData.iL;

    % Sine and cosine of the launch angle
    cP = cosd(launchData.Psi);    sP = sind(launchData.Psi);

    % ARK chamber dimensions
    h = 0.44; % [m] height
    w = 1.22; % [m] width

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
                % Isotropic launch angles
                % Theta: Angle w.r.t the normal vector of the facet
                r = rand; % randNumber for Theta
                %  cos(Theta)            sin(Theta)
                cT = (1-r).^0.5;        sT = r.^0.5;

                % Phi: Angle that rotates around the normal vector to the facet
                Phi = 2*pi*rand; % CDF Phi
                sP_ = sin(Phi); % sin(Phi)

                % Psi: Projected 3D angle into the XY plane
                sTsP = sT.*sP_; % constant
                dX = w;
                N = cT*(dX); % numerator
                D = sTsP*(dX); % denominator
                Psi = atan2d(N,D); % projected 2D launch angles
                
                x1 = xTI; % updated initial X coordinate
                z1 = h; % updated initial Z coordinate
                x2 = x1 + cosd(Psi); % updated final X coordinate
                z2 = z1 + sind(Psi); % updated final Z coordinate

                % Updated line properties
                m = (z1-z2)/(x1-x2); % slopes
                n = z1 - m.*x1; % z-intersects

                % Direct surface intersection
                [xSc,ind] = intersectSurface(m,n,xi,xf,zi,zf,mF,nF,zM,uC);
                if ind > 0
                    ind = travelDist(xSc,ind,x1);
                    location(ii) = ind; % arriving location
                    continue
                end
            end
        end

        % Particle's lost
        location(ii) = 0; % lost particle
    end
end