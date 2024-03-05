function location = planetaryWithGravitySMC(sData,launchData,SMC) %#codegen
%       Ray-traces particles leaving from the surface inside a domain with
%   periodic boundary conditions on the sides and open at the top.

    [xi,xf,zi,zf,mF,nF,zM,uC] = commonVariables(sData);

    % Launch coordinates and facet
    xL = launchData.xL;    zL = launchData.zL;     iL = launchData.iL;

    % Sine and cosine of the launch angle
    cP = cosd(launchData.Psi);    sP = sind(launchData.Psi);

    h = max(sData.zi); % height of the domain
    w = sData.xf(end); % width of the domain

    % Mass conservation method
    if ~isfield(SMC,'conservationMethod')
        method = "Gravity"; % default method
        
        % Launch velocity
        vL = launchData.vL; % [m/s]

        % Newtonian constant of gravitation
        G = 6.67430e-11; % [m^3/kg-s^2] 6.67430e-11 (+- 0.00015)

        % Gravity
        g = G*SMC.planetaryMass/SMC.equatorialRadius^2; % [m/s^2]

        % Escape velocity
        eV = (G*SMC.planetaryMass/SMC.equatorialRadius)^0.5; % [m/s]

        % Flight time
        t = 2*launchData.vZ/g;

        % Horizontal displacement
        d = launchData.vX.*t;

    else
        method = SMC.conservationMethod; % overwrite default method
    end

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

        % Initial sine and cosine of the launch angle
        sPsi = sP(ii); % vertical component
        cPsi = cP(ii); % horizontal component

        % Mass conservation process
        while true
        % Top Boundary Condition
        if sPsi > 0
            xTI = ((h-z1)*(x2-x1) + (z2-z1)*x1)/(z2-z1);
            if 0 < xTI && xTI < w
                if method == "Gravity"

                    % Overcoming the escape velocity
                    if vL(ii) > eV % if molecule travels faster than eV
                        location(ii) = 0; % lost condition
                        break
                    end

                    % Gravity induced ballistic path
                    % Re-entering location
                    if d(ii) >= 0 % +X travel
                        xTI = rem(abs(d(ii))+xTI,w); % X intersect
                    elseif d(ii) < 0 % -X travel
                        xTI = w - rem(w-xTI+abs(d(ii)),w); % X intersect
                    end

                    % Re-entering angle
                    sPsi = -sPsi; % vertical component

                elseif method == "Isotropic"

                    % Isotropic launch angles
                    % Theta: Angle w.r.t the normal vector of the facet
                    r = rand; % randNumber for Theta
                    %  cos(Theta)            sin(Theta)
                    cT = (1-r).^0.5;        sT = r.^0.5;
    
                    % Phi: Phase angle
                    Phi = 2*pi*rand; % CDF Phi
                    sP_ = sin(Phi); % sin(Phi)
    
                    % Psi: Projected 3D angle into the XY plane
                    sTsP = sT.*sP_; % constant
                    dX = w;
                    N = cT*(dX); % numerator
                    D = sTsP*(dX); % denominator
                    Psi = -atan2d(N,D); % projected 2D launch angles
                    sPsi = sind(Psi); % vertical component
                    cPsi = cosd(Psi); % horizontal component

                elseif method == "Random"

                    % Random launch angle
                    Psi = -180*rand; % [deg]
                    sPsi = sind(Psi); % vertical component
                    cPsi = cosd(Psi); % horizontal component

                elseif method == "Reflective"

                    sPsi = -sPsi; % vertical component

                else
                    error('Conservation method not installed.');
                end

                x1 = xTI; % updated initial X coordinate
                z1 = h; % updated initial Z coordinate

                x2 = x1 + cPsi; % updated final X coordinate
                z2 = z1 + sPsi; % updated final Z coordinate

                % Updated line properties
                m = (z1-z2)/(x1-x2); % slopes
                n = z1 - m.*x1; % z-intersects

                % Direct surface intersection
                [xSc,ind] = intersectSurface(m,n,xi,xf,zi,zf,mF,nF,zM,uC);
                if ind > 0
                    ind = travelDist(xSc,ind,x1);
                    location(ii) = ind; % arriving location
                    break
                end
            end
        end

        % Periodic boundary conditions
        if cPsi < 0
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
            break
        end

        x2 = x1 + cPsi; % final X coordinate
        z2 = z1 + sPsi; % final Z coordinate

        end
    end
end