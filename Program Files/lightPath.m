function path = lightPath(orbit,model)
%       Determines the path of the light source. For Planetary simulations
% the light path is adjusted depending on the location in the orbit, the
% axial tilt of the planetary body and the orientation of the surface
% locally.

    % Day, orbit steps, and rotation period
    dS = model.resolution.daySteps;
    oS = model.resolution.orbitSteps;
    dL = rotationPeriod(model.body.rotation);

    % Generic light source path for E-W orientation at 0 [deg] definition
    t = linspace(0, 360, dS+1); % angular distribution
    u = [zeros(1,length(t));... x-coordinates
                  -cosd(t);... y-coordinates
                  sind(t)];... z-coordinates
%     u = [zeros(1,length(t));... x-coordinates
%                    -sind(t);... y-coordinates
%                    -cosd(t)];... z-coordinates

    % Rotation matrix around z axis
    Rz = @(phi)([cosd(phi) sind(phi) 0;-sind(phi) cosd(phi) 0;0 0 1]);
    
    % Local surface orientation angle w.r.t North
    orientation = surfaceOrientation(model.surface.orientation);

    % Light source path by the surface orientation
    ui = Rz(orientation)*u;

    % New body x-axis
    v = ui(:,1); % <xi,yi,zi>

    % Rotation matrix about new body x-axis
    vx = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0]; % skew symmetrix matrix
    Rv = @(a)(cosd(a)*eye(3) + sind(a)*vx + (1-cosd(a))*(v*v'));

    % Orbital factor from true anomaly
    d = sind(orbit.lFactor); % [~]

    % Zenith, azimuth, and cartesian light path
    path.zenith = zeros(dS+1, oS); % [deg] zenith
    path.azimuth = zeros(dS+1, oS); % [deg] azimuth
    path.lightPath = zeros(3, dS+1, oS); % light source path
    path.dayLength = zeros(1,oS); % [s] length of day

    if model.gravity
        % Inclination angle at the latitude and current orbital location
        phi = model.surface.latitude + orbit.lFactor; % [deg]
    end

    % Light source path for each orbital period
    for ii = 1:oS % for each orbital period
        if model.gravity % for mass conservation
            % Light source path at the current orbit location
            path.lightPath(:,:,ii) = Rv(phi(ii))*ui; % [~]

        else % for regions of net sublimation
            % Original light source path definition from E-W orientation
            lightPath = u; % defined on the y-z plane
    
            % Adjusted radius of circular path
            r = 1-abs(d(ii));
    
            % Scaling of the light source path
            lightPath([2,3],:) = rescale(lightPath([2,3],:),-r,r);
    
            % Adjusted light path from orientation (rotation about z axis)
            lightPath = Rz(orientation)*lightPath;
    
            % Light path at the current orbit step
            lightPath(1,:) = lightPath(1,:) - d(ii)*cosd(orientation);
            lightPath(2,:) = lightPath(2,:) + d(ii)*sind(orientation);

            % Light path adjusted for latitude and stored in path struct
            path.lightPath(:,:,ii) = Rv(model.surface.latitude)*lightPath;
        end

        % Length of day [s]
        path.dayLength(ii) = sum(path.lightPath(3,:,ii)>0)/(dS+1)*dL;
        
        % Zenith
        path.zenith(:,ii) = acosd(path.lightPath(3,:,ii)); % [deg]

        % Azimuth
        path.azimuth(:,ii) = atan2d(path.lightPath(2,:,ii),...
                                        path.lightPath(1,:,ii)); % [deg]
    end
end