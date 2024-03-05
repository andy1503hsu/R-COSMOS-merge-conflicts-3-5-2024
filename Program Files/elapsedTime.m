function dt = elapsedTime(radii,theta,T,sma,ecc,GMsun)
%       Calculates the elapsed time since periapsis at every point in the
% orbit (including at periapsis).

    % Eccentric anomaly
    E = acos((sma-radii)/(ecc*sma)); % [rad]

    % Kepler's equation
    M = E - ecc*sin(E);

    % Time since periapsis
    tp = sqrt(sma^3/GMsun)*M; % [s]

    % Time at periapsis
    tp(imag(E) ~= 0) = 0; % [s]

    % Elapsed time between subsequent radii
    dt = abs(tp(2:end)-tp(1:end-1)); % [s]

    % Time singularities
    dt = timeSingularities(dt,theta,T,tp); % [s]

    % Circular orbits
    if ecc == 0
        dt(:) = T/length(dt); % orbital period
    end

end