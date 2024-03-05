function dt = timeSingularities(dt,theta,T,tp)
%       Adjusts the elapsed time between a set of two radii at each side of 
% apoapsis or periapsis.

    % True anomalies and time sets
    %   Initial set             final set
    vi = theta(1:end-1);    vf = theta(2:end); % [deg]
    ti = tp(1:end-1);       tf = tp(2:end); % [s] 
    
    % Apoapsis crossing condition
    apogee = vi < 180 & vf > 180;
    dt(apogee) = T - ti(apogee) - tf(apogee);

    % Periapsis crossing condition
    perigee = vi > 180 & vi > vf;
    dt(perigee) = ti(perigee) + tf(perigee);
    
end