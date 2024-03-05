function orientation = surfaceOrientation(input)
%       Retrieves the orientation of the surface given the user input. The
% surface is defined at 0 [deg] North.

    % Surface orientation
    if isstring(input) % planets
        switch input % predetermined inputs
            case "E-W"
                orientation = 0; % [deg] latitude
            case "S-N"
                orientation = 90; % [deg] latitude
            case "W-E"
                orientation = 180; % [deg] latitude
            case "N-S"
                orientation = 270; % [deg] latitude

            otherwise
                error('Only N-S and E-W pairs are available as a string.')
        end

    elseif isa(input,'double') % surface orientation input by user
        orientation = input; % [deg] latitude

    else
        error('Input type %s is not supported.', class(input))
    end

end