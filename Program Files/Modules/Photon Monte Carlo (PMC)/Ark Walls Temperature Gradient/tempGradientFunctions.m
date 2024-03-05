%{
    Creates a temperature-varying function handle with the given parameters
    Parameters:
        Temp1 -- temperature 1
        Temp2 -- temperature 2
        x1 -- First x coordinate
        x2 -- Second x coordinate
        functionType -- the type of function (as a string)
            Currently, functionType can be "linear" or "sinusoidal"
    
    Returns:
        func -- the function handle to the temperature function

    Methohology:
        For a "linear" functionType, a linear function is created between two points, (x1, temp1) and (x2, temp2)

        For a "sinusoidal" functionType, temp1 and temp2 corresponds to the
        minimum and maximum temperatures respectively, while x1 and x2 correspond to two
        different points where (x1, temp1) and (x2, temp1) are valid
        points on the sinusoid.

    Last edited: 9/2/21 by Andy Hsu
 %}

% TODO: Ensure that x1 != x2, maybe change ordering of parameters, ensure
% x2 > x1?
function func = tempGradientFunctions(temp1, temp2, x1, x2, functionType)
        
    % Point slope form, y = y1 + (x - x1)*m
    if functionType == "linear"
        func = @(x) temp1 + (x - x1) * (temp2 - temp1)/(x2-x1);
    
    % Sinusoid is mapped between (x1, temp1), (0.5*(x1 + x2), temp2), and (x2, temp1)
    elseif functionType == "sinusoidal"
        amplitude = (temp2 - temp1)/2;
        vertical_shift = (temp1 + temp2)/2;

        func = @(x) -amplitude * cos(2*pi/(x2 - x1) * (x - x1)) + vertical_shift;
    else
        fprintf("Unknown temperature gradient function type.\n")
    end
end