function t = timePeriod(input)
%       Converts time from a (value1,type1,value2,type2...) syntax input 
% into seconds, or returns the same input if its class type is double. Note 
% that a time input using this class assumes that the input is in seconds.

    % Time unformat
    if isstring(input) % R-COSMOS time format
        t = 0; % [s]
        for ii = 1:length(input)/2
            timeValue = str2double(input(2*ii-1));
            switch input(2*ii)
                case "seconds" % [seconds]
                    t = timeValue + t; % [s]
                case "minutes" % [minutes]
                    t = 6e1*timeValue + t; % [s]
                case "hours" % [hours]
                    t = 3.6e3*timeValue + t; % [s]
                case "days" % [days]
                    t = 8.64e4*timeValue + t; % [s]
                case "years" % [years]
                    t = 3.15576e7*timeValue + t; % [s]
                case "millionYears" % [million years]
                    t = 3.15576e13*timeValue + t; % [s]
                case "billionYears" % [billion years]
                    t = 3.15576e16*timeValue + t; % [s]
            end
        end

    elseif isa(input,'double') % user direct input in seconds
        t = input; % [s]
    end
end