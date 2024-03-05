function Cp = specificHeat(T,species)
%       Calculates the heat capacity of the species as a function of 
% temperature. T is the input temperature.
    
    % For water
    if species == "H2O"
        % IAPWS formulation (Feistel and Wagner, 2006; IAPWS, 2009)
        
        % Temperature range
        T_ = [-260, -240, -220, -200, -180, -160, -140, -120, -100, -80,...
                           -60, -50, -40, -30, -20, -10, 0] + 273.15; % [K]
        
        % Corresponding specific heat
        Cp_ = [0.036, 0.27, 0.47, 0.65, 0.82, 0.97, 1.11, 1.25, 1.38,...
            1.52, 1.66, 1.73, 1.80, 1.88, 1.95, 2.02, 2.10]*1e3; % [J/Kg-K]
    else
        error('Material type not installed. Select 1) "H2O"')
    end

    % Interpolate the temperature state to find the specific heat
    Cp = interp1(T_,Cp_,T); % [J/Kg-K]

end