function k = thermalConductivity(~, state, rho, I, species)
%       Calculates the thermal conductivity as a function of Cp (which is a
% function of temperature) at constant density and thermal inertia.

    % Initialize the thermal conductivity matrix
    k = zeros(1,numel(state.u));

    % % % Returning a NaN when time=NaN tells
    % % % the solver that the heat source is a function of time.
    if(isnan(state.time))
        k(1,:) = NaN;
      return
    end

    if I > 0 % If constrained by the thermal inertia
        % Heat capacity
        Cp = specificHeat(state.u,species); % [J/kg-K]

        % Thermal conductivity
        k = I^2./(Cp*rho); % [W/m-K]

    else % if only constrained to density
        % Sturm et al. thermal conductivity from density
        if rho <= 156
            k(1,:) = 0.023 + 0.234e-3.*rho; % [W/m-K]
        else
            k(1,:) = 0.138 - 1.01e-3.*rho + 3.233e-6.*rho^2; % [W/m-K]
        end
    end
end