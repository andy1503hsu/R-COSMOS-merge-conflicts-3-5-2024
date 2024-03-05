function mieValues = getMieValues(wavelengths, complexIndex, rMean)
    % Retrieves the extinction effiency Qext, single-scattering albedo
    % omega, and asymmetry parameter g using Mie Theory for the inputted
    % wavelengths, complex refractive indices, and grain radius.

    % Note: Both wavelengths and rMean are in meters

    rVariance = 0.05;
    rCount = 5;
    radii = linspace(rMean*(1-rVariance),rMean*(1+rVariance),rCount); % [microns]

    % Preallocate output arrays
    Qext = zeros(length(wavelengths),rCount);
    omega = zeros(length(wavelengths),rCount);
    g = zeros(length(wavelengths),rCount);
    
    % Get Mie values for all grain radii (will average values afterwards)
    for rId = 1:length(radii) % loop through all grain radii
        
        fprintf("%5s Computing coefficients for r = %.2f microns... ","",radii(rId)*1e6);
        
        x = 2*pi*radii(rId)./wavelengths; % [-], size parameter
        
        for wId = 1:length(wavelengths) % loop through all wavelengths

            [Qext_, omega_, g_] = mie(complexIndex(wId), x(wId));

            Qext(wId,rId) = real(Qext_);
            omega(wId,rId) = real(omega_);
            g(wId,rId) = real(g_);
        end
        
        fprintf("done\n");
    end
    
    fprintf("%5s Averaging coefficients across grain radii... ", "")
    %% Remove ripples 1 - Average across grain radii
    Qext = mean(Qext,2);
    omega = mean(omega,2);
    g = mean(g,2);
    
    % Remove any NaN numbers
    validIndices = ~isnan(Qext) & ~isnan(omega) & ~isnan(g);
    wavelengths = wavelengths(validIndices);
    Qext = Qext(validIndices);
    omega = omega(validIndices);
    g = g(validIndices);
    
    %% Remove ripples 2 - Interpolation of Qext 
    wId = 0 <= wavelengths & wavelengths < 2.78e-6; % [m]
    P = polyfit(wavelengths(wId),Qext(wId),2);
    Qext(wId) = polyval(P,wavelengths(wId));
    mieValues = [wavelengths,Qext,omega,g];

    fprintf("done\n")
end