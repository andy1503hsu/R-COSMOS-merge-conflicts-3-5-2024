function binnedSpectralRadiance = getSpectralRadiance(eData, wavelengthEdges, zenEdges, aziEdges, zenBin, aziBin, temp)
    
    % Bin everything into specified wavelength bins
    wavelengthIndexes = discretize(eData(:, 1), wavelengthEdges);
    binnedEnergy = accumarray(wavelengthIndexes, eData(:, 2), [length(wavelengthEdges)-1, 1]);  % [J]
    normalizedEnergy = binnedEnergy ./ diff(wavelengthEdges)';  % [J/m], normalized by width of wavelength bins

    % Right now, we would be plotting energy per wavelength -- we need to convert to
    % spectral radiance (J/s/m^2/m/sr, Joules per second unit-area wavelength steradian)
    
    binnedPower = normalizedEnergy / temp.time; % Convert from J/m to J/s/m (Watts per wavelength)
    binnedPowerPerSteradian = binnedPower / (-diff(cos(zenEdges(zenBin:zenBin+1))) *diff(aziEdges(aziBin:aziBin+1)) );  % Convert from W/m to W/m/sr
    binnedSpectralRadiance = binnedPowerPerSteradian / max(temp.xS) / cos(mean(zenEdges(zenBin:zenBin+1)));  % Convert from W/m/sr to W/m^2/m/sr
end