clear
clc
close all

% The goal of this file is to:
% (1) Document the use of ASTM-E490-00
% (2) Filter/denoise ASTM-E490-00
% (3) Provide a smoothed spectrum that can be used for PMC sampling

% ASTM-E490-00 is an empirical solar spectrum that is publicly available 
% and commonly used in relevant contexts. The raw data can be found at
% this NREL (National Renewable Energy Laboaratory) link:
% https://www.nrel.gov/grid/solar-resource/spectra-astm-e490.html

%% 1. Why is this file needed?
% Currently (4/10/23), the PMC module uses a denoised ASTM-E490-00
% spectrum, but this is not documented well. Furthermore, the "true"
% ASTM-E490-00 spectrum goes from 0.1 to 1000 microns, while the current
% PMC module's spectrum (see irradSpectrum_sunEarth.mat) inexplicably goes
% from 0.24 to 2.6 microns.

% Furthermore, the integrated .mat file produces a spectral
% irradiance of 1309.3 W/m^2, which is ~4% off from the accepted solar
% constant (1366.1 W/m^2) that ASTM-E490-00 was designed to conform with.

% Although not a critical difference, coupled with the lack of
% documentation and truncation of the spectrum, it was decided to redo the 
% denoising process.

%% 1a. A note about solar constant and spectrums
% As mentioned above, the accepted solar constant is 1366.1 W/m^2. However,
% recent findings suggest that it is closer to 1360.8 W/m^2, a 0.39%
% difference. Nevertheless, the 1366.1 value will be used.

% https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010GL045777
%% 2.1: Load in raw ASTM-E490-00 data
raw_data = readmatrix("offical_data_from_NREL.xls");
raw_wavelengths = raw_data(:, 1) / 1e6; % [m]
raw_specIrrad = raw_data(:, 2) * 1e6; % [W/m^2/m]

clc
fprintf("%35s %.3f W/m^2\n", "Accepted solar constant:", 1366.1)
fprintf("%35s %.3f W/m^2\n", "Integrated solar irradiance:", trapz(raw_wavelengths, raw_specIrrad))
%% 2.2: Sanity check -- plot it
figure
hold on
grid on

plot(raw_wavelengths, raw_specIrrad)
xlim([0, 4e-6]) % Spectrum goes to 1e-3 but most is within this range

legend("Raw ASTM-E490-00 Data")

xlabel("Wavelengths [m]")
ylabel("Spectral Irradiance [W/m^2/m]")

%% 2.3: Weighted Moving Average Filter
% As shown in 2.2, the raw spectrum is rather noise. As a result, we will
% filter it to denoise/smooth it out.

% The wavelengths in the raw spectrum are NOT evenly spaced, so a weighted
% moving average filter (as opposed to a "normal" moving average filter)
% must be used.

window = 5; % Number of points left and right of the actual point to average
            % For example, a window of 3 would mean looking at 7 points (3
            % left of the data point, 3 right of the data point, and the
            % point itself)

windowMethod = "MA";
endpointsMethod = "leave as is";

% 1st column is wavelengths, 2nd column is spectral irradiance
smoothedSpectrum = zeros(length(raw_wavelengths), 2);

for i = 1:length(raw_wavelengths)

    smoothedSpectrum(i, 1) = raw_wavelengths(i);

    % Get wavelengths and specIrrad of within-window points
    indices = i-window:i+window;

    % Endpoint condition
    if any(indices < 1 | indices > length(raw_wavelengths))
        if endpointsMethod == "leave as is"  % Don't run filter on endpoints
            smoothedSpectrum(i, 2) = raw_specIrrad(i);
        elseif endpointsMethod == "MA valid points"
            % Run a moving average filter, but only on the valid indices
            valid = indices(indices > 0 & indices <= length(raw_wavelengths));
            smoothedSpectrum(i, 2) = mean(raw_specIrrad(valid));
        end
        continue
    end

    windowWavelengths = raw_wavelengths(indices);
    windowSpecIrrad = raw_specIrrad(indices);

    % Your standard moving average filter
    % Equal weights throughout the entire window (distance independent)
    if windowMethod == "MA"
        smoothedSpectrum(i, 2) = mean(windowSpecIrrad);
    end
end

%% 2.4 Plot smoothed spectrum with raw spectrum
figure
hold on
grid on

plot(raw_wavelengths, raw_specIrrad)
plot(smoothedSpectrum(:, 1), smoothedSpectrum(:, 2), "LineWidth", 1.5)
xlim([0, 1.5e-6]) % Spectrum goes to 1e-3 but biggest differences are within this range

legend(["Raw ASTM-E490-00 Data", "Smoothed Spectrum"])

xlabel("Wavelengths [m]")
ylabel("Spectral Irradiance [W/m^2]")

%% 2.5 Look at integrated solar irradiance
fprintf("\n%35s %.3f W/m^2/m\n", "Using smoothed spectrum:", trapz(smoothedSpectrum(:, 1), smoothedSpectrum(:, 2)))

%% 2.6: Summary
% After playing around with window size and endpoint methods, it appears
% that a window size of *5* (or 11, depending on how you look at it) and an
% endpoint method of "leave as is" is the best choice.

% A window size of 5 actually denoises the spectrum, but appears to still
% maintain most of the important increases/dips within the ultraviolet, 
% visible, and near-infrared infrared.

% A "leave as is" endpoint method likely works best for the right end of
% the spectrum, as the wavelength points become very spaced part. It is
% worth noting that, because the raw spectrum spans such a large range of
% wavelengths, the endpoint method likely doesn't matter.

% This configuration results in a integrated solar irradiance of 1366.5
% W/m^2, which is only 0.4 W/m^2 greater than the accepted value.
%% 3.1: Saving the denoised spectrum

irrad = smoothedSpectrum;  % specIrrad() function assumes that the retrieved spectrum's variable name is "irrad"
save("Smoothed ASTM-E490-00.mat", "irrad")
