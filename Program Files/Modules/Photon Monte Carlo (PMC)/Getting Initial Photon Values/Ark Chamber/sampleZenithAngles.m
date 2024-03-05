function zenith_angles = sampleZenithAngles(numPhotons, LEDType)

    % NIR_CDF(1,:) contains zenith angles (in degrees)
    % NIR_CDF(2,:) corresponds to CDF at these angles
    if LEDType == "NIR"
        LED_CDF = load("Zenith Angle CDF of NIR and Visible LEDs.mat", 'NIR_CDF');
    else % LEDType == "Visible"
        LED_CDF = load("Zenith Angle CDF of NIR and Visible LEDs.mat", 'Visible_CDF');
    end
    
    rand_nums = rand(1, numPhotons);
    bins = discretize(rand_nums, LED_CDF(2, :));
    
    % Do psuedo-interpolation
    % If zen = 30 degrees has a CDF of 0.20 and zen = 31 has a CDF of 0.22,
    % then a random number of 0.21 will correspond to a zenith angle of 30.5
    % There are sufficient points where the "distance" in-between points
    % can be assumed to be linear
    
    % Guaranteed to be positive or zero, since rand_nums >= left_edge of
    % relevant bin (otherwise discretize() did something wrong)
    CDF_from_left_edge = rand_nums - LED_CDF(2, bins);
    
    % (1,:) stores the delta Zenith Angle between CDF points
    % (2,:) stores the PDF within each corresponding bin
    PDF_each_bin = LED_CDF(:, 2:end) - LED_CDF(:, 1:end-1);
    
    % Get zenith angles using psuedo-interpolation strategy
    zenith_angles = LED_CDF(1, bins) + PDF_each_bin(1, bins) .* (CDF_from_left_edge ./ PDF_each_bin(2, bins));
    zenith_angles = zenith_angles'; % Change into column vector
end