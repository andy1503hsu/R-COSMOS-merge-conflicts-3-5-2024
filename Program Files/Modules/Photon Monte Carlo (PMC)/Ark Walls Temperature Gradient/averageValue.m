%{
    Given a function handle and a vector of lower bounds and upper bounds, calculates the average
    value this function over each of those ranges

    Parameters:
        func -- A function handle that the average value(s) will be calculated with
        lowerBounds -- A vector that corresponds to the lower bounds of the
        specified ranges where the average value will be calculated "on"
        upperBounds -- Same logic as lowerBounds but for upper bounds

        lowerBounds and upperBounds should have the same size! (1 x n)

    Returns:
        avgValues -- A vector of average values, will have size of 1 x (n-1)

    Last edited: 9/1/21 by Andy Hsu
%}

function avgValues = averageValue(func, lowerBounds, upperBounds)
    avgValues = zeros(1, length(lowerBounds));
    for i = 1:length(lowerBounds)
        
        % Trapezoidal integration
        % Would use integral() but this is not compilable
        avgValue = 0.5*(func(upperBounds(i)) - func(lowerBounds(i)))*(upperBounds(i) - lowerBounds(i));
        avgValues(1,i) = avgValue;
    end
end