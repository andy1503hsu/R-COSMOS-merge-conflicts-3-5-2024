function options = checkUserInputs(options) %#codegen
%       Creates the simulation model according to the user-selected options
%   and checks that all variables are present in the model.

    % Angular resolution
    if exist('options.res', 'var') == 0
        options.res = 2; % default resolution up to a hundredth of a degree
    end

    % Surface modeling method
    if exist('options.declare.adjustByDistance', 'var') == 0
        options.declare.adjustByDistance = false; % step 7 in main function
    end
    
    % Number of iterations
    if exist('options.maxiter', 'var') == 0
        options.maxiter = 2; % default number of iterations is 2
    end
end