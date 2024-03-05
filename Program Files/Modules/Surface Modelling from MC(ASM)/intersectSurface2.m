function [xSc,ySc] = intersectSurface2(x,y,xS,yS) %#codegen
%% DESCRIPTION:
%       Function finds the intersection points between a line and a
%   surface.
%           x = x-coordinates of the line
%           y = y-coordinates of the line
%
%           xS = x-coordinates of the surface
%           yS = y-coordinates of the surface
%
%       Note: xS and yS can be vectors representing the initial and final
%   positions of a second line. In this case, the function will return the
%   intercept between the two lines only if they intercept in the segment
%   described by xS and yS.

%% Calculate the slope of the line and the y-axis interception
    %  Formula for a line:  y = m*x+n
    %               Slope:  m = (y0-y)/(x0-x)
    %           Intercept:  n = y - m*x
    m = (y(1)-y(2))/(x(1)-x(2));
    n = y(1)-m*x(1);

%% Find the intersection points
    tempX = zeros(1,length(xS)-1);  % preallocate temporary x-intercepts
    tempY = zeros(1,length(xS)-1);  % preallocate temporary y-intercepts
    index = 0;% initialize index for the position of the interceptions

    %%% Check every surface segment for interceptions with the line.
    for i = 1:length(xS)-1
        mS = (yS(i)-yS(i+1))/(xS(i)-xS(i+1)); %Slope of current surface segment
        nS = yS(i)-mS*xS(i); %y-axis interception of current surface segment

        % Calculate x and y intercepts between the line and the current segment
        x_intercept = (nS-n)/(m-mS); 
        y_intercept = m*x_intercept + n;

        % Check if the interception point is in the current segment
        if (x_intercept > xS(i) && x_intercept < xS(i+1)) || (x_intercept > xS(i+1) && x_intercept < xS(i))
            % Adjust the index and store the intercept in that position
            index = index + 1;
            tempX(1,index) = x_intercept;
            tempY(1,index) = y_intercept;
        end
    end
    
    %%% Check if interceptions exist and create the final arrays.
    % If the index is zero, it means that there are not intersections
    if index == 0 % if that is the case, return empty interception arrays
        xSc = []; ySc = [];    
    else % if the index is not zero (meaning that there are interceptions)
        xSc = tempX(1,1:index); % X coordinates of the interception points
        ySc = tempY(1,1:index); % Y coordinates of the interception points
    end
end