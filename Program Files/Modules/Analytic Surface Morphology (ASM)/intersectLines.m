function [pX,pY] = intersectLines(x1,y1,x2,y2,x3,y3,x4,y4) %#codegen
%       Calculates the intersection point between two infinite,
%   non-parallel lines that are denoted by two points in space.
    
    % Calculate the common denominator
    D = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
    
    % Coordinate of the x-intersection point
    pX = ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/D;
    
    % Coordinate of the y-intersection point
    pY = ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/D;
    
end