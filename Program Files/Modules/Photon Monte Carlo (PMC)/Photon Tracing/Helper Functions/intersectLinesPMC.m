function [pX,pY] = intersectLinesPMC(x1,y1,x2,y2,x34,y34,determinant34) %#codegen
%       Calculates the intersection point between two infinite,
%   non-parallel lines that are denoted by two points in space.
    

    %x34 = xi-xf;
    %z34 = zi-zf; --> called y34 in this function
    %determinant34 = xi.*zf-zi.*xf;
    % Calculate the common denominator
    invD = 1./((x1-x2)*y34-(y1-y2)*x34);
    determinant12 = x1*y2-y1*x2;
    
    pX = (determinant12.*x34-(x1-x2)*determinant34).*invD;
    pY = (determinant12.*y34-(y1-y2)*determinant34).*invD; 
    
end