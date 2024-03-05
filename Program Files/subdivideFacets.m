function [xS,zS] = subdivideFacets(x,z,s)
%
    % Surface coordinates
    %    Initial                      Final
    xi = x(1:end-1);        xf = x(2:end); % temporary X [m]
    zi = z(1:end-1);        zf = z(2:end); % temporary Z [m]

    % Slopes and y-intercepts of the facets (formula for a line)
    %    Formulas                     SMC definitions
    mF = (zf-zi)./(xf-xi);      mF(mF == -Inf) = Inf; % slopes
    nF = zf - mF.*xf;           nF(isnan(nF)) = -Inf; % y-intercepts

    xS = xi(1);
    zS = zi(1);
    for ii = 1:length(x)-1
        if mF == Inf
            X = xi(ii)*ones(1,s+2);
            Z = linspace(zi(ii),zf(ii),s+2);

        else
            X = linspace(xi(ii),xf(ii),s+2);
            Z = mF(ii)*X + nF(ii);
        end

        X(1) = [];
        Z(1) = [];

        xS = [xS X];
        zS = [zS Z];

    end
end