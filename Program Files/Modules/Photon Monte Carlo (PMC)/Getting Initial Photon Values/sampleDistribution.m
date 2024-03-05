function [xxS,xxW] = sampleDistribution(xx,yy,xLim,xxN)
%       Sample a total of xxN values from the distribution described by 
%   datapoints (xx,yy) within the x interval xLim(1) <= x <= xLim(2). The
%   distribution is assumed to be continuous (rather than discrete),
%   contain non-negative values, and have at least one positive value.

    xL = xLim(1); xU = xLim(2); % interval bounds

    xx = xx(:); % force xxs to be a column vector
    yy = yy(:); % force yys to be a column vector

    % [Validation] Plot the data and interval bounds
%     plot(xx,yy,'--r','LineWidth',2); hold on;
%     plot([xL xL],[0 max(yy)],':k','LineWidth',2);
%     plot([xU xU],[0 max(yy)],':k','LineWidth',2);
%     plot([xL xU],[0 0],':k','LineWidth',2);

    % Truncate the data to fit within interval bounds
    if all(xL <= xx & xx <= xU) % if all xxs are in or on bounds
        
        % do nothing, we'll use all the data
        
    elseif all(xx <= xL | xx >= xU) % if all xxs are out of bounds
        
        xxS = []; % return zero sample array
        xxW = 0; % each sample has no weight
        return;

    elseif all(xL <= xx) && any(xU <= xx) % if some xxs are above xU

        yU = interp1(xx,yy,xU); % find the y value for xU

        % truncate xx up to xU
        id = xx < xU;
        xx = [xx(id);xU];
        yy = [yy(id);yU];

    elseif any(xx <= xL) && all(xx <= xU) % if some xxs are below xL

        yL = interp1(xx,yy,xL); % find the y value for xL

        % truncate xx down to xL
        id = xL < xx;
        xx = [xL;xx(id)];
        yy = [yL;yy(id)];

    elseif any(xx <= xL) && any(xU <= xx) % if some xxs are beyond xL & xU

        yL = interp1(xx,yy,xL); % find the y value for xL
        yU = interp1(xx,yy,xU); % find the y value for xU

        % truncate xx down to xL and up to xU
        id = xL < xx & xx < xU;
        xx = [xL;xx(id);xU];
        yy = [yL;yy(id);yU];
    end

    yyC = cumtrapz(xx,yy); % find the cummulative integral
    yyNet = yyC(end); % net integral
    
    % Sample yy values and find corresponding xx values
    yyS = rand(xxN,1)*yyNet; % uniformly sampled yy values
    id = discretize(yyS,yyC); % line segments where yy values fall into
    mm = diff(yyC)./diff(xx); % slope of each line
    xxS = (yyS - yyC(id))./mm(id) + xx(id); % corresponding xx values
    xxW = yyNet/xxN; % the weight of each xx value
    
    % [Validation] Plot in a histogram
%     [N,edges] = histcounts(xxS); % bin the data
%     histogram('BinEdges',edges,'BinCounts',N*xxW,...
%                 'Normalization','CountDensity',...
%                 'DisplayStyle','Stairs',...
%                 'LineWidth',2,...
%                 'EdgeColor','g'); hold on
%     plot(xx,yy);
end

