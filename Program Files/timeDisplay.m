function timeDisplay(inputTime)
%       Displays time in time categories for different time periods.

    % Simulation time
    inputTime = timePeriod(inputTime); % generate time period
    sec    = inputTime;  % [seconds]
    min    = 6e1;               % [minutes]
    hours  = 3.6e3;             % [hours]
    days   = 8.64e4;            % [days]
    years  = 3.15576e7;          % [years]
    Myears = 3.15576e13;         % [million years]
    Byears = 3.15576e16;         % [bilion years]

    % Based on the time categories selected, find the remainder from the 
    % model input (in seconds) and use it to calculate the value of the
    % subsequently smaller time category
    tB     = floor(sec/Byears);     % billion years left
    remB   = rem(sec, Byears);      % remainder from billion years
    tM     = floor(remB/Myears);    % million years left
    remM   = rem(remB, Myears);     % remainder from million years
    ty     = floor(remM/years);     % years left
    remY   = rem(remM, years);      % remainder from years
    td     = floor(remY/days);      % days left
    remD   = rem(remY, days);       % remainder from days
    th     = floor(remD/hours);     % hours left
    remH   = rem(remD, hours);      % remainder from hours
    tm     = floor(remH/min);       % minutes left
    remM   = rem(remH, min);        % remainder from minutes
    ts     = round(remM);           % seconds left

    % Display the projected simulation time in the time categories
    fprintf('%dB %dM %dy %dd %dh %dm %ds\n', tB, tM, ty, td, th, tm, ts)

end