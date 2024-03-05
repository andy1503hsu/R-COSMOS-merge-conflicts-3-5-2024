function s = surfaceArkSnowCake
%       Creates a sinusoidal penitente field with a number of penitentes of
%   a certain width and height, wP and hP, at a given depth with a certain
%   roughness and shifted a certain number of elements to the right.

    % User inputs
    h = 0.1016;             % penitente height
    w = 0.1016;             % penitente width
    d = 2*0.0254;           % depth below trough
    p = 8;                  % number of penitentes
    nP = 20;                % number of surface segments per penitente
    nS = 20;                % number of surface segments per cake side
    t = 0.000;              % frost thickness
    e = 8*0.0254;           % frost length on either side
    nF = 20;                % number of surface segments per frost
    
    
    xx = linspace(0,w*p,p*nP+1); % snow cake x coords
    zz = h/2*sin(xx/w*2*pi-pi/2) + h/2 + d; % snow cake z coords
    
    xx = xx + (1.2-(w*p))/2; % center the snow cake x coords
    
    % Add snow cake side segments
    xxL = xx(1)*ones(1,nS+1); xxR = xx(end)*ones(1,nS+1);
    zzL = linspace(t,zz(1),nS+1); zzR = linspace(zz(end),t,nS+1);
    
    % Add frost to the sides and put it all together
    xx = [linspace(xx(1)-e,xx(1),nF+1),xxL(2:end-1),xx,xxR(2:end-1),...
          linspace(xx(end),xx(end)+e,nF+1)];
    zz = [t*ones(1,nF+1),zzL(2:end-1),zz,zzR(2:end-1),t*ones(1,nF+1)];
    
    s.x = [xx(1),xx,xx(end)];
    s.z = [0    ,zz,0      ];
    
%     plot(s.x,s.z,'r','linewidth',2); hold on; axis equal;
%     scatter(s.x,s.z,40,'k','filled');
end