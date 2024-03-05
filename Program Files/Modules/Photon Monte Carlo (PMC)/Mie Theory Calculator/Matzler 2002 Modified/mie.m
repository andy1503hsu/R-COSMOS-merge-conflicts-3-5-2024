% Code is originally taken from Matzler 2002, which can be accessed here:
% https://omlc.org/software/mie/maetzlermie/Maetzler2002.pdf

% The R-COSMOS code only requires a *subset* of the values that the original
% "mie.m" from Matzler 2002 outputs. Consequently, the following "mie.m"
% is a truncated version of Matzler's and only contains the computations
% required for three outputs:
%   Extinction effiency, Qext (qext in Matzler code)
%   Single-scattering albedo, omega  (qsca/qext in Matzler)
%   Asymmetry Parameter, g  (asy in Matzler)

function [qext, omega, asy] = mie(m, x)

    % Computation of Mie Efficiencies for given
    % complex refractive-index ratio m=m'+im"
    % and size parameter x=k0*a, where k0= wave number in ambient
    % medium, a=sphere radius, using complex Mie Coefficients
    % an and bn for n=1 to nmax,
    % s. Bohren and Huffman (1983) BEWI:TDD122, p. 103,119-122,477.
    % Result: m', m", x, efficiencies for extinction (qext),
    % scattering (qsca), absorption (qabs), backscattering (qb),
    % asymmetry parameter (asy=<costeta>) and (qratio=qb/qsca).
    % Uses the function "mie_abcd" for an and bn, for n=1 to nmax.
    % C. MŠtzler, May 2002.

    nmax=round(2+x+4*x^(1/3));
    n1=nmax-1;
    n=(1:nmax);cn=2*n+1; c1n=n.*(n+2)./(n+1); c2n=cn./n./(n+1);
    x2=x*x;
    f=mie_abcd(m,x);
    anp=(real(f(1,:))); anpp=(imag(f(1,:)));
    bnp=(real(f(2,:))); bnpp=(imag(f(2,:)));

    dn=cn.*(anp+bnp);
    q=sum(dn);
    qext=2*q/x2;  % Extinction Effiency

    en=cn.*(anp.*anp+anpp.*anpp+bnp.*bnp+bnpp.*bnpp);
    q=sum(en);
    qsca=2*q/x2;
    omega = qsca/qext; % Single-Scattering Albedo
    
    g1=zeros(4,nmax); % displaced numbers used for
    g1(1,1:n1)=anp(2:nmax); % asymmetry parameter, p. 120
    g1(2,1:n1)=anpp(2:nmax);
    g1(3,1:n1)=bnp(2:nmax);
    g1(4,1:n1)=bnpp(2:nmax);
    asy1=c1n.*(anp.*g1(1,:)+anpp.*g1(2,:)+bnp.*g1(3,:)+bnpp.*g1(4,:));
    asy2=c2n.*(anp.*bnp+anpp.*bnpp);
    asy=4/x2*sum(asy1+asy2)/qsca;  % Asymmetry Parameter
end
