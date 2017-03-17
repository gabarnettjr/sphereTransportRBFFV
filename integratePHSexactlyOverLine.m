function I = integratePHSexactlyOverLine( rbforder, x1, y1, x2, y2 )

Tx = x2 - x1;  Ty = y2 - y1;                    %tangent vector
normT = sqrt( Tx.^2 + Ty.^2 );                  %norm of tangent vector
Tx = Tx ./ normT;
Ty = Ty ./ normT;                               %unit tangent vector
x1hat = x1.*Tx + y1.*Ty;                        %first new x coordinate
x2hat = x2.*Tx + y2.*Ty;                        %second new x coordinate
a = sqrt( x1.^2+y1.^2 - x1hat.^2 );             %distance from RBF to line
if rbforder == 1
    I = lineIntegralOfPHS1(x2hat,a) - lineIntegralOfPHS1(x1hat,a);
elseif rbforder == 3
    I = lineIntegralOfPHS3(x2hat,a) - lineIntegralOfPHS3(x1hat,a);
elseif rbforder == 5
    I = lineIntegralOfPHS5(x2hat,a) - lineIntegralOfPHS5(x1hat,a);
end
I = (x1.*y2-y1.*x2)./normT./(rbforder+2) .* I;
I( abs(x1.*y2-y1.*x2) <= eps ) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = lineIntegralOfPHS1(x,a)

% int( (x^2+a^2)^(1/2), x ) :

I = 1/2.*(x.*sqrt(a.^2+x.^2)+a.^2.*log(sqrt(a.^2+x.^2)+x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = lineIntegralOfPHS3(x,a)

% int( (x^2+a^2)^(3/2), x ) :

I = 1/8.*(x.*sqrt(a.^2+x.^2).*(5.*a.^2+2.*x.^2)+3.*a.^4.*log(sqrt(a.^2+x.^2)+x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = lineIntegralOfPHS5(x,a)

% int( (x^2+a^2)^(5/2), x ) :

I = 1/48.*(15.*a.^6.*log(sqrt(a.^2+x.^2)+x)+x.*sqrt(a.^2+x.^2).*(33.*a.^4+26.*a.^2.*x.^2+8.*x.^4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%