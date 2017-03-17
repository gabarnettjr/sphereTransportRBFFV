function I = integrateRbfOverPolygon( h, xc, yc, xv, yv, area, rbforder )

x1 = xv(1) - xc;  y1 = yv(1) - yc;
x2 = xv(2) - xc;  y2 = yv(2) - yc;
x3 = xv(3) - xc;  y3 = yv(3) - yc;

if length(xv) == 3
    
    I =   integratePHSexactlyOverLine( rbforder, x1, y1, x2, y2 ) ...
        + integratePHSexactlyOverLine( rbforder, x2, y2, x3, y3 ) ...
        + integratePHSexactlyOverLine( rbforder, x3, y3, x1, y1 );
    
elseif length(xv) == 4
    
    x4 = xv(4) - xc;  y4 = yv(4) - yc;
    
    I =   integratePHSexactlyOverLine( rbforder, x1, y1, x2, y2 ) ...
        + integratePHSexactlyOverLine( rbforder, x2, y2, x3, y3 ) ...
        + integratePHSexactlyOverLine( rbforder, x3, y3, x4, y4 ) ...
        + integratePHSexactlyOverLine( rbforder, x4, y4, x1, y1 );
    
elseif length(xv) == 5
    
    x4 = xv(4) - xc;  y4 = yv(4) - yc;
    x5 = xv(5) - xc;  y5 = yv(5) - yc;
    
    I =   integratePHSexactlyOverLine( rbforder, x1, y1, x2, y2 ) ...
        + integratePHSexactlyOverLine( rbforder, x2, y2, x3, y3 ) ...
        + integratePHSexactlyOverLine( rbforder, x3, y3, x4, y4 ) ...
        + integratePHSexactlyOverLine( rbforder, x4, y4, x5, y5 ) ...
        + integratePHSexactlyOverLine( rbforder, x5, y5, x1, y1 );
    
elseif length(xv) == 6
    
    x4 = xv(4) - xc;  y4 = yv(4) - yc;
    x5 = xv(5) - xc;  y5 = yv(5) - yc;
    x6 = xv(6) - xc;  y6 = yv(6) - yc;
    
    I =   integratePHSexactlyOverLine( rbforder, x1, y1, x2, y2 ) ...
        + integratePHSexactlyOverLine( rbforder, x2, y2, x3, y3 ) ...
        + integratePHSexactlyOverLine( rbforder, x3, y3, x4, y4 ) ...
        + integratePHSexactlyOverLine( rbforder, x4, y4, x5, y5 ) ...
        + integratePHSexactlyOverLine( rbforder, x5, y5, x6, y6 ) ...
        + integratePHSexactlyOverLine( rbforder, x6, y6, x1, y1 );
    
end

I = I / h^rbforder;

I = I ./ area;