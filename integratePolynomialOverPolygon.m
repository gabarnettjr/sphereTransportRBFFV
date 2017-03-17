function I = integratePolynomialOverPolygon( h, m, n, xv, yv, area )

I = moment ( length(xv), xv, yv, m, n );

I = I ./ h^(m+n);

I = I ./ area;