function [ Wq, idx, U, V ] = sphIntegralWeights( xyz, vd )

%approximate integrals over spherical polygons using accurate config: 
%r^5 + poly4 + stencil28, given function values F at center-points of cells

rbforder    = 5;
polyorder   = 5;
stencilSize = 36;

[ idx, h ] = knnsearch( xyz, xyz, 'k', stencilSize );
h = sum(h(:,2:6),2) / 5;

Wq = zeros( size(xyz,1), stencilSize );
U = zeros( size(xyz,1), stencilSize );
V = zeros( size(xyz,1), stencilSize );
for i = 1 : size(xyz,1)
    vv = rotateSphereNodes( xyz(idx(i,:),:), xyz(i,1), xyz(i,2), xyz(i,3) );
    u = vv(:,1) ./ vv(:,3);  v = vv(:,2) ./ vv(:,3);
    U(i,:) = u;  V(i,:) = v;
    uu = meshgrid( u );  vv = meshgrid( v );
    A = 1/h(i)^rbforder * ( (uu.'-uu).^2 + (vv.'-vv).^2 ) .^ (rbforder/2);
    P = [];
    if polyorder > -1
        P = [ P, ones(size(u)) ];
    end
    if polyorder > 0
        P = [ P, [ u, v ] / h(i) ];
    end
    if polyorder > 1
        P = [ P, [ u.^2, u.*v, v.^2 ] / h(i)^2 ];
    end
    if polyorder > 2
        P = [ P, [ u.^3, u.^2.*v, u.*v.^2, v.^3 ] / h(i)^3 ];
    end
    if polyorder > 3
        P = [ P, [ u.^4, u.^3.*v, u.^2.*v.^2, u.*v.^3, v.^4 ] / h(i)^4 ];
    end
    if polyorder > 4
        P = [ P, [ u.^5, u.^4.*v, u.^3.*v.^2, u.^2.*v.^3, u.*v.^4, v.^5 ] / h(i)^5 ];
    end
    A = [ A, P;  P.', zeros(size(P,2)) ];
    lam = A \ [ eye(stencilSize); zeros(size(P,2),stencilSize) ];
    vv = rotateSphereNodes( vd{i}, xyz(i,1), xyz(i,2), xyz(i,3) );
    vertU = vv(:,1) ./ vv(:,3);
    vertV = vv(:,2) ./ vv(:,3);
    B = zeros( 1, length(u)+size(P,2) );
    for j = 1 : length(u)
        B(j) = integrateRbfOverPolygon( h(i), u(j), v(j), vertU, vertV, 1, rbforder );
    end
    if polyorder > -1
        B( length(u)+1 )  = integratePolynomialOverPolygon( h(i), 0, 0, vertU, vertV, 1 );
    end
    if polyorder > 0
        B( length(u)+2 )  = integratePolynomialOverPolygon( h(i), 1, 0, vertU, vertV, 1 );
        B( length(u)+3 )  = integratePolynomialOverPolygon( h(i), 0, 1, vertU, vertV, 1 );
    end
    if polyorder > 1
        B( length(u)+4 )  = integratePolynomialOverPolygon( h(i), 2, 0, vertU, vertV, 1 );
        B( length(u)+5 )  = integratePolynomialOverPolygon( h(i), 1, 1, vertU, vertV, 1 );
        B( length(u)+6 )  = integratePolynomialOverPolygon( h(i), 0, 2, vertU, vertV, 1 );
    end
    if polyorder > 2
        B( length(u)+7 )  = integratePolynomialOverPolygon( h(i), 3, 0, vertU, vertV, 1 );
        B( length(u)+8 )  = integratePolynomialOverPolygon( h(i), 2, 1, vertU, vertV, 1 );
        B( length(u)+9 )  = integratePolynomialOverPolygon( h(i), 1, 2, vertU, vertV, 1 );
        B( length(u)+10 ) = integratePolynomialOverPolygon( h(i), 0, 3, vertU, vertV, 1 );
    end
    if polyorder > 3
        B( length(u)+11 ) = integratePolynomialOverPolygon( h(i), 4, 0, vertU, vertV, 1 );
        B( length(u)+12 ) = integratePolynomialOverPolygon( h(i), 3, 1, vertU, vertV, 1 );
        B( length(u)+13 ) = integratePolynomialOverPolygon( h(i), 2, 2, vertU, vertV, 1 );
        B( length(u)+14 ) = integratePolynomialOverPolygon( h(i), 1, 3, vertU, vertV, 1 );
        B( length(u)+15 ) = integratePolynomialOverPolygon( h(i), 0, 4, vertU, vertV, 1 );
    end
    if polyorder > 4
        B( length(u)+16 ) = integratePolynomialOverPolygon( h(i), 5, 0, vertU, vertV, 1 );
        B( length(u)+17 ) = integratePolynomialOverPolygon( h(i), 4, 1, vertU, vertV, 1 );
        B( length(u)+18 ) = integratePolynomialOverPolygon( h(i), 3, 2, vertU, vertV, 1 );
        B( length(u)+19 ) = integratePolynomialOverPolygon( h(i), 2, 3, vertU, vertV, 1 );
        B( length(u)+20 ) = integratePolynomialOverPolygon( h(i), 1, 4, vertU, vertV, 1 );
        B( length(u)+21 ) = integratePolynomialOverPolygon( h(i), 0, 5, vertU, vertV, 1 );
    end
    Wq(i,:) = B * lam;
end