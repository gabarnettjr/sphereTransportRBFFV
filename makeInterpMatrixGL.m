function [ W, xGL, yGL, zGL, weightsGL, U, V, idx, areas ] = makeInterpMatrixGL( ...
    xyz, edges, vd, ...
    rbforder, polyorder, stencilSize, nGL, ...
    xv1, yv1, zv1, xv2, yv2, zv2 )

[ idx, h ] = knnsearch( xyz, edges, 'k', stencilSize );
h = h(:,1) + h(:,2);

ne = size( edges, 1 );                              %number of edges
np = ( polyorder + 1 ) * ( polyorder + 2 ) / 2;     %number of polynomials

W = zeros( ne, stencilSize, nGL );

xGL = zeros( ne, nGL );
yGL = zeros( ne, nGL );
zGL = zeros( ne, nGL );
U   = zeros( ne, nGL );
V   = zeros( ne, nGL );

weightsGL = zeros( ne, nGL );

areas = zeros( ne, stencilSize );

for i = 1 : ne
    [ xyz, S, T ] = rotateSphereNodes( [xv1(i),yv1(i),zv1(i);xv2(i),yv2(i),zv2(i)], ...
        edges(i,1), edges(i,2), edges(i,3) );
    u = xyz(:,1) ./ xyz(:,3);  v = xyz(:,2) ./ xyz(:,3);
    [ uGL, vGL, weightsGL(i,:) ] = getGLnodesAndWeights( u(1), v(1), u(2), v(2), nGL );
    U(i,:) = uGL;  V(i,:) = vGL;
%     %plot stencil:
%     figure(1),clf
%     sphere
%     shading( 'interp' )
%     hold( 'on' )
%     plot3( uGL, vGL, ones(size(uGL)), 'r.' )
%     %scaling factor so that GL weights will work properly on the sphere:
%     tmp = ( 1 + uGL.^2 + vGL.^2 ) .^ (3/2);
%     weightsGL(i,:) = weightsGL(i,:) .* tmp.';
    %project GL nodes back down to the sphere:
    xyzGL(:,1) = uGL ./ sqrt( 1 + uGL.^2 + vGL.^2 );
    xyzGL(:,2) = vGL ./ sqrt( 1 + uGL.^2 + vGL.^2 );
    xyzGL(:,3) = sqrt( 1 - xyzGL(:,1).^2 - xyzGL(:,2).^2 );
    %rotate GL nodes back to where they would have started:
    xyzGL = xyzGL * S;
    xyzGL = xyzGL / T;
    xyzGL = xyzGL / S;
    xGL(i,:) = xyzGL(:,1);
    yGL(i,:) = xyzGL(:,2);
    zGL(i,:) = xyzGL(:,3);
    %initialize matrices for calculating RBF-poly interpolation weights:
    A = zeros( stencilSize, stencilSize );
    P = zeros( stencilSize, np );
    B = zeros( nGL, stencilSize );
    C = [ ones(nGL,1), zeros( nGL, np-1 ) ];
    for j = 1 : stencilSize
        xyz = vd{ idx(i,j) };
        xyz = rotateSphereNodes( xyz, edges(i,1), edges(i,2), edges(i,3) );
        u = xyz(:,1) ./ xyz(:,3);
        v = xyz(:,2) ./ xyz(:,3);
        uc = repmat( sum(u)/length(u), nGL, 1 );
        vc = repmat( sum(v)/length(v), nGL, 1 );
        B(:,j) = phi( h(i), uGL-uc, vGL-vc, rbforder, [] );
        areas(i,j) = moment( length(u), u, v, 0, 0 );
        for k = 1 : stencilSize
            xyz = vd{ idx(i,k) };
            xyz = rotateSphereNodes( xyz, edges(i,1), edges(i,2), edges(i,3) );
            tmpu = xyz(:,1) ./ xyz(:,3);
            tmpv = xyz(:,2) ./ xyz(:,3);
            uc = sum(tmpu) ./ length(tmpu);
            vc = sum(tmpv) ./ length(tmpv);
            A(j,k) = integrateRbfOverPolygon( h(i), uc, vc, u, v, areas(i,j), rbforder );
            
%             xyz = [ u, v, ones(size(u)) ];
%             plot3( xyz([1:end,1],1), xyz([1:end,1],2), xyz([1:end,1],3), 'k' )
%             plot3( uc, vc, 1, 'ko' )
        end
        if polyorder > -1
            P(j,1) = integratePolynomialOverPolygon( h(i), 0, 0, u, v, areas(i,j) );
        end
        if polyorder > 0
            P(j,2) = integratePolynomialOverPolygon( h(i), 1, 0, u, v, areas(i,j) );
            P(j,3) = integratePolynomialOverPolygon( h(i), 0, 1, u, v, areas(i,j) );
        end
        if polyorder > 1
            P(j,4) = integratePolynomialOverPolygon( h(i), 2, 0, u, v, areas(i,j) );
            P(j,5) = integratePolynomialOverPolygon( h(i), 1, 1, u, v, areas(i,j) );
            P(j,6) = integratePolynomialOverPolygon( h(i), 0, 2, u, v, areas(i,j) );
        end
        if polyorder > 2
            P(j,7)  = integratePolynomialOverPolygon( h(i), 3, 0, u, v, areas(i,j) );
            P(j,8)  = integratePolynomialOverPolygon( h(i), 2, 1, u, v, areas(i,j) );
            P(j,9)  = integratePolynomialOverPolygon( h(i), 1, 2, u, v, areas(i,j) );
            P(j,10) = integratePolynomialOverPolygon( h(i), 0, 3, u, v, areas(i,j) );
        end
        if polyorder > 3
            P(j,11) = integratePolynomialOverPolygon( h(i), 4, 0, u, v, areas(i,j) );
            P(j,12) = integratePolynomialOverPolygon( h(i), 3, 1, u, v, areas(i,j) );
            P(j,13) = integratePolynomialOverPolygon( h(i), 2, 2, u, v, areas(i,j) );
            P(j,14) = integratePolynomialOverPolygon( h(i), 1, 3, u, v, areas(i,j) );
            P(j,15) = integratePolynomialOverPolygon( h(i), 0, 4, u, v, areas(i,j) );
        end
    end
    if polyorder > 0
        C(:,2) = uGL / h(i);
        C(:,3) = vGL / h(i);
    end
    if polyorder > 1
        C(:,4) = uGL.^2           / h(i)^2;
        C(:,5) = uGL    .* vGL    / h(i)^2;
        C(:,6) =           vGL.^2 / h(i)^2;
    end
    if polyorder > 2
        C(:,7)  = uGL.^3           / h(i)^3;
        C(:,8)  = uGL.^2 .* vGL    / h(i)^3;
        C(:,9)  = uGL    .* vGL.^2 / h(i)^3;
        C(:,10) =           vGL.^3 / h(i)^3;
    end
    if polyorder > 3
        C(:,11) = uGL.^4           / h(i)^4;
        C(:,12) = uGL.^3 .* vGL    / h(i)^4;
        C(:,13) = uGL.^2 .* vGL.^2 / h(i)^4;
        C(:,14) = uGL    .* vGL.^3 / h(i)^4;
        C(:,15) =           vGL.^4 / h(i)^4;
    end
    w = [ B, C ] * ( [ A, P;  P.', zeros(np) ] \ [ eye(stencilSize); zeros(np,stencilSize) ] );
    for j = 1 : nGL
        W(i,:,j) = w(j,:);
    end
    
%     axis( 'equal' )
%     drawnow,pause
end