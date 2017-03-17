function [ xGL, yGL, wGL ] = getGLnodesAndWeights( x1, y1, x2, y2, nGL )

midpt = [ (x1+x2)/2, (y1+y2)/2 ];

tangent = [x2,y2] - [x1,y1];

nodesGL = zeros( nGL, 2 );

if nGL == 5

    nodesGL(1,:) = midpt - .9061798459 * tangent/2;
    nodesGL(5,:) = midpt + .9061798459 * tangent/2;

    nodesGL(2,:) = midpt - .5384693101 * tangent/2;
    nodesGL(4,:) = midpt + .5384693101 * tangent/2;

    nodesGL(3,:) = midpt;

    wGL = [ .2369268851; .4786286705; .5688888889; .4786286705; .2369268851 ];
    
elseif nGL == 4
    
    nodesGL(1,:) = midpt - .8611363116 * tangent/2;
    nodesGL(4,:) = midpt + .8611363116 * tangent/2;
    
    nodesGL(2,:) = midpt - .3399810436 * tangent/2;
    nodesGL(3,:) = midpt + .3399810436 * tangent/2;
    
    wGL = [ .3478548451; .6521451549; .6521451549; .3478548451 ];
    
elseif nGL == 3

    nodesGL(1,:) = midpt - .7745966692 * tangent/2;
    nodesGL(3,:) = midpt + .7745966692 * tangent/2;

    nodesGL(2,:) = midpt;

    wGL = [ 5/9; 8/9; 5/9 ];
    
elseif nGL == 2
    
    nodesGL(1,:) = midpt - .5773502692 * tangent/2;
    nodesGL(2,:) = midpt + .5773502692 * tangent/2;
    
    wGL = [ 1; 1 ];
    
elseif nGL == 1
    
    nodesGL(1,:) = midpt;
    
    wGL = 2;
    
end

xGL = nodesGL(:,1);  yGL = nodesGL(:,2);

wGL = wGL * norm(tangent)/2;