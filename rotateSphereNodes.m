%Assuming that xyz is a N-by-3 matrix giving coordinates on the unit
%sphere, this function returns a new N-by-3 matrix in which all points have
%been rotated so that (x0,y0,z0) ends up at the north pole, (0,0,1).

function [ xyz, S, T ] = rotateSphereNodes( xyz, x0, y0, z0 )

ep = 1e-12;

if norm( [x0,y0,z0]-[0,0,-1], 2 ) <= ep
    xyz(:,[2,3]) = -xyz(:,[2,3]);
elseif norm( [x0,y0,z0]-[0,0,1], 2 ) <= ep
else
    %get angle between (x0,y0,z0) and (0,0,1):
    th = acos( z0 ./ sqrt(x0^2+y0^2+z0^2) );
    %transform to new cartesian coordinate system:
    S = [ [x0;y0;z0], [-z0*x0;-z0*y0;1-z0^2]/sqrt(1-z0^2), [y0;-x0;0]./sqrt(1-z0^2) ];
    xyz = xyz * S;
    %rotate so that (x0,y0,z0) is at the north pole:
    T = [ [cos(th);-sin(th);0], [sin(th);cos(th);0], [0;0;1] ];
    xyz = xyz * T;
    %go back to original cartesian coordinates:
    xyz = xyz / S;
end