clc,clear

n = 4;                      %number of bifurcations of icosahedral triangle
k = 1/10/2^(n-2);           %delta t
T = 6;                      %final time
t = 0 : k : T;              %vector of time values
rbforder    = 3;
polyorder   = 2;
stencilSize = 10;
nGL         = 3;
plotSolution   = 0;
computeWeights = 0;
computeParams  = 0;

%case 0 (solid body rotation, gaussian hill):
u_lam = @(lam,phi,t)  2*pi*cos(phi)/T;
u_phi = @(lam,phi,t)  zeros( size(lam) );
rho = @(xyz)  ones( size(xyz,1), 1 );
q = @(xyz)  exp( -5*( (xyz(:,1)+1/sqrt(2)).^2 + (xyz(:,2)+1/sqrt(2)).^2 + (xyz(:,3)-0).^2 ) );

% %case 1 (just deforming, nondivergent, cosine bells):
% K = 2.4;
% u_lam = @(lam,phi,t) K * sin(lam/2).^2 .* sin(2*phi) .* cos(pi*t/T);
% u_phi = @(lam,phi,t) K/2 * sin(lam) .* cos(phi) .* cos(pi*t/T);
% q = @(xyz)  exp( -5*( (xyz(:,1)+1/2).^2 + (xyz(:,2)-0).^2 + (xyz(:,3)-sqrt(3)/2).^2 ) ) + ...
%             exp( -5*( (xyz(:,1)+1/2).^2 + (xyz(:,2)-0).^2 + (xyz(:,3)+sqrt(3)/2).^2 ) );
% rho = @(xyz)  ones( size(xyz,1), 1 );

% %case 2 (just deforming, nondivergent, gaussian hills):
% K = 2;
% u_lam = @(lam,phi,t)  K * sin(lam).^2 .* sin(2*phi) .* cos(pi*t/T);
% u_phi = @(lam,phi,t)  K * sin(2*lam) .* cos(phi) .* cos(pi*t/T);
% q = @(xyz)  exp( -5*( (xyz(:,1)+sqrt(3)/2).^2 + (xyz(:,2)-1/2).^2 + (xyz(:,3)-0).^2 ) ) + ...
%             exp( -5*( (xyz(:,1)+sqrt(3)/2).^2 + (xyz(:,2)+1/2).^2 + (xyz(:,3)-0).^2 ) );
% rho = @(xyz)  ones( size(xyz,1), 1 );

% %case 3 (deformational, divergent, gaussian hills):
% K = 1;
% u_lam = @(lam,phi,t)  -K * sin(lam/2).^2 .* sin(2*phi) .* cos(phi).^2 .* cos(pi*t/T);
% u_phi = @(lam,phi,t)  K/2 * sin(lam) .* cos(phi).^3 .* cos(pi*t/T);
% q = @(xyz)  exp( -5*( (xyz(:,1)+sqrt(2)/2).^2 + (xyz(:,2)-sqrt(2)/2).^2 + (xyz(:,3)-0).^2 ) ) + ...
%             exp( -5*( (xyz(:,1)+sqrt(2)/2).^2 + (xyz(:,2)+sqrt(2)/2).^2 + (xyz(:,3)-0).^2 ) );
% rho = @(xyz)  ones( size(xyz,1), 1 );

% %case 4 (translating while it deforms, nondivergent, gaussian hills):
% K = 2;
% u_lam = @(lam,phi,t) K * sin(lam-2*pi*t/T).^2 .* sin(2*phi) .* cos(pi*t/T) + 2*pi*cos(phi)/T;
% u_phi = @(lam,phi,t) K * sin(2*(lam-2*pi*t/T)) .* cos(phi) .* cos(pi*t/T);
% q = @(xyz)  exp( -5*( (xyz(:,1)+sqrt(3)/2).^2 + (xyz(:,2)-1/2).^2 + (xyz(:,3)-0).^2 ) ) + ...
%             exp( -5*( (xyz(:,1)+sqrt(3)/2).^2 + (xyz(:,2)+1/2).^2 + (xyz(:,3)-0).^2 ) );
% rho = @(xyz)  ones( size(xyz,1), 1 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xyz = getIcosNodes( n, 0 );
if computeParams == 1
    vd = voronoiDiagramSph( xyz );
    vd = makeCellsCounterClockwise( vd );
    tic
    [ Wq, idx, U, V ] = sphIntegralWeights( xyz, vd );
    toc
    [ edges, vertices ] = getEdgesAndVertices( vd );
    save( ['n',num2str(n),'.mat'], 'vd', 'Wq', 'idx', 'U', 'V', 'edges', 'vertices' )
else
    load( ['n',num2str(n),'.mat'], 'vd', 'Wq', 'idx', 'U', 'V', 'edges', 'vertices' )
end
sphAreas2 = sum( Wq .* ( 1./(1+U.^2+V.^2).^(3/2) ), 2 );
F = rho(xyz) .* q(xyz);
avgVals = sum( Wq .* ( F(idx)./(1+U.^2+V.^2).^(3/2) ), 2 );
rhoq = avgVals ./ sphAreas2;
F = rho(xyz);
avgVals = sum( Wq .* ( F(idx)./(1+U.^2+V.^2).^(3/2) ), 2 );
rho = avgVals ./ sphAreas2;
exactAns = rhoq ./ rho;

indexGL = knnsearch( vertices, edges, 'k', 2 );
xv1 = vertices( indexGL(:,1), 1 );
yv1 = vertices( indexGL(:,1), 2 );
zv1 = vertices( indexGL(:,1), 3 );
xv2 = vertices( indexGL(:,2), 1 );
yv2 = vertices( indexGL(:,2), 2 );
zv2 = vertices( indexGL(:,2), 3 );

if computeWeights == 1
    tic
    [ W, xGL, yGL, zGL, weightsGL, U, V, idx, areas ] = makeInterpMatrixGL( ...
        xyz, edges, vd, ...
        rbforder, polyorder, stencilSize, nGL, ...
        xv1, yv1, zv1, xv2, yv2, zv2 );
    toc
    save( ['./weights/n',num2str(n),'r',num2str(rbforder), ...
        'p',num2str(polyorder),'stencil',num2str(stencilSize), ...
        'nGL',num2str(nGL), '.mat'], 'W', 'xGL', 'yGL', 'zGL', 'weightsGL', 'U', 'V', 'idx', 'areas' );
else
    load( ['./weights/n',num2str(n),'r',num2str(rbforder), ...
        'p',num2str(polyorder),'stencil',num2str(stencilSize), ...
        'nGL',num2str(nGL), '.mat'], 'W', 'xGL', 'yGL', 'zGL', 'weightsGL', 'U', 'V', 'idx', 'areas' );
end

sphAreas = sphAreas2(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ lam, phi, ~ ] = cart2sph( xGL, yGL, zGL );
uGL = @(t)  -u_lam(lam,phi,t) .* sin(lam) - u_phi(lam,phi,t) .* cos(lam) .* sin(phi);
vGL = @(t)   u_lam(lam,phi,t) .* cos(lam) - u_phi(lam,phi,t) .* sin(lam) .* sin(phi);
wGL = @(t)   u_phi(lam,phi,t) .* cos(phi);

cellIndex = knnsearch( xyz, edges, 'k', 2 );
cellIndexA = cellIndex(:,1);
cellIndexB = cellIndex(:,2);

normals = xyz( cellIndexB, : ) - xyz( cellIndexA, : );
normals = normals ./ repmat( sqrt(sum(normals.^2,2)), 1, 3 );
normals1 = repmat( normals(:,1), 1, nGL );
normals2 = repmat( normals(:,2), 1, nGL );
normals3 = repmat( normals(:,3), 1, nGL );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % for k = 1 : size(edges,1)
%     figure(1),clf
%     for j = 1 : size(vd,1)
%         vv = vd{j};
%         patch( vv(:,1), vv(:,2), vv(:,3), size(vv,1) );
%     end
%     hold( 'on' )
% %     plot3( edges(:,1), edges(:,2), edges(:,3), 'k*', 'markerSize', 5 )
%     plot3( xyz(:,1), xyz(:,2), xyz(:,3), 'ko' )
% %     plot3( xyz(cellIndex(k,:),1), xyz(cellIndex(k,:),2), xyz(cellIndex(k,:),3), 'ro' )
% %     plot3( edges(k,1), edges(k,2), edges(k,3), 'r*' )
% %     plot3( edges(k,1)+normals(k,1), edges(k,2)+normals(k,2), edges(k,3)+normals(k,3), 'r*' )
%     view(3), axis( 'square', 'off' ), colormap(copper);
%     title(sprintf('%g cells',size(xyz,1)))
%     hold( 'off' )
% %     pause
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tri = delaunaySph( xyz );
[ lam, phi, ~ ] = cart2sph( xyz(:,1), xyz(:,2), xyz(:,3) );
lam(lam<0) = lam(lam<0) + 2*pi;
triFlat = delaunay( lam, phi );

for i = 1 : length(t)-1
    %rho, stage 1:
    rhoS = rk3substep( rho, rho, t(i), k/3, uGL, vGL, wGL, normals1, normals2, normals3, ...
        cellIndexA, cellIndexB, idx, areas, sphAreas, sphAreas2, W, weightsGL, U, V );
    %rho, stage 2:
    rhoS = rk3substep( rho, rhoS, t(i)+k/3, k/2, uGL, vGL, wGL, normals1, normals2, normals3, ...
        cellIndexA, cellIndexB, idx, areas, sphAreas, sphAreas2, W, weightsGL, U, V );
    %rho, stage 3:
    rho = rk3substep( rho, rhoS, t(i)+k/2, k, uGL, vGL, wGL, normals1, normals2, normals3, ...
        cellIndexA, cellIndexB, idx, areas, sphAreas, sphAreas2, W, weightsGL, U, V );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %rhoq, stage 1:
    rhoS = rk3substep( rhoq, rhoq, t(i), k/3, uGL, vGL, wGL, normals1, normals2, normals3, ...
        cellIndexA, cellIndexB, idx, areas, sphAreas, sphAreas2, W, weightsGL, U, V );
    %rho, stage 2:
    rhoS = rk3substep( rhoq, rhoS, t(i)+k/3, k/2, uGL, vGL, wGL, normals1, normals2, normals3, ...
        cellIndexA, cellIndexB, idx, areas, sphAreas, sphAreas2, W, weightsGL, U, V );
    %rho, stage 3:
    rhoq = rk3substep( rhoq, rhoS, t(i)+k/2, k, uGL, vGL, wGL, normals1, normals2, normals3, ...
        cellIndexA, cellIndexB, idx, areas, sphAreas, sphAreas2, W, weightsGL, U, V );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tmp = abs( t(i+1) - (0:1:T) ) <= eps;
    if nnz(tmp) > 0
        if t(i+1) == T
            surf1 = rhoq./rho - exactAns;
            plotSolution = 1;
        else
            surf1 = rhoq ./ rho;
        end
        if plotSolution == 1
            figure(8),clf
            trisurf( tri, xyz(:,1), xyz(:,2), xyz(:,3), surf1 )
            shading( 'interp' )
            axis( 'tight' )
            daspect( [1,1,1] )
            colorbar
            title(sprintf('rhoq/rho, t=%1.2f, min=%1.2e, max=%1.2e, mass=%1.14f, tracerMass=%1.14f', ...
                t(i+1),min(surf1),max(surf1),sum(rho.*sphAreas2),sum(rhoq.*sphAreas2)))
            drawnow
            figure(9),clf
            trisurf( triFlat, lam, phi, surf1 )
            view(2)
            shading( 'interp' )
            axis( 'equal', 'tight' )
            colorbar
            xlabel( 'lambda' )
            ylabel( 'theta' )
            drawnow
        end
    end
end
