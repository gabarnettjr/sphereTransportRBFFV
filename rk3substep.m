function rho = rk3substep( rho0, rho, t, k, uGL, vGL, wGL, normals1, normals2, normals3, ...
    cellIndexA, cellIndexB, idx, areas, sphAreas, sphAreas2, W, weightsGL, U, V )

uDotN = uGL(t) .* normals1 + ...
        vGL(t) .* normals2 + ...
        wGL(t) .* normals3;

%rbf interp to GL nodes:
F = zeros( size(uDotN) );
rho = sphAreas ./ areas .* rho(idx);
for j = 1 : size(uDotN,2)
    F(:,j) = sum( W(:,:,j).*rho, 2 );
end
% F = F .* ( 1 + U.^2 + V.^2 ) .^ (3/2);
% %multiply by u dot n:
% F = F .* uDotN;
% %line integral of flux on the spherical arc:
% F = F ./ ( 1 + U.^2 + V.^2 );
F = F .* uDotN .* ( 1 + U.^2 + V.^2 ) .^ (1/2);
F = sum( F .* weightsGL, 2 );

% % low order centered:
% F = ( rho(cellIndexA) + rho(cellIndexB) ) ./ 2 .* uDotN .* edgeLengths;

% %low order upwind (donor cell):
% F = rho(cellIndexA) .* uDotN .* edgeLengths;
% ii = uDotN < 0;
% F(ii) = rho(cellIndexB(ii)) .* uDotN(ii) .* edgeLengths(ii);

%divergence:
D = zeros( size(rho0) );
for j = 1 : length(F)
    D(cellIndexA(j)) = D(cellIndexA(j)) + F(j);
    D(cellIndexB(j)) = D(cellIndexB(j)) - F(j);
end
D = D ./ sphAreas2;

rho = rho0 - k * D;