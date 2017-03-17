function [ edges, vertices ] = getEdgesAndVertices( vd )

edges = [];
vertices = [];
for i = 1 : size(vd,1)
    vv = vd{i};
    tmp = ( vv + vv([2:end,1],:) ) ./ 2;
    tmp = tmp ./ repmat( sqrt(tmp(:,1).^2+tmp(:,2).^2+tmp(:,3).^2), 1, 3 );
    edges = [ edges; tmp ];
    vertices = [ vertices; vv ];
end
[ edges, ii ] = uniquetol( edges, 1e-12, 'ByRows', true );
vertices = uniquetol( vertices, 1e-12, 'ByRows', true );