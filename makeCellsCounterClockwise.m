function vd = makeCellsCounterClockwise( vd )

for i = 1 : size(vd,1)
    vv = vd{i};
    mid = sum(vv) ./ size(vv,1);
    midSph = mid ./ norm(mid,2);
    vv = rotateSphereNodes( vv, midSph(1), midSph(2), midSph(3) );
    mid = sum(vv) ./ size(vv,1);
    v1 = vv(1,:) - mid;
    v2 = vv(2,:) - mid;
    k = cross( v1, v2 );
    if k < 0
        vd{i} = flipud( vd{i} );
    end
end