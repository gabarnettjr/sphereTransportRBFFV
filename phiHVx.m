function z = phiHVx( h, x, y, rbforder, K )

if K == 1
    z = rbforder * x .* (x.^2+y.^2) .^ ( (rbforder-2)/2 );
elseif K == 2
    z = (rbforder-2)*rbforder^2 * x .* (x.^2+y.^2) .^ ( (rbforder-4)/2 );
elseif K == 3
    z = (rbforder-4)*(rbforder-2)^2*rbforder^2 * x .* (x.^2 + y.^2) .^ ( (rbforder-6)/2 );
end

z = 1/h^rbforder * z;