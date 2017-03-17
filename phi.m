function z = phi( h, x, y, rbforder, ~ )

z = 1/h^rbforder * ( x.^2 + y.^2 ) .^ (rbforder/2);