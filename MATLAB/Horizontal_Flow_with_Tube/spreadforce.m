function f_euler = spreadforce(F, X, dtheta, dx, Num_b, Nx, Ny)
    % This function spreads the force from lagrangian coordinate to the
    % grid point for calculation in the Naview Stokes Equation
    constant = dtheta / (dx * dx);

    % Initilize the Force
    size_grid = size(X);
    f_euler = zeros(Nx, Ny, 2); % 2 for two directions
    for k = 1 : Num_b
        s = X(k, :) / dx; % i-th point on the boundary
        i = floor(s);
        r = s - i;

       % Mod handles periodic boundary condition. 
       i1_temp = (i(1) - 1) : (i(1) + 2);
       i1 = mod(i1_temp - 1, Nx) + 1; % Ensure indices are within bounds 

       i2_temp = (i(2) - 1) : (i(2) + 2);
       i2 = mod(i2_temp - 1, Ny) + 1;

       % Discrete Delta Functions
       w = tube_phi1(r(1)).*tube_phi2(r(2));

       % Spread in the x direction
       f_euler(i1, i2, 1) = (constant * F(k, 1)) * w;
       f_euler(i1, i2, 2) = (constant * F(k, 2)) * w;
    end
end



