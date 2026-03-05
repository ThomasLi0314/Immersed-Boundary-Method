% Spread the force, similar logic


function f_euler = spreadforce(F, X, dtheta, dx, dy, Num_b, Nx, Ny)
    % factor
    c = dtheta / (dx * dy);
    % Initiate force
    f_euler = zeros(Nx, Ny, 2);

    for k = 1 : Num_b
        s_x = X(k, 1) / dx;
        s_y = X(k, 2) / dy;
        % Floor
        i_x = floor(s_x);
        i_y = floor(s_y);
        % Fractal difference
        r_x = s_x - i_x;
        r_y = s_y - i_y;

        % Indices
        ix_temp = (i_x - 1) : (i_x + 2);
        ix = mod(ix_temp - 1, Nx) + 1;

        iy_temp = (i_y - 1) : (i_y + 2);
        iy = mod(iy_temp - 1, Ny) + 1;

        % Discrete Delta Function
        w = twoD_delta_discrete_delta(r_x, r_y);
        f_euler(ix, iy, 1) = f_euler(ix, iy, 1) + c * F(k, 1) * w; % F(k, 1) is a scalar, it just got spread to grid points
        f_euler(ix, iy, 2) = f_euler(ix, iy, 2) + c * F(k, 2) * w;
    end
end