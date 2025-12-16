% Shouldn't -1 in line 16 19.
% Checked Dec 13th. 
function U = interpolation(u, X, Num_b, Nx, Ny, dx, dy)
    U = zeros(Num_b, 2);
    for k = 1 : Num_b
        s_x = X(k, 1) / dx;
        s_y = X(k, 2) / dy;
        i_x = floor(s_x);
        i_y = floor(s_y);
        % Fractal Difference
        r_x = s_x - i_x;
        r_y = s_y - i_y;

        % Mod handles 
        ix_temp = (i_x - 1) : (i_x + 2);
        ix = mod(ix_temp, Nx) + 1; % Error Found!!!!!! 251213

        iy_temp = (i_y - 1) : (i_y + 2);
        iy = mod(iy_temp, Ny) + 1; % Error Found!!!!!!! 251213

        % Discrete Delta Function
        w = twoD_delta_discrete_delta(r_x, r_y);
        % Interpolate x velocity
        U(k, 1) = sum(sum(w .* u(ix, iy, 1)));
        % Interpolate y velocity
        U(k, 2) = sum(sum(w .* u(ix, iy, 2)));
    end
end

% Interpolation writen in Vector form (Should be faster). 