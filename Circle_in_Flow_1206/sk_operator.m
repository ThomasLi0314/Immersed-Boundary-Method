function S_u = sk_operator(u, g, el_fx, el_bx, el_fy, el_by, dx, dy)
    % u has size [Nx, Ny, 2] g has size [Nx, Ny] 

    % Recall the skew symmetric operator is defined by  
    % 0.5 u * Dg + 0.5 * D(u g). We calculate term by term

    % First part call I_11 and I_12 for two direction
    I_11 = u(:, :, 1) .* (g(el_fx, :) - g(el_bx, :));
    I_12 = u(:, :, 2) .* (g(:, el_fy) - g(:, el_by));

    % Second part call I_21 and I_22
    I_21 = u(el_fx, :, 1) .* g(el_fx, :) - u(el_bx, :, 1) .* g(el_bx, :);
    I_22 = u(:, el_fy, 2) .* g(:, el_fy) - u(:, el_by, 2) .* g(:, el_by);

    % Add up
    S_u = 0.5 * (I_11 + I_21) / (2 * dx) + 0.5 * (I_12 + I_22) / (2 * dy); % should consider dx and dy seperately, but here dx = dy

end
