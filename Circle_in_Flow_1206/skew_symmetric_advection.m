function w = skew_symmetric_advection(u,  el_fx, el_bx, el_fy, el_by, dx, dy)
    % The operator is defined by 
    % S(u) g = 0.5 ( grad.u) g + 0.5 ( div(u * g))

    % Here u is a [Nx, Ny, 2] matrix of vector field with two directions. 
    w = zeros(size(u));

    % First direction (x)
    w(:, :, 1) = sk_operator(u, u(:, :, 1), el_fx, el_bx, el_fy, el_by, dx, dy);

    % Second direction (y)
    w(:, :, 2) = sk_operator(u, u(:, :, 2), el_fx, el_bx, el_fy, el_by, dx, dy);

    end