function w = laplacian_tube(u, el_fx, el_bx, el_fy, el_by, dx)
    w_temp = u(el_fx, :, :) + u(el_bx, :, :) + u(:, el_fy, :) + u(:, el_by, :) - 4 * u;
    w = w_temp / (dx ^2);
end