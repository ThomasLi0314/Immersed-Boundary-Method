function omega = vorticity_cal(u, ind_fx, ind_bx, ind_fy, ind_by, dx, dy)
    p_xv = (u(ind_fx, :, 2) - u(ind_bx, :, 2)) / (2 * dx);
    p_yu = (u(:, ind_fy, 1) - u(:, ind_by, 1)) / (2 * dy);

    omega = p_xv - p_yu;
end
