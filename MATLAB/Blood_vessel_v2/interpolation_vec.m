% Vectorized interpolation — processes all Num_b boundary nodes at once.
% Drop-in replacement for interpolation(u, X, Num_b, Nx, Ny, dx, dy).

function U = interpolation_vec(u, X, Num_b, Nx, Ny, dx, dy)
    % Scaled positions (Num_b x 1)
    s_x = X(:, 1) / dx;
    s_y = X(:, 2) / dy;
    i_x = floor(s_x);
    i_y = floor(s_y);
    r_x = s_x - i_x;       % fractional part in x
    r_y = s_y - i_y;       % fractional part in y

    % --- 1D delta weights (Num_b x 4) ---
    root_x = sqrt(1 + 4*r_x - 4*r_x.^2);
    root_y = sqrt(1 + 4*r_y - 4*r_y.^2);

    phi_x = zeros(Num_b, 4);
    phi_x(:, 1) = (3 - 2*r_x - root_x) / 8;
    phi_x(:, 2) = (3 - 2*r_x + root_x) / 8;
    phi_x(:, 3) = (1 + 2*r_x + root_x) / 8;
    phi_x(:, 4) = (1 + 2*r_x - root_x) / 8;

    phi_y = zeros(Num_b, 4);
    phi_y(:, 1) = (3 - 2*r_y - root_y) / 8;
    phi_y(:, 2) = (3 - 2*r_y + root_y) / 8;
    phi_y(:, 3) = (1 + 2*r_y + root_y) / 8;
    phi_y(:, 4) = (1 + 2*r_y - root_y) / 8;

    % --- Grid indices for the 4x4 stencil (Num_b x 4) ---
    offsets = int32(-1:2);  % stencil offsets: -1, 0, +1, +2
    ix = mod(int32(i_x) + offsets, Nx) + 1;   % Num_b x 4
    iy = mod(int32(i_y) + offsets, Ny) + 1;   % Num_b x 4

    % --- Gather velocity values at all 16 stencil points for all nodes ---
    %   For each node k, we need u(ix(k,:), iy(k,:), :) — a 4x4x2 block.
    %   We flatten this into linear indexing for speed.

    % Build linear index arrays: 16 indices per node (Num_b x 16)
    % Grid layout: u is (Nx, Ny, 2), linear index = ix + (iy-1)*Nx
    ix_rep = repmat(ix, 1, 4);                     % Num_b x 16: each ix repeated 4 times
    iy_rep = repelem(iy, 1, 4);                     % Num_b x 16: each iy element repeated 4 times
    lin_idx = double(ix_rep) + (double(iy_rep) - 1) * Nx;   % Num_b x 16

    % 2D delta weights: outer product phi_x(k,:) ⊗ phi_y(k,:) → Num_b x 16
    w_x_rep = repmat(phi_x, 1, 4);                 % Num_b x 16
    w_y_rep = repelem(phi_y, 1, 4);                 % Num_b x 16
    w = w_x_rep .* w_y_rep;                         % Num_b x 16

    % Gather u_x and u_y at all stencil points
    u_x_flat = u(:, :, 1);   % Nx x Ny → accessed via lin_idx
    u_y_flat = u(:, :, 2);

    u_x_vals = u_x_flat(lin_idx);   % Num_b x 16
    u_y_vals = u_y_flat(lin_idx);   % Num_b x 16

    % Weighted sum over the 16 stencil points → Num_b x 1
    U = [sum(w .* u_x_vals, 2), sum(w .* u_y_vals, 2)];
end
