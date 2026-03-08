% Vectorized force spreading — processes all Num_b boundary nodes at once.
% Drop-in replacement for the loop-based spreadforce.
%
% Uses accumarray to scatter-accumulate the 4x4 stencil contributions
% from every boundary node onto the Eulerian grid in one batch.

function f_euler = spreadforce_vec(F, X, dtheta, dx, dy, Num_b, Nx, Ny)
    c = dtheta / (dx * dy);

    % Scaled positions (Num_b x 1)
    s_x = X(:, 1) / dx;
    s_y = X(:, 2) / dy;
    i_x = floor(s_x);
    i_y = floor(s_y);
    r_x = s_x - i_x;       % fractional part in x
    r_y = s_y - i_y;       % fractional part in y

    % --- 1D delta weights (Num_b x 4), same as interpolation_vec ---
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
    offsets = int32(-1:2);
    ix = mod(int32(i_x) + offsets, Nx) + 1;   % Num_b x 4
    iy = mod(int32(i_y) + offsets, Ny) + 1;   % Num_b x 4

    % --- Expand to 16 stencil points per node (Num_b x 16) ---
    ix_rep = repmat(ix, 1, 4);                 % each ix repeated 4 times
    iy_rep = repelem(iy, 1, 4);                % each iy element repeated 4 times

    % 2D delta weights: outer product phi_x ⊗ phi_y → Num_b x 16
    w_x_rep = repmat(phi_x, 1, 4);
    w_y_rep = repelem(phi_y, 1, 4);
    w = w_x_rep .* w_y_rep;                    % Num_b x 16

    % --- Scatter-accumulate onto the grid using accumarray ---
    % Linear index into the (Nx x Ny) grid
    lin_idx = double(ix_rep) + (double(iy_rep) - 1) * Nx;   % Num_b x 16

    % Weighted force contributions (Num_b x 16) for each component
    vals_x = c * F(:, 1) .* w;   % broadcast: F(:,1) is Num_b x 1, w is Num_b x 16
    vals_y = c * F(:, 2) .* w;

    f_x = accumarray(lin_idx(:), vals_x(:), [Nx * Ny, 1], @sum, 0);
    f_y = accumarray(lin_idx(:), vals_y(:), [Nx * Ny, 1], @sum, 0);

    f_euler = cat(3, reshape(f_x, Nx, Ny), reshape(f_y, Nx, Ny));
end
