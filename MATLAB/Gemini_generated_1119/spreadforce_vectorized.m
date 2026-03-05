function f_euler = spreadforce_vectorized(F, X, dtheta, dx, Nx, Ny)
    % Spreads Lagrangian forces to Eulerian grid
    
    Num_b = size(X, 1);
    constant = dtheta / (dx * dx);
    
    s = X / dx;
    i_floor = floor(s);
    r = s - i_floor;
    
    base_x = i_floor(:, 1);
    base_y = i_floor(:, 2);
    
    % Arrays to hold data for accumarray
    % We expect at most Num_b * 16 entries (4x4 stencil)
    total_entries = Num_b * 16;
    ii_list = zeros(total_entries, 1); % X indices
    jj_list = zeros(total_entries, 1); % Y indices
    val_x_list = zeros(total_entries, 1); % Force X values
    val_y_list = zeros(total_entries, 1); % Force Y values
    
    count = 0;
    
    for m = -1:2
        for n = -1:2
            idx_x = mod(base_x + m - 1, Nx) + 1;
            idx_y = mod(base_y + n - 1, Ny) + 1;
            
            w_x = phi_func(r(:,1), m); % Same helper as interpolation
            w_y = phi_func(r(:,2), n);
            weights = w_x .* w_y * constant;
            
            % Store indices and values
            start_idx = count + 1;
            end_idx = count + Num_b;
            
            ii_list(start_idx:end_idx) = idx_x;
            jj_list(start_idx:end_idx) = idx_y;
            
            val_x_list(start_idx:end_idx) = F(:,1) .* weights;
            val_y_list(start_idx:end_idx) = F(:,2) .* weights;
            
            count = end_idx;
        end
    end
    
    % Accumulate duplicates into the grid
    fx_grid = accumarray([ii_list, jj_list], val_x_list, [Nx, Ny]);
    fy_grid = accumarray([ii_list, jj_list], val_y_list, [Nx, Ny]);
    
    f_euler = cat(3, fx_grid, fy_grid);
end

function val = phi_func(r, offset)
    % Same phi helper as interpolation
    root = sqrt(1 + 4*r - 4*r.^2);
    if offset == -1, val = (2 - 2*r - root) / 8;
    elseif offset == 0, val = (3 - 2*r + root) / 8;
    elseif offset == 1, val = (1 + 2*r + root) / 8;
    elseif offset == 2, val = (1 + 2*r - root) / 8;
    else, val = 0; end
end