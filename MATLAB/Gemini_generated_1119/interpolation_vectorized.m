function U = interpolation_vectorized(u, X, Nx, Ny, dx)
    % Vectorized interpolation from Eulerian grid to Lagrangian points
    % u: [Nx, Ny, 2]
    % X: [Num_b, 2]
    
    Num_b = size(X, 1);
    U = zeros(Num_b, 2);
    
    % Normalize coordinates
    s = X / dx;
    i_floor = floor(s);
    r = s - i_floor;
    
    % Base indices (1-based for Matlab)
    base_x = i_floor(:, 1); 
    base_y = i_floor(:, 2);
    
    % Loop over the 4x4 kernel stencil
    for m = -1:2
        for n = -1:2
            % Grid indices with periodic wrapping
            % Formula: mod(index - 1, N) + 1
            % We use round() to ensure they are treated as integers
            idx_x = round(mod(base_x + m - 1, Nx) + 1);
            idx_y = round(mod(base_y + n - 1, Ny) + 1);
            
            % Manual Linear Indexing (Replaces sub2ind)
            % Matlab stores matrices column-wise. 
            % Index = Row + (Col - 1) * NumRows
            lin_idx = idx_x + (idx_y - 1) * Nx;
            
            % Weights
            w_x = phi_func(r(:,1), m);
            w_y = phi_func(r(:,2), n);
            weights = w_x .* w_y;
            
            % Accumulate
            U(:, 1) = U(:, 1) + weights .* u(lin_idx); % u_x component
            
            % For u_y, we jump Nx*Ny elements to get to the 2nd page of the 3D array
            U(:, 2) = U(:, 2) + weights .* u(lin_idx + Nx*Ny); 
        end
    end
end

function val = phi_func(r, offset)
    % Standard Peskin 4-point function approximation
    root = sqrt(1 + 4*r - 4*r.^2);
    
    if offset == -1
        val = (2 - 2*r - root) / 8;
    elseif offset == 0
        val = (3 - 2*r + root) / 8;
    elseif offset == 1
        val = (1 + 2*r + root) / 8;
    elseif offset == 2
        val = (1 + 2*r - root) / 8;
    else
        val = 0;
    end
end