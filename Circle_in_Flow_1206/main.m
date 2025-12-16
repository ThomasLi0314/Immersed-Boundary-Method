initialize;

init_a_circ;

% Set up for storing
u_history = zeros(Nx, Ny, 2, n_saved);
X_history = zeros(Num_b, 2, n_saved);

%% Main Loop
for clock = 1 : n_steps
    %% Reset Initial COndition
    % Speed at left side
    u(1 : Nx / 32, :, 1) = u_left;
    u(1 : Nx / 32, :, 2) = 0;

    % Speed at Right side
    % u(end - 2 : end, :, 1) = 0;
    
    %% Premilary Steps
    % half step for X
    X_half = X + dt / 2 * interpolation(u, X, Num_b, Nx, Ny, dx, dy);

    % Calculate Force
    F_temp = K * (Y - X_half);

    % Spread force to Eucleadian Space
    f = spreadforce(F_temp, X_half, dtheta, dx, dy, Num_b, Nx, Ny);

    % Add Perturbation
    if abs(clock - per_step) < per_time_steps
        f = f + f_perturbed;
    end


    % Navier Stokes Solver
    [u_final, u_intermidiate] = navier_stokes_solver_pIB(u, f, a, ind_fx, ind_bx, ind_fy, ind_by, dx, dy, dt, mu, rho);

    %% Final Step
    X_final = X + dt * interpolation(u_final, X_half, Num_b, Nx, Ny, dx, dy);

    % Update Data
    u = u_final;
    X = X_final;

    % Store data
    if mod(clock, sample_rate) == 0
        save_idx = floor(clock / sample_rate); % Should be Integer. 
        u_history(:, :, :, save_idx) = u;
        X_history(:, :, save_idx) = X;
    end

    % Display Progress
    if mod(clock, 500) == 0
        fprintf('Step %d / %d completed. \n', clock, n_steps)
    end   
end

%% Save Data 
fprintf('Simulation complete. Saving data...\n');
output_folder = 'Simulation_Results';

% 1. Generate a filename-safe timestamp (No colons or spaces)
% Format: YearMonthDay_HourMinuteSecond (e.g., 20251210_134500)
timestamp = string(datetime('now'), 'yyyyMMdd_HHmmss'); 

% 2. Ensure output folder exists (prevents "Directory not found" errors)
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 3. Create filenames
% Use char(timestamp) to ensure compatibility with sprintf if using older MATLAB versions
filename_sim = sprintf('simulation_result_u_rho%g_K%g_mu%g_t%g_%s.mat', rho, K, mu, Tmax, timestamp);
full_path_sim = fullfile(output_folder, filename_sim);

% 4. Save data
save(full_path_sim, 'u_history', 'X_history', 'Y', 'Nx', 'Ny', 'Lx', 'Ly', 'dx', 'dy', 'n_steps', 'dt', 'sample_rate');

fprintf('Data saved to: %s\n', full_path_sim);