initialize;

init_a_circ;

%% Parallel Pool Setup
if use_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local');
        fprintf('Parallel pool started with %d workers.\n', gcp('nocreate').NumWorkers);
    else
        fprintf('Reusing existing parallel pool (%d workers).\n', pool.NumWorkers);
    end
else
    fprintf('Serial mode  (set use_parallel = true in initialize.m to enable parallel pool).\n');
end

%% Profiling Accumulators
time_ns    = 0;   % cumulative time in NS solver
time_interp = 0;  % cumulative time in interpolation calls
time_spread = 0;  % cumulative time in spreadforce calls
total_tic  = tic; % wall-clock start for the whole simulation

% Set up for storage
u_history = zeros(Nx, Ny, 2, n_saved);
X_history_top = zeros(Num_b, 2, n_saved);
X_history_bot = X_history_top;

%% Display Reynolds number of this setting
% Wall positions from saved ghost/target points
y_top_wall = Y_top(1, 2);
y_bot_wall = Y_bot(1, 2);
H_channel  = y_top_wall - y_bot_wall;

U_mean_theory = G * H_channel^2 / (8 * mu);
Re_channel    = rho * U_mean_theory * H_channel / mu;

disp(Re_channel);

%% Main Loop
for clock = 1 : n_steps

    %% Premilary Steps
    t0 = tic;
    X_half_top = X_top + dt / 2 * interpolation(u, X_top, Num_b, Nx, Ny, dx, dy);
    X_half_bot = X_bot + dt / 2 * interpolation(u, X_bot, Num_b, Nx, Ny, dx, dy);
    time_interp = time_interp + toc(t0);

    % Calculate the total force. This force contains two part
    % Part 1, ghost mass and real mass
    F_penalty_top = K_tar * (Y_top - X_half_top);
    F_penalty_bot = K_tar * (Y_bot - X_half_bot);
    % set force at end point and start point to be 0 ( is this necessary?)
    % For F_penalty_top
    % F_penalty_top([1, size(F_penalty_top, 1)], :) = 0;
    % 
    % % For F_penalty_bot
    % F_penalty_bot([1, size(F_penalty_bot, 1)], :) = 0;

    % Part 2, force between adjoint points
    F_spring_top = zeros(Num_b, 2);
    F_spring_bot = zeros(Num_b, 2);

    % Vectorized spring force calculation
    i_vec = 2:Num_b-1;
    F_spring_top(i_vec, :) = K_mem * (2 * X_top(i_vec, :) - X_top(i_vec - 1, :) - X_top(i_vec + 1, :));
    F_spring_bot(i_vec, :) = K_mem * (2 * X_bot(i_vec, :) - X_bot(i_vec - 1, :) - X_bot(i_vec + 1, :));
    
    % Manually set the stiffness for connecting the first and second points, 
    % and the last and second last points very large. 
    % K_large_spring = K_mem * 1000; % You can adjust this large stiffness value
    K_large_spring = K_mem * 1;
    
    % Recalculate forces for the 2nd point (connecting 1st and 2nd)
    F_spring_top(2, :) = K_large_spring * (X_top(2, :) - X_top(1, :)) + K_mem * (X_top(2, :) - X_top(3, :));
    F_spring_bot(2, :) = K_large_spring * (X_bot(2, :) - X_bot(1, :)) + K_mem * (X_bot(2, :) - X_bot(3, :));
    
    % Recalculate forces for the second to last point (connecting last and second last)
    F_spring_top(Num_b-1, :) = K_mem * (X_top(Num_b-1, :) - X_top(Num_b-2, :)) + K_large_spring * (X_top(Num_b-1, :) - X_top(Num_b, :));
    F_spring_bot(Num_b-1, :) = K_mem * (X_bot(Num_b-1, :) - X_bot(Num_b-2, :)) + K_large_spring * (X_bot(Num_b-1, :) - X_bot(Num_b, :));

    F_tot_top = F_penalty_top + F_spring_top;
    F_tot_bot = F_penalty_bot + F_spring_bot;

    % Spread the force to Eucledian space
    t0 = tic;
    f_top = spreadforce(F_tot_top, X_half_top, dtheta, dx, dy, Num_b, Nx, Ny);
    f_bot = spreadforce(F_tot_bot, X_half_bot, dtheta, dx, dy, Num_b, Nx, Ny);
    time_spread = time_spread + toc(t0);

    f_boundary = f_top + f_bot; % Force from boundary to fluid;
    % Total force
    f = f_boundary + f_drive;

    % Add perturbation force to total force
    if abs(clock - per_step) < per_time_steps / 2
       f = f + f_perturbation;
    end

    % Navier Stokes Solver
    t0 = tic;
    [u_final, u_intermidiate] = navier_stokes_solver_pIB(u, f, a, ind_fx, ind_bx, ind_fy, ind_by, dx, dy, dt, mu, rho);
    time_ns = time_ns + toc(t0);

    %% Final step
    t0 = tic;
    X_final_top = X_top + dt * interpolation(u_final, X_half_top, Num_b, Nx, Ny, dx, dy);
    X_final_bot = X_bot + dt * interpolation(u_final, X_half_bot, Num_b, Nx, Ny, dx, dy);
    time_interp = time_interp + toc(t0);

    % Update the data
    u = u_final;
    X_top = X_final_top;
    X_bot = X_final_bot;

    % At final step, fix the start and end boundary points and set the end
    % velocity to be zero. 
    X_top(fix_ind, :) = Y_top(fix_ind, :);
    X_bot(fix_ind, :) = Y_bot(fix_ind, :);

    % 不能把尾端拽回来，把spring force做大
    % Note: Removed u(Nx-2:Nx,:,:) = 0 — the FFT NS solver assumes periodic BC;
    % zeroing outlet columns every step breaks periodicity and corrupts the solution.

    % Store the data
    if mod(clock, sample_rate) == 0
        save_idx = floor(clock / sample_rate); % Should be Integer.
        u_history(:, :, :, save_idx) = u;
        X_history_top(:, :, save_idx) = X_top;
        X_history_bot(:, :, save_idx) = X_bot;
    end

    % Display Progress
    if mod(clock, 2000) == 0
        fprintf('Step %d / %d completed. \n', clock, n_steps)
    end
end

%% Timing Report
total_wall = toc(total_tic);
fprintf('\n========= Timing Report =========\n');
fprintf('Total wall time       : %8.2f s\n', total_wall);
fprintf('NS solver             : %8.2f s  (%5.1f%%)\n', time_ns,    100*time_ns/total_wall);
fprintf('Spread force          : %8.2f s  (%5.1f%%)\n', time_spread, 100*time_spread/total_wall);
fprintf('Interpolation         : %8.2f s  (%5.1f%%)\n', time_interp, 100*time_interp/total_wall);
other_time = total_wall - time_ns - time_spread - time_interp;
fprintf('Other (spring, I/O..) : %8.2f s  (%5.1f%%)\n', other_time,  100*other_time/total_wall);
fprintf('=================================\n\n');

%% Save data
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
filename_sim = sprintf('simulation_result_u_rho%g_Ktar%g_Kmem%g_mu%g_per%g_G%g_t%g_%s.mat', rho, K_tar, K_mem, mu, per_force, G, Tmax, timestamp);
full_path_sim = fullfile(output_folder, filename_sim);

% 4. Save data
save(full_path_sim, 'u_history', 'X_history_top','X_history_bot', 'Y_top','Y_bot', 'Nx', 'Ny', 'Lx', 'Ly', 'dx', 'dy', 'n_steps', 'dt', 'sample_rate', 'per_force', 'G', 'f_drive', 'mu', 'rho');

fprintf('Data saved to: %s\n', full_path_sim);