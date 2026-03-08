initialize;

init_a_circ;
init_b;

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
p_history = zeros(Nx, Ny, n_saved);
X_history_top = zeros(Num_b, 2, n_saved);
X_history_bot = X_history_top;
F_history_top = zeros(Num_b, 2, n_saved);  % total Lagrangian force on top wall
F_history_bot = zeros(Num_b, 2, n_saved);  % total Lagrangian force on bot wall

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
    X_half_top = X_top + dt / 2 * interpolation_vec(u, X_top, Num_b, Nx, Ny, dx, dy);
    X_half_bot = X_bot + dt / 2 * interpolation_vec(u, X_bot, Num_b, Nx, Ny, dx, dy);
    time_interp = time_interp + toc(t0);

    %% Force part

    %% Penalty (tether) force — applied ONLY at anchor nodes (fix_ind)
    F_penalty_top = zeros(Num_b, 2);
    F_penalty_bot = zeros(Num_b, 2);

    % F_penalty_top(fix_ind, :) = K_tar.* (Y_top(fix_ind, :) - X_half_top(fix_ind, :));
    % F_penalty_bot(fix_ind, :) = K_tar.* (Y_bot(fix_ind, :) - X_half_bot(fix_ind, :));

    % Use the exponential decay version
    F_penalty_top(:, :) = K_tar_vec' .* (Y_top(:, :) - X_half_top(:, :));
    F_penalty_bot(:, :) = K_tar_vec' .* (Y_bot(:, :) - X_half_bot(:, :));


    %% Part 2, spring force between adjacent points (Hooke's law, open ends)
    %  Num_b points -> (Num_b - 1) segments
    %  Segment k connects node k to node k+1
    dX_top = X_top(2:end, :) - X_top(1:end-1, :);            
    dX_bot = X_bot(2:end, :) - X_bot(1:end-1, :);
    
    % Segment lengths
    L_top  = sqrt(sum(dX_top.^2, 2));       
    L_bot  = sqrt(sum(dX_bot.^2, 2));
    
    % Extension beyond rest length
    stretch_top = L_top - ds0;              
    stretch_bot = L_bot - ds0;
    
    % Unit direction vectors (pointing from k to k+1)
    dir_top = dX_top ./ L_top;             
    dir_bot = dX_bot ./ L_bot;

    % Segment force magnitude & direction (Num_b - 1 x 2)
    % +seg_f pulls node k toward k+1
    % -seg_f pulls node k+1 toward k
    seg_f_top = (K_mem * stretch_top) .* dir_top;
    seg_f_bot = (K_mem * stretch_bot) .* dir_bot;

    % Accumulate forces onto nodes
    F_spring_top = zeros(Num_b, 2);
    F_spring_bot = zeros(Num_b, 2);
    
    % Node k experiences force from segment k (pulling forward)
    F_spring_top(1:end-1, :) = F_spring_top(1:end-1, :) + seg_f_top;
    F_spring_bot(1:end-1, :) = F_spring_bot(1:end-1, :) + seg_f_bot;
    
    % Node k+1 experiences negative force from segment k (pulling backward)
    F_spring_top(2:end, :) = F_spring_top(2:end, :) - seg_f_top;
    F_spring_bot(2:end, :) = F_spring_bot(2:end, :) - seg_f_bot;

    F_tot_top = F_penalty_top + F_spring_top;
    F_tot_bot = F_penalty_bot + F_spring_bot;

    %% Part 3 Myogenic (EC) force — vertical spring between opposing wall nodes
    %  M(s,t) = K_M * D/D0  +  gamma_M * (dD/dt) / D0    (Arthurs Eq. 60)
    D_current = abs(X_top(:,2) - X_bot(:,2));          % current diameter at each node
    dDdt      = (D_current - D_prev) / dt;              % finite-difference ∂_t D

    M_magnitude = K_M * (D_current ./ D0) + gamma_M * (dDdt ./ D0);

    % Direction: from each wall toward the vessel centreline (τ_M)
    %   upper wall → pushed down,  lower wall → pushed up
    F_myo_top        = zeros(Num_b, 2);
    F_myo_bot        = zeros(Num_b, 2);
    F_myo_top(:, 2)  = -M_magnitude;   % push upper wall downward
    F_myo_bot(:, 2)  =  M_magnitude;   % push lower wall upward

    % Do NOT apply myogenic force at anchor nodes (they are tethered)
    F_myo_top(fix_ind, :) = 0;
    F_myo_bot(fix_ind, :) = 0;

    D_prev = D_current;   % save for next time-step

    F_tot_top = F_tot_top + F_myo_top;
    F_tot_bot = F_tot_bot + F_myo_bot;

    % %% Part 4, Horizontal restoring force (x-direction only)
    % %  Prevents downstream drift; nodes stay near initial x-positions
    % F_hor_top = zeros(Num_b, 2);
    % F_hor_bot = zeros(Num_b, 2);
    % F_hor_top(:, 1) = K_hor * (Y_top(:, 1) - X_top(:, 1));
    % F_hor_bot(:, 1) = K_hor * (Y_bot(:, 1) - X_bot(:, 1));

    % F_tot_top = F_tot_top + F_hor_top;
    % F_tot_bot = F_tot_bot + F_hor_bot;

    %% Spread the force to Eucledian space
    t0 = tic;
    f_top = spreadforce_vec(F_tot_top, X_half_top, dtheta, dx, dy, Num_b, Nx, Ny);
    f_bot = spreadforce_vec(F_tot_bot, X_half_bot, dtheta, dx, dy, Num_b, Nx, Ny);
    time_spread = time_spread + toc(t0);

    f_boundary = f_top + f_bot; % Force from boundary to fluid;
    % Total force
    f = f_boundary + f_drive;

    % % Add perturbation force to total force
    % if abs(clock - per_step) < per_time_steps / 2
    %    f = f + f_perturbation;
    % end

    % Navier Stokes Solver
    t0 = tic;
    [u_final, u_intermidiate, p] = navier_stokes_solver_pIB(u, f, a, b, ind_fx, ind_bx, ind_fy, ind_by, dx, dy, dt, mu, rho);
    time_ns = time_ns + toc(t0);

    %% Final step
    t0 = tic;
    X_final_top = X_top + dt * interpolation_vec(u_final, X_half_top, Num_b, Nx, Ny, dx, dy);
    X_final_bot = X_bot + dt * interpolation_vec(u_final, X_half_bot, Num_b, Nx, Ny, dx, dy);
    time_interp = time_interp + toc(t0);

    % Update the data
    u = u_final;
    X_top = X_final_top;
    X_bot = X_final_bot;

    % At final step, fix the start and end boundary points and set the end
    % % velocity to be zero. 
    % X_top(fix_ind, :) = Y_top(fix_ind, :);
    % X_bot(fix_ind, :) = Y_bot(fix_ind, :);

    % 不能把尾端拽回来，把spring force做大
    % Note: Removed u(Nx-2:Nx,:,:) = 0 — the FFT NS solver assumes periodic BC;
    % zeroing outlet columns every step breaks periodicity and corrupts the solution.

    % Store the data
    if mod(clock, sample_rate) == 0
        save_idx = floor(clock / sample_rate); % Should be Integer.
        u_history(:, :, :, save_idx) = u;
        p_history(:, :, save_idx) = p;
        X_history_top(:, :, save_idx) = X_top;
        X_history_bot(:, :, save_idx) = X_bot;
        F_history_top(:, :, save_idx) = F_tot_top;
        F_history_bot(:, :, save_idx) = F_tot_bot;
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
% Build output path: Simulation_Results/<project_name>/<date>
results_base = 'D:\Documents\College\Research\Immersed Boundary Methods\Simulation_Results';
[project_dir, project_name] = fileparts(fileparts(mfilename('fullpath')));
date_folder = string(datetime('now'), 'yyyy-MM-dd');
output_folder = fullfile(results_base, project_name, date_folder);

% 1. Generate a filename-safe timestamp (No colons or spaces)
% Format: YearMonthDay_HourMinuteSecond (e.g., 20251210_134500)
timestamp = string(datetime('now'), 'yyyyMMdd_HHmmss');

% 2. Ensure output folder exists (prevents "Directory not found" errors)
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 3. Create filenames
% Use char(timestamp) to ensure compatibility with sprintf if using older MATLAB versions
filename_sim = sprintf('simulation_result_u_rho%g_Ktar%g_Kmem%g_mu%g_G%g_t%g_%s.mat', rho, K_tar, K_mem, mu, G, Tmax, timestamp);
full_path_sim = fullfile(output_folder, filename_sim);

% 4. Save data
save(full_path_sim, 'u_history', 'p_history', 'X_history_top','X_history_bot', 'Y_top','Y_bot', 'Nx', 'Ny', 'Lx', 'Ly', 'dx', 'dy', 'n_steps', 'dt', 'sample_rate', 'G', 'f_drive', 'mu', 'rho', ...
    'K_M', 'gamma_M', 'D0', 'fix_ind', 'F_history_top', 'F_history_bot');

fprintf('Data saved to: %s\n', full_path_sim);


% thest