% clear;
% clc;
% close all;
% %% Parameters
% % Grid size (High Resolution)
% Nx = 256; 
% Ny = 128; % Maintain 2:1 aspect ratio to match geometry
% Lx = 2.0;
% Ly = 1.0;
% 
% % Initial Velocity
% u_0 = 3.0;
% 
% % Grid spacing
% dx = Lx / Nx;
% dy = Ly / Ny;
% 
% % --- Finite Difference Indices ---
% % Forward index corresponds to (i+1), Backward to (i-1)
% el_fx = [(2 : Nx), 1];        % Forward index x
% el_bx = [Nx, (1 : (Nx - 1))]; % Backward index x
% el_fy = [(2 : Ny), 1];        % Forward index y
% el_by = [Ny, (1 : (Ny - 1))]; % Backward index y
% 
% % Fluid Parameters
% rho = 1; 
% mu = 0.02; % Viscosity
% 
% % Material Parameters
% % Num_b increased to prevent leakage at high resolution
% Num_b = 100; 
% K = 1.0;   % Stiffness (Elastic Modulus)
% index_forward = [Num_b, (1 : (Num_b - 1))];
% index_backward = [(2 : Num_b), 1];
% dtheta = Lx / Num_b; 
% 
% % Time spacing
% Tmax = 5.0; 
% dt = 0.0001; % Reduced time step for stability (CFL condition)
% n_steps = floor(Tmax / dt);
% 
% %% Initial Conditions
% 
% % Initial Boundary Geometry
% % We use a straight tube initially
% X_coord = linspace(0, Lx, Num_b)'; % Column vector
% y_ini = 0.15; 
% 
% % Top and Bottom wall initial positions
% p_top = zeros(Num_b, 2);
% p_bot = zeros(Num_b, 2);
% 
% p_top(:, 1) = X_coord;
% p_top(:, 2) = (Ly / 2 + y_ini); 
% 
% p_bot(:, 1) = X_coord;
% p_bot(:, 2) = (Ly / 2 - y_ini);
% 
% % --- FIXED POINTS SETUP ---
% % We identify the first and last points as the fixed anchors
% fixed_indices = [1, Num_b]; 
% 
% % Store the initial positions of these fixed points
% Top_fix = p_top(fixed_indices, :);
% Bot_fix = p_bot(fixed_indices, :);
% 
% % Initialize Fluid Fields
% u = zeros(Nx, Ny, 2); 
% 
% % Grid coordinates for plotting
% x_grid = linspace(dx/2, Lx-dx/2, Nx);
% y_grid = linspace(dy/2, Ly-dy/2, Ny);
% [PX, PY] = meshgrid(x_grid, y_grid);
% 
% % Initialize Velocity (Poiseuille-like profile or Plug flow inside)
% for j = 1 : Ny
%     y = j * dy;
%     if y < (Ly / 2 + y_ini) && y > (Ly / 2 - y_ini)
%         u(:, j, 1) = u_0; % Flow inside the tube
%     end
% end
% 
% %% Initialize Solver Operators
% init_a_tube 
% 
% %% Main Time Loop
% fprintf('Starting simulation...\n');
% figure('Units','normalized','Position',[0.2 0.2 0.5 0.6]);
% 
% gif_filename = 'simulation_result.gif';
% 
% for clock = 1 : n_steps
%     % Inside the loop, right after 'for clock = 1 : n_steps'
%     if any(isnan(u(:))) || any(isinf(u(:))) || max(abs(u(:))) > 100
%         error('Simulation Exploded at step %d! Max U: %f', clock, max(abs(u(:))));
%     end
%     %% 1. Half Step (Predictor)
%     % Interpolate velocity to boundary
%     u_b_top = interpolation_vectorized(u, p_top, Nx, Ny, dx);
%     u_b_bot = interpolation_vectorized(u, p_bot, Nx, Ny, dx);
% 
%     % Move boundary to half-step position
%     p_half_top = p_top + (dt / 2) * u_b_top;
%     p_half_bot = p_bot + (dt / 2) * u_b_bot;
% 
%     % --- ENFORCE FIXED BOUNDARY (Half Step) ---
%     p_half_top(fixed_indices, :) = Top_fix;
%     p_half_bot(fixed_indices, :) = Bot_fix;
% 
%     %% 2. Calculate Forces
%     % Calculate elastic force on the boundary
%     F_top_lag = force(p_half_top, K, index_forward, index_backward, dtheta);
%     F_bot_lag = force(p_half_bot, K, index_forward, index_backward, dtheta);
% 
%     % Spread force to Eulerian grid
%     F_top_eu = spreadforce_vectorized(F_top_lag, p_top, dtheta, dx, Nx, Ny);
%     F_bot_eu = spreadforce_vectorized(F_bot_lag, p_bot, dtheta, dx, Nx, Ny);
% 
%     Force = F_top_eu + F_bot_eu;
% 
%     %% 3. Fluid Solver (Navier-Stokes)
%     [u_final, u_intermediate] = navier_stokes_solver(u, Force, a, el_fx, el_bx, el_fy, el_by, dx, dt, mu, rho);
% 
%     %% 4. Final Step (Corrector)
%     % Interpolate new velocity to half-step boundary position
%     u_b_top_final = interpolation_vectorized(u_final, p_half_top, Nx, Ny, dx);
%     u_b_bot_final = interpolation_vectorized(u_final, p_half_bot, Nx, Ny, dx);
% 
%     % Update boundary position
%     p_top = p_top + dt * u_b_top_final;
%     p_bot = p_bot + dt * u_b_bot_final;
% 
%     % --- ENFORCE FIXED BOUNDARY (Final Step) ---
%     % Overwrite the positions of the endpoints with the fixed coordinates
%     p_top(fixed_indices, :) = Top_fix;
%     p_bot(fixed_indices, :) = Bot_fix;
% 
%     % Update fluid for next step
%     u = u_final;
% 
%     %% Animation
%     if mod(clock, 100) == 0
%         clf;
%         hold on;
% 
%         % Visualization: Velocity Magnitude
%         u_mag = sqrt(u(:,:,1).^2 + u(:,:,2).^2);
% 
%         % Pcolor for fluid speed
%         pcolor(PX, PY, permute(u_mag, [2 1])); 
%         shading interp;
%         colormap(jet);
%         cb = colorbar;
%         ylabel(cb, 'Velocity Magnitude');
%         clim([0 1.2]); 
% 
%         % Plot Vectors (subsampled)
%         stride = 2;
%         quiver(PX(1:stride:end, 1:stride:end), ...
%                PY(1:stride:end, 1:stride:end), ...
%                permute(u(1:stride:end, 1:stride:end, 1), [2 1]), ...
%                permute(u(1:stride:end, 1:stride:end, 2), [2 1]), ...
%                1.0, 'k');
% 
%         % Plot Material Boundaries
%         % Top Wall
%         plot(p_top(:, 1), p_top(:, 2), 'w-', 'LineWidth', 2);
%         plot(p_top(fixed_indices, 1), p_top(fixed_indices, 2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Fixed');
% 
%         % Bottom Wall
%         plot(p_bot(:, 1), p_bot(:, 2), 'w-', 'LineWidth', 2);
%         plot(p_bot(fixed_indices, 1), p_bot(fixed_indices, 2), 'ro', 'MarkerFaceColor', 'r');
% 
%         axis equal;
%         xlim([0 Lx]);
%         ylim([0 Ly]);
%         title(sprintf('Immersed Boundary | Time: %.3f', clock*dt));
%         drawnow;
% 
%         % --- GIF SAVING CODE ---
%         frame = getframe(gcf); 
%         im = frame2im(frame); 
%         [imind, cm] = rgb2ind(im, 256); 
% 
%         % Write to the GIF file 
%         if clock == 100 % Change this number to match your mod(clock, X)
%             % On the first frame, create the file and set LoopCount
%             imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1); 
%         else 
%             % On subsequent frames, append to the existing file
%             imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1); 
%         end
%     end
% end
% fprintf('Simulation Complete.\n');

%% Newest 1120
clear;
clc;

%% Check GPU
try
    g = gpuDevice;
    fprintf('Using GPU: %s\n', g.Name);
catch
    error('No compatible GPU found. Please run the CPU version.');
end

%% Parameters
Nx = 256; 
Ny = 128;
Lx = 2.0;
Ly = 1.0;
dx = Lx / Nx;
dy = Ly / Ny;

% --- GPU Array Creation ---
% Move indices to GPU immediately
el_fx = gpuArray([(2 : Nx), 1]);
el_bx = gpuArray([Nx, (1 : (Nx - 1))]);
el_fy = gpuArray([(2 : Ny), 1]);
el_by = gpuArray([Ny, (1 : (Ny - 1))]);

rho = 1; 
mu = 0.02; 

Num_b = 300; 
K = 0.1;   
index_forward = gpuArray([Num_b, (1 : (Num_b - 1))]);
index_backward = gpuArray([(2 : Num_b), 1]);
dtheta = Lx / Num_b; 

Tmax = 1.0; 
dt = 0.0001; 
n_steps = floor(Tmax / dt);

%% Initial Conditions (On GPU)

% Create coordinates directly on GPU
X_coord = gpuArray.linspace(0, Lx, Num_b)'; 
y_ini = 0.15; 

p_top = gpuArray.zeros(Num_b, 2);
p_bot = gpuArray.zeros(Num_b, 2);

p_top(:, 1) = X_coord;
p_top(:, 2) = (Ly / 2 + y_ini); 
p_bot(:, 1) = X_coord;
p_bot(:, 2) = (Ly / 2 - y_ini);

% Fixed indices
fixed_indices = gpuArray([1, Num_b]); 
Top_fix = p_top(fixed_indices, :);
Bot_fix = p_bot(fixed_indices, :);

% Initialize Fluid Fields on GPU
u = gpuArray.zeros(Nx, Ny, 2); 

% Coordinates for initialization
% We do this on CPU first then move to GPU to avoid complexity with meshgrid on GPU
x_grid = linspace(dx/2, Lx-dx/2, Nx);
y_grid = linspace(dy/2, Ly-dy/2, Ny);
[PX_cpu, PY_cpu] = meshgrid(x_grid, y_grid);
PX = gpuArray(PX_cpu);
PY = gpuArray(PY_cpu);

% Initialize velocity
for j = 1 : Ny
    y = j * dy;
    if y < (Ly / 2 + y_ini) && y > (Ly / 2 - y_ini)
        u(:, j, 1) = 1.0; 
    end
end

%% Initialize Solver Operators (Must be GPU adapted)
% We need to ensure 'init_a_tube' creates 'a' as a gpuArray.
% The easiest way is to run the script, then move 'a' to GPU.
init_a_tube; 
a = gpuArray(a); % Move the operator matrix to GPU

% GIF setup
currentDateTime = datetime('now');
timeStr = datestr(currentDateTime, 'yyyymmdd_HHMMSS'); %#ok<DATST>
fileExtension = '.gif';
gif_filename = ['simulation_gpu', timeStr, fileExtension];
if exist(gif_filename, 'file'), delete(gif_filename); end
first_frame = true;

%% Main Time Loop
fprintf('Starting GPU simulation...\n');

for clock = 1 : n_steps
    
    % Operations here happen on GPU automatically because inputs are gpuArrays
    
    % 1. Half Step
    u_b_top = interpolation_vectorized(u, p_top, Nx, Ny, dx);
    u_b_bot = interpolation_vectorized(u, p_bot, Nx, Ny, dx);
    
    p_half_top = p_top + (dt / 2) * u_b_top;
    p_half_bot = p_bot + (dt / 2) * u_b_bot;
    
    p_half_top(fixed_indices, :) = Top_fix;
    p_half_bot(fixed_indices, :) = Bot_fix;

    % 2. Forces
    F_top_lag = force(p_half_top, K, index_forward, index_backward, dtheta);
    F_bot_lag = force(p_half_bot, K, index_forward, index_backward, dtheta);

    F_top_eu = spreadforce_vectorized(F_top_lag, p_top, dtheta, dx, Nx, Ny);
    F_bot_eu = spreadforce_vectorized(F_bot_lag, p_bot, dtheta, dx, Nx, Ny);
    
    Force = F_top_eu + F_bot_eu;

    % 3. Fluid Solver (FFT on GPU)
    [u_final, u_intermediate] = navier_stokes_solver(u, Force, a, el_fx, el_bx, el_fy, el_by, dx, dt, mu, rho);

    % 4. Final Step
    u_b_top_final = interpolation_vectorized(u_final, p_half_top, Nx, Ny, dx);
    u_b_bot_final = interpolation_vectorized(u_final, p_half_bot, Nx, Ny, dx);

    p_top = p_top + dt * u_b_top_final;
    p_bot = p_bot + dt * u_b_bot_final;

    p_top(fixed_indices, :) = Top_fix;
    p_bot(fixed_indices, :) = Bot_fix;

    u = u_final;

    %% Animation (Must GATHER to CPU)
    if mod(clock, 20) == 0
        % --- Transfer Data back to CPU for plotting ---
        u_cpu = gather(u);
        p_top_cpu = gather(p_top);
        p_bot_cpu = gather(p_bot);
        fixed_indices_cpu = gather(fixed_indices);
        
        clf; hold on;
        
        u_mag = sqrt(u_cpu(:,:,1).^2 + u_cpu(:,:,2).^2);
        pcolor(PX_cpu, PY_cpu, permute(u_mag, [2 1])); 
        shading interp; colormap(jet); clim([0 1.2]);
        
        stride = 2;
        quiver(PX_cpu(1:stride:end, 1:stride:end), ...
               PY_cpu(1:stride:end, 1:stride:end), ...
               permute(u_cpu(1:stride:end, 1:stride:end, 1), [2 1]), ...
               permute(u_cpu(1:stride:end, 1:stride:end, 2), [2 1]), ...
               1.0, 'k');

        plot(p_top_cpu(:, 1), p_top_cpu(:, 2), 'w-', 'LineWidth', 2);
        plot(p_top_cpu(fixed_indices_cpu, 1), p_top_cpu(fixed_indices_cpu, 2), 'ro', 'MarkerFaceColor', 'r');
        plot(p_bot_cpu(:, 1), p_bot_cpu(:, 2), 'w-', 'LineWidth', 2);
        plot(p_bot_cpu(fixed_indices_cpu, 1), p_bot_cpu(fixed_indices_cpu, 2), 'ro', 'MarkerFaceColor', 'r');

        axis equal; xlim([0 Lx]); ylim([0 Ly]);
        title_full = sprintf(['GPU Simulation | Time: %.3f s\n\n', ...
        'Parameters: \n', ...
        '  Density (rho): %.2f, Kinematic Viscosity (mu): %.2e\n', ...
        '  Elastic Modulus (K): %.2e, Time Step (dt): %.2e\n', ...
        '  Max Time (Tmax): %.1f s'], ...
        clock * dt, rho, mu, K, dt, Tmax);

        title(title_full);
        drawnow;
        
        % Save GIF (Same as before)
        frame = getframe(gcf); 
        im = frame2im(frame); 
        [imind, cm] = rgb2ind(im, 256); 
        if first_frame
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1); 
            first_frame = false;
        else 
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1); 
        end 
    end
end