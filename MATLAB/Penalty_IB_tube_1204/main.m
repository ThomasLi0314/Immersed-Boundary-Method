% The set up is as follow: 

% The tube is a rigid body, we use the coordinate X to denote the location
% of the tube( doesn't depend on time obviously). 

% There is another set of Y to denote a ghost mass connected to the system
% by a spring. Here Y depends on time. 
% 
% The force is spread based on location X, however the force
% add to the fluid at each point of the tube is propotional to the distance
% from X to Y. 

% We run through all the animation first, store the velocity at each time
% in each step and output a mat file storing all the data. Then we have
% another function (writen with GPT) that plot the vorticity contour. 


clear;
clc;
close all;

%% Parameters and Setup
% Grid size
Nx = 256;
Ny = 128;
Lx = 2.0;
Ly = 1.0;

% Grid spacing
dx = Lx / Nx;
dy = Ly / Ny;
ind_fx = [Nx, (1 : (Nx - 1))]; % forward index for x
ind_bx = [(2 : Nx), 1]; % backward index for x
ind_fy = [Ny, (1 : (Ny - 1))]; % same as above.
ind_by = [(2 : Ny), 1];

% Fluid Parameters
rho = 1; % Density
mu = 0.01; % Kinematic Viscosity

% Material Parameters
Num_b = 100;
K = 1000; % Modulus of Elasticity
L_tube = Lx;
Tube_ind_f = [Num_b, (1 : Num_b - 1)];
Tubd_ind_b = [(2 : Num_b), 1];
dtheta = L_tube / Num_b; % Distance for lagrange coordinate

% Time spacing
Tmax = 2.5; % Total time of simulation
dt = 0.001; % Time interval
n_steps = floor(Tmax / dt);

%% Initial and Boundary Conditions

y_ini = 0.1; % placement from the center line to the top or bottom of tube
tube_top = Ly / 2 + y_ini;
tube_bot = Ly / 2 - y_ini;

% X coordinate of tube (fix as time)
% Top
X_tube_top(:, 1) = linspace(0, Lx, Num_b);
X_tube_top(:, 2) = tube_top;
% Bottom
X_tube_bot(:, 1) = linspace(0, Lx, Num_b);
X_tube_bot(:, 2) = tube_bot;

% Initial Y coordinate of tube
Y_tube_top = X_tube_top;
Y_tube_bot = X_tube_bot;

% Initial Fluid Fields
u = zeros(Nx, Ny, 2); % u(:, :, 1) = velocity in x direction

% Rate of saving 
sample_rate = 20;

% The initial velocity only exists at the left boundary. 
for j = 1 : Ny
    y = j * dy;
    if y < tube_top - 2 * dy && y > tube_bot + 2 * dy
        u(1 : 2, j , 1) = 1; % Velocity only for the first and second grid points
    end
end

%% matrix a in navier stokes solver
init_a_tube

size(a)

%% Plot
%% Plot
figure('Name', 'Initial Vorticity Field', 'Color', 'w');
hold on;

% dv/dx
v = u(:, :, 2);
dv_dx = (v(ind_bx, :) - v(ind_fx, :)) / (2 * dx);

% du/dy
u_vel = u(:, :, 1);
du_dy = (u_vel(:, ind_by) - u_vel(:, ind_fy)) / (2 * dy);

vorticity = dv_dx - du_dy;

% 2. Setup Grid for Plotting
% Create meshgrid for visualization (Transpose required because meshgrid is Y,X)
[X_grid, Y_grid] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));

% 3. Plot the Vorticity Contour
% We transpose vorticity' to match the meshgrid dimensions (Ny, Nx)
contourf(X_grid, Y_grid, vorticity', 20, 'LineColor', 'none'); 
colormap('jet'); 
colorbar;
clim([-5 5]); % Adjust color limits based on expected vorticity magnitude
title('Initial Vorticity Contour');
xlabel('X');
ylabel('Y');

% 4. Overlay the Tube Structure
% Plot Top Wall
plot(X_tube_top(:, 1), X_tube_top(:, 2), 'k-', 'LineWidth', 2);
% Plot Bottom Wall
plot(X_tube_bot(:, 1), X_tube_bot(:, 2), 'k-', 'LineWidth', 2);

% 5. Plot the Ghost Mass / Spring connection (Optional visualization)
% Plot the initial y_ini reference lines if helpful
yline(tube_top, 'k--', 'Alpha', 0.3);
yline(tube_bot, 'k--', 'Alpha', 0.3);

axis equal;
axis([0 Lx 0 Ly]);
box on;
hold off;

%% Data Storage Initialization
% Dimensions: [Nx, Ny, components(2), time_steps]
n_saved = floor(n_steps / sample_rate);
u_history = zeros(Nx, Ny, 2, n_saved);

%% Main Time Loop
fprintf('Starting Simulation .../n');

for clock = 1 : n_steps
    %% Premilary Step

    % Update the ghost mass using velocity
    X_half_top = X_tube_top + dt / 2 * interpolation(u, X_tube_top, Num_b, Nx, Ny, dx, dy);
    X_half_bot = X_tube_bot + dt / 2 * interpolation(u, X_tube_bot, Num_b, Nx, Ny, dx ,dy);
    
    % Calculate the Force
    F_top_temp = K * (Y_tube_top - X_half_top);
    F_bot_temp = K * (Y_tube_bot - X_half_bot);

    % Spread the force to Eucleadian Coordinate
    Force_top = spreadforce(F_top_temp, X_half_top, dtheta, dx, dy, Num_b, Nx, Ny);
    Force_bot = spreadforce(F_bot_temp, X_half_bot, dtheta, dx, dy, Num_b, Nx, Ny);

    % Total force
    Force = Force_top + Force_bot;

    % Update the Fluid Velocity
    [u_final, u_intermidiate] = navier_stokes_solver_pIB(u, Force, a, ind_fx, ind_bx, ind_fy, ind_by, dx, dy, dt, mu, rho); % Navier stokes solver. 

    %% Final Step

    X_final_top = X_tube_top + dt * interpolation(u_final, X_half_top, Num_b, Nx, Ny, dx, dy);
    X_final_bot = X_tube_bot + dt * interpolation(u_final, X_half_bot, Num_b, Nx, Ny, dx, dy);

    % Update data
    u = u_final;
    X_tube_top = X_final_top;
    X_tube_bot = X_final_bot;

    % Reset the input velocity,
    for j = 1 : Ny
        y = j * dy;
        if y < tube_top - 2 * dy && y > tube_bot + 2 * dy
            u(1 : 2, j , 1) = 1; % Velocity only for the first and second grid points
        end
    end

    % Store the data for plot
    % Here since the Position of the Tube is fixed. We just need one file
    % for storing the material boundary. Another file to store the Velocity field
    if mod(clock, sample_rate) == 0
        save_idx = clock / sample_rate;
        u_history(:, :, :, save_idx) = u;
    end

    % Display Progress
    if mod(clock, 100) == 0
        fprintf('Step %d / %d completed.\n', clock, n_steps)
    end
end

%% Save Data
fprintf('Simulation complete. Saving data...\n');

% Define output folder
output_folder = 'Simulation_Results';

% Create the folder if it does not exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Generate timestamped filename
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
filename = sprintf('simulation_results_rho%g_K%g_mu%g_%s.mat', rho, K, mu, timestamp);
full_path = fullfile(output_folder, filename);

% Save variables needed for plotting
save(full_path, 'u_history', 'X_tube_top', 'X_tube_bot', ...
     'Nx', 'Ny', 'Lx', 'Ly', 'dx', 'dy', 'n_steps', 'dt', 'sample_rate');

fprintf('Data saved to: %s\n', full_path);


