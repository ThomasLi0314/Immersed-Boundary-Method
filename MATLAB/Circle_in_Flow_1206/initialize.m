% Initialization


clear;
clc;
% close all;



%% Parameters and Setup

%% Parameters to be Changed
% Time spacing
Tmax = 12; % Total time of simulation (in second)
dt = 0.00001; % Time interval
n_steps = floor(Tmax / dt);
K = 100000; % Modulus of Elasticity

% Perturbation Time (add perturbation at t_per (s)
t_per = 1;
per_step = t_per / dt;

% Duration of the perturbation (This is the time the perturbation would
% last)
per_time = 0.2;
per_time_steps = per_time / dt;



% Grid size
Nx = 128;
Lx = 4.0;
Ly = 2.0;
Ny = Nx * Ly / Lx;

% Grid spacing
dx = Lx / Nx;
dy = Ly / Ny;

% Forward index (fetches i+1)
ind_fx = [(2 : Nx), 1];      

% Backward index (fetches i-1)
ind_bx = [Nx, (1 : (Nx - 1))]; 

% Same for Y
ind_fy = [(2 : Ny), 1];
ind_by = [Ny, (1 : (Ny - 1))];

% Fluid Parameters
rho = 1; % Density
mu = 0.01; % Kinematic Viscosity

% Material Parameters
Num_b = 100;
Tube_ind_f = [Num_b, (1 : Num_b - 1)];
Tubd_ind_b = [(2 : Num_b), 1];
dtheta = 2 * pi / Num_b; % Distance for lagrange coordinate

% Material Location
X = zeros(Num_b, 2);
L_geom = 1.0;
circ_center = [0.5, Ly /2];
radius = L_geom / 16;
for k = 1 : Num_b
    theta = k * dtheta;
    X(k, 1) = circ_center(1) + (radius) * cos(theta);
    X(k, 2) = circ_center(2) + (radius) * sin(theta);
end


% Ghost points
Y = X;

frame = 50; % (number of frame per second)

sample_rate = 1 / frame / dt;

% Initial Velocity
u = zeros(Nx, Ny, 2);

% Boundary velocity
u_left = 1.5;

u(1 : Nx / 32, :, 1) = u_left; % Velocity in the second and third column. 

% Initialized Storing
n_saved = floor(n_steps / sample_rate);

%% Force Perturbation (Comment this part if not needed
% Perturbation Location
f_perturbed = zeros(Nx, Ny);
per_ind_y = circ_center(2) / Ly * Ny;
per_ind_x_start = (circ_center(1) + radius) / Lx * Nx + 10;
per_ind_x_end = per_ind_x_start + 20;
per_force = 100.0;



% Add this perturbed matrix to the initial force
f_perturbed(per_ind_x_start : per_ind_x_end, per_ind_y) = per_force;
%% Plotting Initial Plot


% Vorticity and Voricity Contour
vorticity = (u(ind_fx, :, 2) - u(ind_bx, :, 2)) / (2 * dx) - (u(:, ind_fy, 1) - u(:, ind_by, 1)) / (2 * dy);

% Create Grid for plotting
[x_grid, y_grid] = meshgrid((0.5 : Nx - 0.5) * dx, (0.5 : Ny - 0.5) * dy);
% [xx, yy] = meshgrid(linspace(dx / 2, Lx - dx / 2, Nx), linspace(dy / 2 , Ly - dy / 2, Ny));

% 3. Generate Plot
figure(1);
clf;
hold on;

% Plot Vorticity Contour
% We transpose omega (omega') because MATLAB plots usually treat 
% the first dimension as Y (rows) and second as X (columns).
contourf(x_grid, y_grid, vorticity', 20, 'LineColor', 'none'); 
colormap(jet);
c = colorbar;
c.Label.String = 'Vorticity';

% Plot Material Boundary (The Tube)
% Connect the last point to the first to close the loop visually
plot([X(:, 1); X(1, 1)], [X(:, 2); X(1, 2)], 'k-', 'LineWidth', 2);
plot(X(:, 1), X(:, 2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 4);

% Formatting
axis equal;
axis([0 Lx 0 Ly]);
title('Vorticity Contour and Material Boundary');
xlabel('X');
ylabel('Y');
box on;
set(gca, 'Layer', 'top'); % Puts tick marks on top of the contour fill

hold off;
