clear;
clc;

%% Parameters
% Grid size
Nx = 128; 
Ny = 64;
Lx = 2.0;
Ly = 1.0;
% Grid spacing
dx = Lx / Nx;
dy = Ly / Ny;
el_fx = [Nx, (1 : (Nx - 1))]; % forward index for x
el_bx = [(2 : Nx), 1]; % backward index for x
el_fy = [Ny, (1 : (Ny - 1))]; % same as above.
el_by = [(2 : Ny), 1];

% Fluid Parameters
rho = 1; % density
mu = 0.01; % Kinematic Viscosity

% Material Parameters
Num_b = 100; % Number of Lagrangian points on the boundary
K = 1; % Modulus of Elasticity 
index_forward = [Num_b, (1 : (Num_b - 1))];
index_backward = [(2 : Num_b), 1];
% dtheta = 2 * pi / Num_b; % Lagrange Coordinate distance
dtheta = 1 / Num_b;

% Time spacing
Tmax = 0.1;
dt = 0.001;
n_steps = floor(Tmax / dt);

%% Initial Conditions

% Initial Boundary
y_ini = 0.1; % define the displacement from center line for the tube. This is the radius of the tube at initial time
X_coord = linspace(0, Lx, Num_b);

% Top and Bottom wall
X_top = X_coord;
X_bot = X_coord;
Y_top = ones(1, Num_b) * (Ly / 2 + y_ini);
Y_bot = ones(1, Num_b) * (Ly / 2 - y_ini);

p_top(:, 1) = X_top;
p_top(:, 2) = Y_top;
p_bot(:, 1) = X_bot;
p_bot(:, 2) = Y_bot;
% % Now the 1st, Num_b -th, Num_b + 1 - th and 2 * Num_b -th points should be
% % fixed. 
fixed_ind = [1, Num_b];
Top_fix(:, 1) = X_top(fixed_ind);
Top_fix(:, 2) = Y_top(fixed_ind);
Bot_fix(:, 1) = X_bot(fixed_ind);
Bot_fix(:, 2) = Y_bot(fixed_ind);
% Initialize Fluid Fields
u = zeros(Nx, Ny, 2); % 2 for two directions

% Grid coordinates for cell centers
x_grid = linspace(dx/2, Lx-dx/2, Nx);
y_grid = linspace(dy/2, Ly-dy/2, Ny);
[PX, PY] = meshgrid(x_grid, y_grid);

%% Initial Velocity 
% If greater than Ly / 2 - y_ini or smaller than Ly / 2 + y_ini the
% velocity in x direction is 1. 
for j = 1 : Ny
    y = j * dy;
    if y < (Ly / 2 + y_ini) && y > (Ly / 2 - y_ini)
        u(:, j, 1) = 1; % Set velocity in x direction to 1 for all x positions
    end
end

%% a matrix used in solver
init_a_tube

size(a)

%% Plot 
figure;
hold on
% --- Plot Velocity Field ---
% quiver needs (y, x) ordering from meshgrid, so we permute u
u_x_plot = permute(u(:,:,1), [2 1]);
u_y_plot = permute(u(:,:,2), [2 1]);

% Plot a subset of vectors so the plot is readable
plot_stride = 4; % Plot every 4th vector
quiver(PX(1:plot_stride:end, 1:plot_stride:end), ...
       PY(1:plot_stride:end, 1:plot_stride:end), ...
       u_x_plot(1:plot_stride:end, 1:plot_stride:end), ...
       u_y_plot(1:plot_stride:end, 1:plot_stride:end), ...
       0.5, 'k'); % 0.5 scales the arrows

% --- Plot Material Boundary ---
plot(X_top, Y_top, 'r-', 'LineWidth', 2, 'DisplayName', 'Top Wall');
plot(X_bot, Y_bot, 'b-', 'LineWidth', 2, 'DisplayName', 'Bottom Wall');

% --- Plot Fixed Points ---
% plot(X_fixed, Y_fixed, 'ko', 'MarkerFaceColor', 'black', 'MarkerSize', 8, 'DisplayName', 'Fixed Ends');

% --- Add Labels and Final Touches ---
hold off;
axis equal;
xlim([0 Lx]);
ylim([0 Ly]);
xlabel('X Position');
ylabel('Y Position');
title('Initial State');
legend('show', 'Location', 'northwest');
grid on;


%% Main Time Loop
fprintf('Starting simulation...\n');

for clock = 1 : n_steps
    % Half Step (Prelinary step)
    p_half_top = p_top + dt / 2 * interpolation(u, p_top, Num_b, Nx, Ny, dx);
    p_half_bot = p_bot + dt / 2 * interpolation(u, p_bot, Num_b, Nx, Ny, dx);

    % Obtain the force adding on the boundary
    Force_top_temp = force(p_half_top, K, index_forward, index_backward, dtheta);
    Force_bot_temp = force(p_half_bot, K, index_forward, index_backward, dtheta);

    % Spread the force to Euclidean Coordinate
    Force_top = spreadforce(Force_top_temp, p_top, dtheta, dx, Num_b, Nx, Ny);
    Force_bot = spreadforce(Force_bot_temp, p_bot, dtheta, dx, Num_b, Nx, Ny);

    % The total force is just the sum of these two forces
    Force = Force_top + Force_bot;

    % Update the Fluid Velocity
    [u_final, u_intermediate] = navier_stokes_solver(u, Force, a, el_fx, el_bx, el_fy, el_by, dx, dt, mu, rho);

    % Final step
    p_final_top = p_top + dt * interpolation(u_final, p_half_top, Num_b, Nx, Ny, dx);
    p_final_bot = p_bot + dt * interpolation(u_final, p_half_bot, Num_b, Nx, Ny, dx);

    % Update data
    u = u_final;
    p_top = p_final_top;
    p_bot = p_final_bot;

    % Reset the boundary points to fixed points
    p_top(fixed_ind, :) = Top_fix;
    p_bot(fixed_ind, :) = Bot_fix;

    %% Animation
    % Calculate vorticity
    u_x = u(:, :, 1);
    u_y = u(:, :, 2);

    vorticity = (u_x(:, el_fy) - u_x(:, el_by)) / (2 * dy) - (u_y(el_fx, :) - u_y(el_bx, :)) / (2 * dx);
    
    clf;

    contour(PX, PY, permute(vorticity, [2 1]), 20, 'LineStyle', 'none');

    % h_bar = colorbar;
    % ylabel(h_bar, 'Vorticity (1/s)');
    % clim([-5, 5]); % Lock color axis for better comparison

    hold on;
    % Plot boundary
    plot(p_top(:, 1), p_top(:, 2), 'ko', 'LineWidth', 2.5, 'DisplayName', 'Top Boundary');
    plot(p_bot(:, 1), p_bot(:, 2), 'ko', 'LineWidth', 2.5, 'DisplayName', 'Bottom Boundary');
    title(sprintf('Immersed Boundary Simulation | Time = %.3f s', clock * dt));

    % Others
    axis equal;
    xlabel(' X position');
    ylabel(' Y position');
    % Update the force and Fluid, Fluid calculate using N-S solver. 

    drawnow
    hold off

end

fprintf(' Simulation Complete.\n');