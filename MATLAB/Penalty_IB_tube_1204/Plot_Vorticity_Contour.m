clear;
clc;
close all;

%% Load Data with File Selection
output_folder = 'Simulation_Results';

% Check if the folder exists to set initial path for file dialog
if exist(output_folder, 'dir')
    initial_path = fullfile(pwd, output_folder);
else
    initial_path = pwd;
end

fprintf('Please select the simulation .mat file to animate...\n');
[file, path] = uigetfile(fullfile(initial_path, '*.mat'), 'Select Simulation Data');

if isequal(file, 0)
    disp('User canceled file selection.');
    return;
else
    full_file_path = fullfile(path, file);
    fprintf('Loading %s...\n', full_file_path);
    load(full_file_path);
end

%% Handle Legacy Data (Check for sample_rate)
if ~exist('sample_rate', 'var')
    warning('variable sample_rate not found in data. Assuming 1 (saved every step).');
    sample_rate = 1;
end

% Determine number of saved frames
n_saved_frames = size(u_history, 4);

%% Setup Grid for Plotting
% Create meshgrid (Note: meshgrid returns Y as rows, X as columns)
x_coords = linspace(0, Lx, Nx);
y_coords = linspace(0, Ly, Ny);
[X_grid, Y_grid] = meshgrid(x_coords, y_coords);

% Indices for finite difference
ind_fx = [Nx, (1 : (Nx - 1))]; 
ind_bx = [(2 : Nx), 1]; 
ind_fy = [Ny, (1 : (Ny - 1))];
ind_by = [(2 : Ny), 1];

%% Animation Settings
figure('Name', 'Vorticity Animation', 'Color', 'w');
axis equal;
axis([0 Lx 0 Ly]);
box on;
set(gca, 'FontSize', 12);
xlabel('X (m)');
ylabel('Y (m)');

% Define fixed contour levels for consistent animation
contour_levels = linspace(-5, 5, 20); 

% Animation loop
% Plot every 'stride' frames (relative to the SAVED frames)
stride = 1; 

fprintf('Starting animation...\n');

for k = 1 : stride : n_saved_frames
    % 1. Extract velocity at the current SAVED frame
    u_current = u_history(:, :, :, k);
    u_vel = u_current(:, :, 1); 
    v_vel = u_current(:, :, 2); 

    % 2. Compute Vorticity
    dv_dx = (v_vel(ind_bx, :) - v_vel(ind_fx, :)) / (2 * dx);
    du_dy = (u_vel(:, ind_by) - u_vel(:, ind_fy)) / (2 * dy);
    vorticity = dv_dx - du_dy;

    % 3. Clear previous frame content but keep axis settings
    cla; 
    hold on;

    % 4. Plot Vorticity Contours (Lines only, Black color)
    contour(X_grid, Y_grid, vorticity', contour_levels, 'LineColor', 'k'); 

    % 5. Plot Material Boundary (Fixed Tube)
    plot(X_tube_top(:, 1), X_tube_top(:, 2), 'k-', 'LineWidth', 2);
    plot(X_tube_bot(:, 1), X_tube_bot(:, 2), 'k-', 'LineWidth', 2);
    
    % 6. Title and Formatting
    % Calculate actual simulation time: frame_index * sample_rate * dt
    current_time = k * sample_rate * dt;
    title(sprintf('Vorticity Contour at t = %.3f s', current_time));
    axis([0 Lx 0 Ly]); 
    
    hold off;
    
    % Force update
    drawnow;
end

fprintf('Animation finished.\n');