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

%% Setup Grid for Plotting
% Create meshgrid (Note: meshgrid returns Y as rows, X as columns)
x_coords = linspace(0, Lx, Nx);
y_coords = linspace(0, Ly, Ny);
[X_grid, Y_grid] = meshgrid(x_coords, y_coords);

%% Animation Settings
figure('Name', 'Velocity Field Animation', 'Color', 'w');
axis equal;
axis([0 Lx 0 Ly]);
box on;
set(gca, 'FontSize', 12);
xlabel('X (m)');
ylabel('Y (m)');

% 1. Animation Speed Control
stride = 10; % Plot every 'stride' time steps

% 2. Arrow Density Control
% Plotting every single point (128x64) makes the plot too crowded.
% We skip points to make the arrows visible.
decimation_factor = 4; 

% Create indices for downsampling
ix = 1:decimation_factor:Nx;
iy = 1:decimation_factor:Ny;

% Create downsampled grids for quiver
X_sub = X_grid(iy, ix);
Y_sub = Y_grid(iy, ix);

% Calculate scaling for arrows to keep them consistent
% You might need to adjust 'scale_factor' depending on your flow magnitude
scale_factor = 0.5; 

fprintf('Starting velocity animation...\n');

for k = 1 : stride : n_steps
    % 1. Extract velocity at the current time step
    u_current = u_history(:, :, :, k);
    u_vel = u_current(:, :, 1); 
    v_vel = u_current(:, :, 2); 

    % 2. Clear previous frame
    cla; 
    hold on;

    % 3. Plot Velocity Magnitude Background (Optional, helps visualization)
    % Calculate magnitude
    mag = sqrt(u_vel.^2 + v_vel.^2);
    % Transpose to match meshgrid (Nx, Ny) -> (Ny, Nx)
    % Using 'pcolor' for a smooth background
    h = pcolor(X_grid, Y_grid, mag'); 
    set(h, 'EdgeColor', 'none'); % Turn off grid lines
    colormap('parula');
    % cb = colorbar; % Uncomment if you want to see the scale bar
    % clim([0 2]);   % Fix color limits if the colors flash too much

    % 4. Plot Velocity Vectors (Quiver)
    % Transpose and downsample data to match X_sub, Y_sub
    u_sub = u_vel(ix, iy)'; 
    v_sub = v_vel(ix, iy)';
    
    quiver(X_sub, Y_sub, u_sub, v_sub, scale_factor, 'k', 'LineWidth', 1);

    % 5. Plot Material Boundary (Fixed Tube)
    plot(X_tube_top(:, 1), X_tube_top(:, 2), 'r-', 'LineWidth', 2);
    plot(X_tube_bot(:, 1), X_tube_bot(:, 2), 'r-', 'LineWidth', 2);
    
    % 6. Title and Formatting
    current_time = k * sample_rate * dt;
    title(sprintf('Vorticity Contour at t = %.3f s', current_time));
    axis([0 Lx 0 Ly]);
    
    hold off;
    
    % Force update
    drawnow;
end

fprintf('Animation finished.\n');