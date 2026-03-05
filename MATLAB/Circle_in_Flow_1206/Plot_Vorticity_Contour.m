clear;
clc;
close all;

%% Load Data with File Selection
output_folder = 'Simulation_Results';
animation_folder = 'Animations';

if ~exist(animation_folder, 'dir')
    mkdir(animation_folder);
end

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

if ~exist('sample_rate', 'var')
    warning('variable sample_rate not found. Assuming 1.');
    sample_rate = 1;
end

%% Video Writer Setup (Real-Time)
[~, name, ~] = fileparts(file);
video_filename = fullfile(animation_folder, [name, '_Vorticity_Contours.mp4']);
v = VideoWriter(video_filename, 'MPEG-4');
stride = 1; 
dt_frame = stride * sample_rate * dt; 
v.FrameRate = 1 / dt_frame; 
% Cap framerate if dt is too small, otherwise video players might struggle
if v.FrameRate > 60
    v.FrameRate = 60;
end
open(v);
fprintf('Video writer initialized. Frame rate: %.2f fps.\n', v.FrameRate);

%% Setup Grid
x_coords = linspace(0, Lx, Nx);
y_coords = linspace(0, Ly, Ny);
[X_grid, Y_grid] = meshgrid(x_coords, y_coords);

% Indices for periodic boundaries (or standard central difference)
ind_fx = [Nx, (1 : (Nx - 1))]; 
ind_bx = [(2 : Nx), 1]; 
ind_fy = [Ny, (1 : (Ny - 1))];
ind_by = [(2 : Ny), 1];

%% Animation Settings
figure('Name', 'Vorticity Contour Animation', 'Color', 'w');
set(gcf, 'Position', [100, 100, 900, 500]); 

fprintf('Starting animation...\n');

for k = 1 : stride : n_steps
    % 1. Extract Data
    u_current = u_history(:, :, :, k);
    X_current = X_history(:, :, k);
    
    u_vel = u_current(:, :, 1); 
    v_vel = u_current(:, :, 2); 
    
    % 2. Compute Vorticity
    % Central difference method
    dv_dx = (v_vel(ind_bx, :) - v_vel(ind_fx, :)) / (2 * dx);
    du_dy = (u_vel(:, ind_by) - u_vel(:, ind_fy)) / (2 * dy);
    vorticity = dv_dx - du_dy;
    
    % Transpose to match meshgrid (Ny x Nx)
    vort_plot = vorticity';
    
    % 3. Clear frame
    cla; 
    hold on;
    axis equal;
    axis([0 Lx 0 Ly]);
    box on;
    
    % 4. Plot Vorticity Contours
    % '40' is the number of contour levels. Increase this for finer detail.
    % 'LineWidth' is kept thin to differentiate fine structures.
    [C, h] = contour(X_grid, Y_grid, vort_plot, 40, 'LineWidth', 1.2); 
    
    % 5. Optimize Visibility
    % Calculate dynamic range for the current frame to highlight structures
    max_vort = max(abs(vort_plot(:)));
    if max_vort == 0
        max_vort = 1; % Prevent error on empty field
    end
    
    % Use a symmetric color limit so 0 is in the middle
    clim([-max_vort, max_vort]);
    
    % Use 'jet' or a diverging map. 'jet' is often best for high contrast 
    % physics data, while 'redblue' (if installed) is better for +/- separation.
    colormap(jet); 
    colorbar;
    
    % 6. Plot Material Boundary
    % Plotted in thick black on top of the contours
    plot(X_current(:, 1), X_current(:, 2), 'k-', 'LineWidth', 2.5); 
    
    % 7. Title and Formatting
    current_time = k * sample_rate * dt;
    title(sprintf('Vorticity Contours (t = %.3f s) | Max |\\omega| = %.2f', ...
        current_time, max_vort), 'Interpreter', 'tex');
    xlabel('X (m)');
    ylabel('Y (m)');
    set(gca, 'FontSize', 12);
    
    hold off;
    
    % 8. Capture
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Animation finished. Video saved to: %s\n', video_filename);