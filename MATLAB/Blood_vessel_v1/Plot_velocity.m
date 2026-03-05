clear;
clc;
close all;

%% Load Data with File Selection
output_folder = 'Simulation_Results';
animation_folder = 'Animations'; % New folder for saving videos

% Create animation folder if it doesn't exist
if ~exist(animation_folder, 'dir')
    mkdir(animation_folder);
end

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

%% Video Writer Setup (Real-Time)
% Construct output filename based on input filename
[~, name, ~] = fileparts(file);
video_filename = fullfile(animation_folder, [name, '_animation.mp4']);

% Create VideoWriter object
v = VideoWriter(video_filename, 'MPEG-4');

% --- Real-Time Calculation ---
% To make 1s animation = 1s real life, the FrameRate must match the
% reciprocal of the time elapsed between plotted frames.
% Time per step = dt * sample_rate
% Time per plotted frame = stride * dt * sample_rate
stride = 1; % Plot every 'stride' time steps
dt_frame = stride * sample_rate * dt; 

% Set FrameRate (1 / time_per_frame)
% Note: If dt_frame is very small, FrameRate might be very high. 
% Standard video is 24-60 fps. If FrameRate > 60, consider increasing stride.
v.FrameRate = 1 / dt_frame; 
% 
open(v);
fprintf('Video writer initialized. Frame rate set to %.2f fps for real-time playback.\n', v.FrameRate);

%% Setup Grid for Plotting
% Create meshgrid (Note: meshgrid returns Y as rows, X as columns)
x_coords = linspace(0, Lx, Nx);
y_coords = linspace(0, Ly, Ny);
[X_grid, Y_grid] = meshgrid(x_coords, y_coords);

%% Animation Settings
figure('Name', 'Velocity Field Animation', 'Color', 'w');
% Set figure size to ensure consistent video resolution (optional but recommended)
set(gcf, 'Position', [100, 100, 800, 400]); 

axis equal;
axis([0 Lx 0 Ly]);
box on;
set(gca, 'FontSize', 12);
xlabel('X (m)');
ylabel('Y (m)');

% 2. Arrow Density Control
% Plotting every single point makes the plot too crowded.
decimation_factor = 1; 

% Create indices for downsampling
ix = 1:decimation_factor:Nx;
iy = 1:decimation_factor:Ny;

% Create downsampled grids for quiver
X_sub = X_grid(iy, ix);
Y_sub = Y_grid(iy, ix);

% Calculate scaling for arrows to keep them consistent
scale_factor = 1; 

fprintf('Starting velocity animation and recording...\n');

% 1. Determine the number of saved frames
n_frames = size(u_history, 4);

for k = 1 : stride : n_frames
    % 1. Extract velocity at the current time step
    u_current = u_history(:, :, :, k);
    X_current_top = X_history_top(:, :, k);
    X_current_bot = X_history_bot(:, :, k);
    u_vel = u_current(:, :, 1); 
    v_vel = u_current(:, :, 2); 

    % 2. Clear previous frame
    cla; 
    hold on;

    % 3. Plot Velocity Magnitude Background
    mag = sqrt(u_vel.^2 + v_vel.^2)'; 
    h = pcolor(X_grid, Y_grid, mag); 
    set(h, 'EdgeColor', 'none'); 
    colormap('parula');
    shading interp
    c = colorbar;
    c.Label.String = 'Velocity Magnitude';

    % 4. Plot Velocity Vectors (Quiver)
    % Transpose and downsample data to match X_sub, Y_sub
    u_sub = u_vel(ix, iy)'; 
    v_sub = v_vel(ix, iy)';
    
    quiver(X_sub, Y_sub, u_sub, v_sub, scale_factor, 'k', 'LineWidth', 1);

    % 5. Plot Material Boundary (Fixed Tube)
    plot(X_current_top(:, 1), X_current_top(:, 2), 'r.', 'MarkerSize', 10);
    plot(X_current_bot(:, 1), X_current_bot(:, 2), 'r.', 'MarkerSize', 10);

    
    % 6. Title and Formatting
    current_time = k * sample_rate * dt;
    
    % Provide default values in case they were not saved in the .mat file
    % if ~exist('mu', 'var'), mu = 0.001; end
    % if ~exist('rho', 'var'), rho = 1; end
    % if ~exist('per_force', 'var'), per_force = 0.4; end
    % if ~exist('G', 'var'), G = 4; end
    
    title(sprintf('Velocity at t = %.3f s | \\mu = %g | \\rho = %g | N_x\\timesN_y = %d\\times%d | F_{per} = %g | F_{drive} = %g', ...
          current_time, mu, rho, Nx, Ny, per_force, G));
    axis([0 Lx 0 Ly]);
    
    hold off;
    
    % 7. Capture and write frame
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('Animation finished. Video saved to: %s\n', video_filename);