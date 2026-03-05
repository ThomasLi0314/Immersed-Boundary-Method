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

% Cross-section video writer
video_filename_cross = fullfile(animation_folder, [name, '_cross_section.mp4']);
v_cross = VideoWriter(video_filename_cross, 'MPEG-4');
v_cross.FrameRate = v.FrameRate;
open(v_cross);

%% Setup Grid for Plotting
% Create meshgrid (Note: meshgrid returns Y as rows, X as columns)
x_coords = linspace(0, Lx, Nx);
y_coords = linspace(0, Ly, Ny);
[X_grid, Y_grid] = meshgrid(x_coords, y_coords);

%% Animation Settings — Pre-create all plot objects once
fig_main = figure('Name', 'Velocity Field Animation', 'Color', 'w', 'Visible', 'off');
set(fig_main, 'Position', [100, 100, 800, 400]);
ax_main = axes(fig_main);
hold(ax_main, 'on');
axis(ax_main, 'equal');
axis(ax_main, [0 Lx 0 Ly]);
box(ax_main, 'on');
set(ax_main, 'FontSize', 12);
xlabel(ax_main, 'X (m)');
ylabel(ax_main, 'Y (m)');
colormap(ax_main, 'parula');

% Arrow Density Control
decimation_factor = 2;
ix = 1:decimation_factor:Nx;
iy = 1:decimation_factor:Ny;
X_sub = X_grid(iy, ix);
Y_sub = Y_grid(iy, ix);
scale_factor = 1;

% --- Create plot objects with dummy data (will be updated in loop) ---
% Velocity magnitude background (pcolor)
mag0 = zeros(Ny, Nx);
h_pcolor = pcolor(ax_main, X_grid, Y_grid, mag0);
set(h_pcolor, 'EdgeColor', 'none');
shading(ax_main, 'interp');
c = colorbar(ax_main);
c.Label.String = 'Velocity Magnitude';

% Quiver arrows
u_sub0 = zeros(size(X_sub));
v_sub0 = zeros(size(X_sub));
h_quiver = quiver(ax_main, X_sub, Y_sub, u_sub0, v_sub0, scale_factor, 'k', 'LineWidth', 1);

% Boundary markers
h_top = plot(ax_main, NaN, NaN, 'r.', 'MarkerSize', 10);
h_bot = plot(ax_main, NaN, NaN, 'r.', 'MarkerSize', 10);

% Title handle
h_title_main = title(ax_main, '');

fprintf('Starting velocity animation and recording...\n');

%% Poiseuille Cross-Section Setup
y_top_wall = Y_top(1, 2);
y_bot_wall = Y_bot(1, 2);
H_channel  = y_top_wall - y_bot_wall;

u_theory = 2 * G ./ (2 * mu) .* max(y_coords - y_bot_wall, 0) .* max(y_top_wall - y_coords, 0);

U_mean_theory = G * H_channel^2 / (8 * mu);
Re_channel    = rho * U_mean_theory * H_channel / mu;
fprintf('Reynolds number Re = %.4f\n', Re_channel);

cross_x_idx = round(Nx / 2);

% --- Cross-section figure: pre-create plot objects ---
fig_cross = figure('Name', 'Poiseuille Velocity Profile — Cross-Section', 'Color', 'w', 'Visible', 'off');
set(fig_cross, 'Position', [950, 100, 420, 520], 'Resize', 'off');
ax_cross = axes(fig_cross);
hold(ax_cross, 'on');
grid(ax_cross, 'on');

h_sim   = plot(ax_cross, NaN(1,Ny), y_coords, 'b-',  'LineWidth', 2.5, 'DisplayName', 'Simulated u_x');
h_theo  = plot(ax_cross, u_theory,  y_coords, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Theory: G/(2\mu)(y-y_b)(y_t-y)');
yline(ax_cross, y_top_wall, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(ax_cross, y_bot_wall, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel(ax_cross, 'u_x  (m/s)');
ylabel(ax_cross, 'y  (m)');
ylim(ax_cross, [0 Ly]);
legend(ax_cross, 'Location', 'northeast');   % fixed location (faster than 'best')
h_title_cross = title(ax_cross, '');

% Pre-compute grid spacings for vorticity
dy = Ly / (Ny - 1);
dx = Lx / (Nx - 1);

% Force one initial render so getframe works on invisible figures
drawnow;

%% Main Rendering Loop — update data only, no cla/recreate
n_frames = size(u_history, 4);
frame_indices = 1 : stride : n_frames;
n_total = numel(frame_indices);

for ii = 1 : n_total
    k = frame_indices(ii);
    fprintf('\rRendering frame %d / %d ...', ii, n_total);

    % Extract velocity at the current time step
    u_vel = u_history(:, :, 1, k);
    v_vel = u_history(:, :, 2, k);

    % ---- Update main figure objects (no cla) ----
    mag = sqrt(u_vel.^2 + v_vel.^2)';
    set(h_pcolor, 'CData', mag);

    set(h_quiver, 'UData', u_vel(ix, iy)', 'VData', v_vel(ix, iy)');

    X_current_top = X_history_top(:, :, k);
    X_current_bot = X_history_bot(:, :, k);
    set(h_top, 'XData', X_current_top(:,1), 'YData', X_current_top(:,2));
    set(h_bot, 'XData', X_current_bot(:,1), 'YData', X_current_bot(:,2));

    current_time = k * sample_rate * dt;
    set(h_title_main, 'String', sprintf('Velocity at t = %.3f s | Re = %.2f | \\mu = %g | G = %g | \\rho = %g', ...
          current_time, Re_channel, mu, G, rho));

    % Capture and write frame
    drawnow;
    frame = getframe(fig_main);
    writeVideo(v, frame);

    % ---- Update cross-section objects (no cla) ----
    u_cross = u_vel(cross_x_idx, :);
    set(h_sim, 'XData', u_cross);

    % Vorticity at cross-section: omega = dv/dx - du/dy
    dvdx_cross  = gradient(v_vel(cross_x_idx, :), dx);
    dudy_cross  = gradient(u_vel(cross_x_idx, :), dy);
    omega_cross = dvdx_cross - dudy_cross;
    omega_peak  = max(abs(omega_cross));

    set(h_title_cross, 'String', sprintf('x = %.2f m,  t = %.3f s  |  Re = %.2f  |  |\\omega|_{max} = %.3g', ...
          x_coords(cross_x_idx), current_time, Re_channel, omega_peak));

    u_max_plot = max([max(u_theory), max(abs(u_cross)), 1e-9]) * 1.15;
    xlim(ax_cross, [-u_max_plot * 0.05, u_max_plot]);

    drawnow;
    frame_cross = getframe(fig_cross);
    writeVideo(v_cross, frame_cross);
end

fprintf('\n');   % finish the progress line
close(v);
close(v_cross);
fprintf('Animation finished. Videos saved to:\n  %s\n  %s\n', video_filename, video_filename_cross);