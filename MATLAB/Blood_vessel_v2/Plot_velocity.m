clear;
clc;
close all;

%% Plot Switches
plot_velocity_profile = true;  % Set true to include u_x cross-section profile at mid-x

%% Load Data with File Selection
% Build paths: Simulation_Results/<project_name>/<date>
results_base = 'D:\Documents\College\Research\Immersed Boundary Methods\Simulation_Results';
[~, project_name] = fileparts(fileparts(mfilename('fullpath')));
date_folder = string(datetime('now'), 'yyyy-MM-dd');
output_folder = fullfile(results_base, project_name);
animation_folder = fullfile(results_base, project_name, date_folder);

% Create animation folder (and parents) if it doesn't exist
if ~exist(animation_folder, 'dir')
    mkdir(animation_folder);
end

% Set initial path for file dialog to the project results folder
if exist(output_folder, 'dir')
    initial_path = output_folder;
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
% video_filename_cross = fullfile(animation_folder, [name, '_cross_section.mp4']);
% v_cross = VideoWriter(video_filename_cross, 'MPEG-4');
% v_cross.FrameRate = v.FrameRate;
% open(v_cross);

%% Setup Grid for Plotting
% Create meshgrid (Note: meshgrid returns Y as rows, X as columns)
x_coords = linspace(0, Lx, Nx);
y_coords = linspace(0, Ly, Ny);
[X_grid, Y_grid] = meshgrid(x_coords, y_coords);

%% Animation Settings — Pre-create all plot objects once
fig_main = figure('Name', 'Velocity & Pressure Animation', 'Color', 'w', 'Visible', 'off');
if plot_velocity_profile
    set(fig_main, 'Position', [100, 100, 1100, 700]);
    ax_main = subplot(2,2,[1,2], 'Parent', fig_main);
else
    set(fig_main, 'Position', [100, 100, 800, 700]);
    ax_main = subplot(2,1,1, 'Parent', fig_main);
end
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

% Quiver arrows (at z=1 to stay above pcolor surface)
u_sub0 = zeros(size(X_sub));
v_sub0 = zeros(size(X_sub));
h_quiver = quiver3(ax_main, X_sub, Y_sub, ones(size(X_sub)), u_sub0, v_sub0, zeros(size(X_sub)), scale_factor, 'k', 'LineWidth', 1);

% Boundary markers (at z=1)
h_top = plot3(ax_main, NaN, NaN, NaN, 'r.', 'MarkerSize', 10);
h_bot = plot3(ax_main, NaN, NaN, NaN, 'r.', 'MarkerSize', 10);

% EC (elastic-contractile) connection lines between opposing nodes (at z=1)
h_ec = plot3(ax_main, NaN, NaN, NaN, 'm-', 'LineWidth', 0.5);

% Distinct anchor-node markers (green squares, larger, at z=1)
h_anchor_top = plot3(ax_main, NaN, NaN, NaN, 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
h_anchor_bot = plot3(ax_main, NaN, NaN, NaN, 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

% Title handle
h_title_main = title(ax_main, '');

%% Pressure Centerline Subplot Setup
centerline_y_idx = round(Ny / 2);
if plot_velocity_profile
    ax_pressure = subplot(2,2,3, 'Parent', fig_main);
else
    ax_pressure = subplot(2,1,2, 'Parent', fig_main);
end
hold(ax_pressure, 'on');
grid(ax_pressure, 'on');
box(ax_pressure, 'on');
set(ax_pressure, 'FontSize', 12);
xlabel(ax_pressure, 'X (\mum)');
ylabel(ax_pressure, 'Pressure');
h_pressure_line = plot(ax_pressure, x_coords, zeros(1, Nx), 'b-', 'LineWidth', 2);
h_title_pressure = title(ax_pressure, 'Pressure along vessel centerline (\mum^2/s^2)');
xlim(ax_pressure, [0 Lx]);

fprintf('Starting velocity & pressure animation and recording...\n');

%% Poiseuille Cross-Section Setup
% y_top_wall = Y_top(1, 2);
% y_bot_wall = Y_bot(1, 2);
% H_channel  = y_top_wall - y_bot_wall;
% 
% u_theory = 2 * G ./ (2 * mu) .* max(y_coords - y_bot_wall, 0) .* max(y_top_wall - y_coords, 0);
% 
% U_mean_theory = G * H_channel^2 / (8 * mu);
% Re_channel    = rho * U_mean_theory * H_channel / mu;
% fprintf('Reynolds number Re = %.4f\n', Re_channel);
% 
% cross_x_idx = round(Nx / 2);
% 
% % --- Cross-section figure: pre-create plot objects ---
% fig_cross = figure('Name', 'Poiseuille Velocity Profile — Cross-Section', 'Color', 'w', 'Visible', 'off');
% set(fig_cross, 'Position', [950, 100, 420, 520], 'Resize', 'off');
% ax_cross = axes(fig_cross);
% hold(ax_cross, 'on');
% grid(ax_cross, 'on');
% 
% h_sim   = plot(ax_cross, NaN(1,Ny), y_coords, 'b-',  'LineWidth', 2.5, 'DisplayName', 'Simulated u_x');
% h_theo  = plot(ax_cross, u_theory,  y_coords, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Theory: G/(2\mu)(y-y_b)(y_t-y)');
% yline(ax_cross, y_top_wall, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% yline(ax_cross, y_bot_wall, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% xlabel(ax_cross, 'u_x  (m/s)');
% ylabel(ax_cross, 'y  (m)');
% ylim(ax_cross, [0 Ly]);
% legend(ax_cross, 'Location', 'northeast');   % fixed location (faster than 'best')
% h_title_cross = title(ax_cross, '');

%% Velocity Cross-Section Setup (conditional)
if plot_velocity_profile
    cross_x_idx = round(Nx / 2);
    y_top_wall_cs = Y_top(1, 2);
    y_bot_wall_cs = Y_bot(1, 2);
    H_channel_cs  = y_top_wall_cs - y_bot_wall_cs;
    u_theory = G ./ (2 * mu) .* max(y_coords - y_bot_wall_cs, 0) .* max(y_top_wall_cs - y_coords, 0);

    ax_cross = subplot(2,2,4, 'Parent', fig_main);
    hold(ax_cross, 'on');
    grid(ax_cross, 'on');
    box(ax_cross, 'on');
    set(ax_cross, 'FontSize', 12);
    h_cross_sim  = plot(ax_cross, NaN(1,Ny), y_coords, 'b-',  'LineWidth', 2.5, 'DisplayName', 'Simulated u_x');
    h_cross_theo = plot(ax_cross, u_theory,  y_coords, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Poiseuille theory');
    yline(ax_cross, y_top_wall_cs, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(ax_cross, y_bot_wall_cs, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel(ax_cross, 'u_x (\mum/s)');
    ylabel(ax_cross, 'y (\mum)');
    ylim(ax_cross, [0 Ly]);
    legend(ax_cross, 'Location', 'northeast');
    h_title_cross = title(ax_cross, 'Velocity profile at mid-x');
end

% Pre-compute grid spacings for vorticity
dy = Ly / (Ny - 1);
dx = Lx / (Nx - 1);

% Force one initial render so getframe works on invisible figures
drawnow;

%% Compute Reynolds number from saved data
y_top_wall = Y_top(1, 2);
y_bot_wall = Y_bot(1, 2);
H_channel  = y_top_wall - y_bot_wall;
U_mean_theory = G * H_channel^2 / (8 * mu);
Re_channel    = rho * U_mean_theory * H_channel / mu;
fprintf('Reynolds number Re = %.4f\n', Re_channel);

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
    set(h_top, 'XData', X_current_top(:,1), 'YData', X_current_top(:,2), 'ZData', ones(size(X_current_top,1),1));
    set(h_bot, 'XData', X_current_bot(:,1), 'YData', X_current_bot(:,2), 'ZData', ones(size(X_current_bot,1),1));

    % ---- Update EC connection lines ----
    ec_stride = 5;   % draw every 5th connection to avoid clutter
    ec_idx = 1:ec_stride:size(X_current_top, 1);
    x_ec = [X_current_top(ec_idx,1), X_current_bot(ec_idx,1), NaN(numel(ec_idx),1)]';
    y_ec = [X_current_top(ec_idx,2), X_current_bot(ec_idx,2), NaN(numel(ec_idx),1)]';
    z_ec = [ones(numel(ec_idx),1), ones(numel(ec_idx),1), NaN(numel(ec_idx),1)]';
    set(h_ec, 'XData', x_ec(:), 'YData', y_ec(:), 'ZData', z_ec(:));

    % ---- Update anchor markers ----
    if exist('fix_ind', 'var')
        n_fix = numel(fix_ind);
        set(h_anchor_top, 'XData', X_current_top(fix_ind,1), 'YData', X_current_top(fix_ind,2), 'ZData', ones(n_fix,1));
        set(h_anchor_bot, 'XData', X_current_bot(fix_ind,1), 'YData', X_current_bot(fix_ind,2), 'ZData', ones(n_fix,1));
    end

    current_time = k * sample_rate * dt;
    set(h_title_main, 'String', sprintf('t = %.3f s | Re = %.2f | \\mu = %g | G = %g | K_M = %g', ...
          current_time, Re_channel, mu, G, K_M));

    % ---- Update pressure centerline plot ----
    if exist('p_history', 'var')
        p_centerline = p_history(:, centerline_y_idx, k);
        set(h_pressure_line, 'YData', p_centerline);
        p_range = max(p_centerline) - min(p_centerline);
        if p_range > 0
            ylim(ax_pressure, [min(p_centerline) - 0.1*p_range, max(p_centerline) + 0.1*p_range]);
        end
    end

    % ---- Update velocity cross-section plot ----
    if plot_velocity_profile
        u_cross = u_vel(cross_x_idx, :);
        set(h_cross_sim, 'XData', u_cross);
        u_max_plot = max([max(u_theory), max(abs(u_cross)), 1e-9]) * 1.15;
        xlim(ax_cross, [-u_max_plot * 0.05, u_max_plot]);
        set(h_title_cross, 'String', sprintf('u_x at x = %.1f \\mum', x_coords(cross_x_idx)));
    end

    % Capture and write frame
    drawnow limitrate;
    frame = getframe(fig_main);
    writeVideo(v, frame);

    % ---- Update cross-section objects (no cla) ----
    % u_cross = u_vel(cross_x_idx, :);
    % set(h_sim, 'XData', u_cross);
    %
    % % Vorticity at cross-section: omega = dv/dx - du/dy
    % dvdx_cross  = gradient(v_vel(cross_x_idx, :), dx);
    % dudy_cross  = gradient(u_vel(cross_x_idx, :), dy);
    % omega_cross = dvdx_cross - dudy_cross;
    % omega_peak  = max(abs(omega_cross));
    %
    % set(h_title_cross, 'String', sprintf('x = %.2f m,  t = %.3f s  |  Re = %.2f  |  |\\omega|_{max} = %.3g', ...
    %       x_coords(cross_x_idx), current_time, Re_channel, omega_peak));
    %
    % u_max_plot = max([max(u_theory), max(abs(u_cross)), 1e-9]) * 1.15;
    % xlim(ax_cross, [-u_max_plot * 0.05, u_max_plot]);
    %
    % drawnow;
    % frame_cross = getframe(fig_cross);
    % writeVideo(v_cross, frame_cross);
end

fprintf('\n');   % finish the progress line
close(v);
% close(v_cross);
fprintf('Animation finished. Video saved to:\n  %s\n', video_filename);