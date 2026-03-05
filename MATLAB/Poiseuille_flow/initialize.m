clear;
clc;
close all;

%% Parameters

% Change this part
% Time spacing
Tmax = 15; % total time (s)
dt = 0.0005; % Time spacing (s)
n_steps = floor(Tmax / dt); % Number of steps
K_tar = 100; % Stiffnes constant between ghost and real mass (large for rigid walls)
K_mem = 100;   % Stiffness constant between mass points (0 = flat rigid wall, no membrane elasticity)
use_parallel = true; % Set true to start a parallel pool (speeds up internal MATLAB ops via implicit parallelism)

% Add a perturbation force
t_per = Tmax /2; % perturbation time
per_step = t_per / dt;
per_dur = 0.05; % Perturbation duration
per_time_steps = per_dur / dt;


ratio = K_tar / max(K_mem, eps); % the stiffness ratio (guard against K_mem = 0)

% Fluid Parameters
rho = 1;   % density
mu = 0.02;  % Dynamic viscosity (large for low-Re laminar Poiseuille flow)

% Driving force 
G = 1;  % Driving body force (pressure gradient surrogate); smaller G with high mu gives clean parabola

% Grid size
Lx = 4.0; % x length (m)
Ly = 2.0; % y length (m)
Nx = 128; % Number of grid points
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

%% Material Parameters
Num_b = 200; % Number of grids points on boundary
Tube_ind_f = [Num_b, (1 : Num_b - 1)];
Tubd_ind_b = [(2 : Num_b), 1];

% 这里要改！！！！！！
% dtheta = 2 * pi / Num_b;
dtheta = Ly / Num_b;

% Initial Material Position
X_top = zeros(Num_b, 2); % Initialization
X_bot = X_top;

y_displace = 0.4; % Displacement from central line
% Set initial position
X_top(:, 1) = linspace(0, Lx, Num_b);
X_bot(:, 1) = X_top(:, 1);
X_top(:, 2) = Ly / 2 + y_displace;
X_bot(:, 2) = Ly / 2 - y_displace;

% We fix the start and end of the tube points. 
fix_ind = [1, Num_b];
top_fix = X_top(fix_ind, :);
bot_fix = X_bot(fix_ind, :);

% Fluid field initialization
u = zeros(Nx, Ny, 2); % In this simulation, we prescribe a total body force on the fluid. So the initial velocity is zero. 

% Driving force definition
% Uniform x-body force across the entire domain (Poiseuille: driven by constant pressure gradient G)
f_drive = zeros(Nx, Ny, 2);
% f_drive(:, :, 1) = G;  % Apply G uniformly in x-direction 

% Apply G only between the walls (tube interior)
for j = 1:Ny
    y = (j-1)*dy;
    if (y > X_bot(1,2)) && (y < X_top(1,2))
        f_drive(:, j, 1) = G;   % uniform across all x-columns, but only inside the tube
    end
end


% Ghost points
Y_top = X_top;
Y_bot = X_bot;

% Animation set up
frame = 30;
sample_rate = 1 / frame / dt;
% Initialize storing
n_saved = floor(n_steps / sample_rate) + 1;


%% Force Perturbation
f_perturbation = zeros(Nx, Ny, 2);  % Bug fix: must be (Nx, Ny, 2) to match f dimensions
per_ind_y = 1 / 2 * Ny;
per_x_range = 10;
per_ind_x_start = 1 / 2 * Nx - per_x_range / 2;
per_ind_x_end = per_ind_x_start + per_x_range;
per_force = 0.4;

% Add this force to the matrix
f_perturbation(per_ind_x_start:per_ind_x_end, per_ind_y) = per_force;


