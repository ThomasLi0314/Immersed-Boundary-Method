clear;
clc;
close all;

%% Parameters

% Change this part
% Time spacing
Tmax = 20; % total time (s)
dt = 0.001; % Time spacing (s)
n_steps = floor(Tmax / dt); % Number of steps
K_tar = 1e3; % Tether stiffness at anchor nodes (very large to pin endpoints)
K_mem = 1e3;   % Stiffness constant between mass points (0 = flat rigid wall, no membrane elasticity)
use_parallel = true; % Set true to start a parallel pool (speeds up internal MATLAB ops via implicit parallelism)

ratio = K_tar / max(K_mem, eps); % the stiffness ratio (guard against K_mem = 0)

%% Myogenic (Elastic-Contractile) Parameters — Arthurs et al. (1998), Eq. 60
K_M     = 0;     % Myogenic spring stiffness (constant for initial impl.)
gamma_M = 0;   % Viscous damping coefficient for diameter rate-of-change
K_hor   = 0;   % Weak horizontal restoring stiffness (prevents downstream drift)

% K_M = 0;
% gamma_M = 0
% K_hor = 0;
% D0 and D_prev are set after y_displace is defined (see below)

%% Add a perturbation force
% t_per = Tmax /2; % perturbation time
% per_step = t_per / dt;
% per_dur = 0.05; % Perturbation duration
% per_time_steps = per_dur / dt;

% Force Perturbation
% f_perturbation = zeros(Nx, Ny, 2);  % Bug fix: must be (Nx, Ny, 2) to match f dimensions
% per_ind_y = 1 / 2 * Ny;
% per_x_range = 10;
% per_ind_x_start = 1 / 2 * Nx - per_x_range / 2;
% per_ind_x_end = per_ind_x_start + per_x_range;
% per_force = 0.4;

% % Add this force to the matrix
% f_perturbation(per_ind_x_start:per_ind_x_end, per_ind_y) = per_force;

% Fluid Parameters (code units: 1 length unit = 1 µm, time in seconds)
%   Physical: ρ_blood ≈ 1060 kg/m³, µ_blood ≈ 3.5e-3 Pa·s
%   Scaled so that Re ≈ 0.1 and T_viscous ≈ 1.25 s (fits in Tmax)
rho = 1;      % density (code units)
mu  = 320;    % dynamic viscosity (code units; ν = µ/ρ = 320 µm²/s)

% Driving force
%   G chosen for Re ≈ 0.1:  Re = ρ·G·D³/(8·µ²)
G = 10;       % body force (code units)

% Grid size — brain arteriole: ~200 µm long, domain 100 µm wide
Lx = 200;   % vessel length (µm)
Ly = 100;   % domain width  (µm)
Nx = 256;   % grid points in x
Ny = Nx * Ly / Lx;  % = 64

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
ds0    = Lx / (Num_b - 1);  % Rest length of each spring segment (initial node spacing)

% Initial Material Position
X_top = zeros(Num_b, 2); % Initialization
X_bot = X_top;

y_displace = 10;  % half-channel width (µm) → D = 20 µm (brain arteriole diameter)
% Set initial position
X_top(:, 1) = linspace(0, Lx, Num_b);
X_bot(:, 1) = X_top(:, 1);
X_top(:, 2) = Ly / 2 + y_displace;
X_bot(:, 2) = Ly / 2 - y_displace;

%% Fixed index

% Test1 first and last two
% We fix the start and end of the tube points. 
% fix_ind = [1, 2, Num_b-1, Num_b];   % Tether anchors at first/last 2 nodes only
% fix_ind = [1, 2];

% Test2 all
% fix index for possieul flow
% fix_ind = (1 : Num_b);

% Test3 middle
% fix_ind = [Num_b / 2, Num_b / 2 + 1];

% Test 4 Fix 10 evenly spaced indices along the full length
fix_ind = round(linspace(1, Num_b, 10));

% Test 5: Exponential decay tether from K_tar (endpoints) to K_wall (interior)
%   - K_tar  = 1e6 at node 1 and node Num_b (strong anchoring)
%   - K_wall = 1e3 in the deep interior   (baseline tissue rigidity)
%   - Smooth exponential transition between the two
K_wall = 1e1;        % interior (minimum) tether stiffness
decay_len = 15;      % characteristic decay length (in node indices)

% Distance from each end (0-based)
dist_from_start = (0 : Num_b - 1)';
dist_from_end   = (Num_b - 1 : -1 : 0)';

% Exponential profiles from each end, take the stronger of the two
K_from_start = (K_tar - K_wall) * exp(-dist_from_start / decay_len) + K_wall;
K_from_end   = (K_tar - K_wall) * exp(-dist_from_end   / decay_len) + K_wall;
K_tar_vec    = max(K_from_start, K_from_end)';   % 1 x Num_b


%% Spring force part
% set spring force only at the potision whe

%% Myogenic reference diameter & state (depends on y_displace, defined above)
D0     = 2 * y_displace;             % Reference diameter (initial wall separation)
D_prev = D0 * ones(Num_b, 1);        % Previous-step diameter for ∂_t D approx.
top_fix = X_top(fix_ind, :);
bot_fix = X_bot(fix_ind, :);

%% Fluid field initialization
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

% Test, add drive force at the middle
for j = 1:Ny
    y = (j-1)*dy;
    if (y > X_bot(1,2)) && (y < X_top(1,2))
        f_drive(Num_b / 2, j, 1) = G;   % only at the middle
    end
end


% Ghost points
Y_top = X_top;
Y_bot = X_bot;

% Animation set up
frame = 30;
sample_rate = round(1 / frame / dt);  % must be integer for mod() check in main loop
% Initialize storing
n_saved = floor(n_steps / sample_rate) + 1;



