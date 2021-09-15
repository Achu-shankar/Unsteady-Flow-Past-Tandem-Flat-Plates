%% User defined variables
%%
U            = [1 0]; % freestream velocity (in global coordinate)
time         = 10;     % time of simulation in seconds
dt           = 0.01;  % time interval size
N_panels    = 100;   
%% Plate A
chord_A        = 1;              %chord length 
X_0A          = [0 1];          % Leading edge position (in global coordinate)   
P_rotA        = [0.5*chord_A,0]; %Point of rotation (in body coordinates)

% kinematics (global coordinates)
omega_rotA   = 1;
alpha_func_A  = @(t)  5*pi*cos(omega_rotA*t)/180; 
omega_alpha_A = @(t) -5*omega_rotA*pi*sin(omega_rotA*t)/180; 

%% Plate B
chord_B        = 1;             %chord length 
X_0B          = [6 1];         % Leading edge position (in global coordinate)
P_rotB        = [0.5*chord_B,0];%Point of rotation (in body coordinates)

% kinematics (global coordinates)
omega_rotB   = 1;
alpha_func_B  = @(t)  5*pi*cos(omega_rotB*t)/180; 
omega_alpha_B = @(t) -5*omega_rotB*pi*sin(omega_rotB*t)/180; 

