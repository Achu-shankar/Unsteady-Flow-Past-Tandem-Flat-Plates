t = 0:dt:time;
% x and y coordinates of the panels(in body coordinates) 
x_0 = linspace(0,chord_A,N_panels+1)';

%Vorticity points and collocation points (in body coordinates)
vort_A = [x_0(1:end-1) + (x_0(2:end)-x_0(1:end-1))*0.25 , zeros(N_panels,1)];
coll_A = [x_0(1:end-1) + (x_0(2:end)-x_0(1:end-1))*0.75 , zeros(N_panels,1)];


vort_B = [x_0(1:end-1) + (x_0(2:end)-x_0(1:end-1))*0.25 , zeros(N_panels,1)];
coll_B = [x_0(1:end-1) + (x_0(2:end)-x_0(1:end-1))*0.75 , zeros(N_panels,1)];

% Wakestrength and global location 
w_str_A    = zeros(length(t),1); 
wake_gA    = zeros(length(t),2);
wake_xgtA  = zeros(length(t),length(t));
wake_ygtA  = zeros(length(t),length(t));

w_str_B    = zeros(length(t),1); 
wake_gB    = zeros(length(t),2);
wake_xgtB  = zeros(length(t),length(t));
wake_ygtB  = zeros(length(t),length(t));
                    
norm_A  = [0,1]; % normal of each panel  (body fixed)
norm_B  = [0,1];

gamma_pre_A = 0; % total gamma of due to the node vortices in the previous time step 
gamma_pre_B = 0;

A           = zeros(2*N_panels+2,2*N_panels+2); %influence coefficient matrix
B           = zeros(2*N_panels+2,1);

n_wake_A = [chord_A+ 0.1*U(1)*dt , 0]; %location of the newest wake (body fixed)
n_wake_B = [chord_B+ 0.1*U(1)*dt , 0];

A(2*N_panels+1,1:N_panels) = 1;
A(2*N_panels+1,2*N_panels+1) = 1;
A(2*N_panels+2,N_panels+1:2*N_panels) = 1;
A(2*N_panels+2,2*N_panels+2) = 1;

