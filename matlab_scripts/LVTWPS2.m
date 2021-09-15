% Code to solve flow past an oscill. flat plate using Lumped vortex method
% Unsteady Aerodynamics - Course Project
% Written by : Alan John Varghese and Achu Shankar
% Oct 2019
clearvars;
clc;
tic
%% Importing user defined variables
USER_VAR

%% Importing private variables
PRIVATE_VAR

%% Time marching

%progress bar
disp("computing...")
% h = waitbar(0,'computing...');
% steps = length(t);


for tstep = 1:length(t)
    % body-global transformations     
    B_G_var_trans
    %influence coefficient
    Influence_coeff_gen
    %RHS B
    RHS
    
    % solving for the unknown circulation values
    X = mldivide(A,B) ;
    
    
    %updating  wake strength
    w_str_A(tstep) = X(2*N_panels+1,1);
    w_str_B(tstep) = X(2*N_panels+2,1);
    wake_gA(tstep,:) = n_wake_gA;
    wake_gB(tstep,:) = n_wake_gB;

    
    %updating wake position
    u = zeros(tstep,1);
    v = zeros(tstep,1);
    for i=1:tstep
        [u_g,v_g]    = VORT_VEL(w_str_A,wake_gA(i,1),wake_gA(i,2),wake_gA(:,1),wake_gA(:,2));
        [u_g1,v_g1]  = VORT_VEL(X(1:N_panels),wake_gA(i,1),wake_gA(i,2),vort_gA(:,1),vort_gA(:,2));
        [u_g2,v_g2]  = VORT_VEL(w_str_B,wake_gA(i,1),wake_gA(i,2),wake_gB(:,1),wake_gB(:,2));
        [u_g3,v_g3]  = VORT_VEL(X(N_panels+1:2*N_panels),wake_gA(i,1),wake_gA(i,2),vort_gB(:,1),vort_gB(:,2));
                u(i)   = sum([u_g;u_g1;u_g2;u_g3]);
                v(i)   = sum([v_g;v_g1;v_g2;v_g3]);
    end
    wake_gA(1:tstep-1,1) = wake_gA(1:tstep-1,1) + (u(1:tstep-1,1)+U(1))*dt;
    wake_gA(1:tstep-1,2) = wake_gA(1:tstep-1,2) + (v(1:tstep-1,1)+U(2))*dt;
    wake_xgtA(1:tstep-1,tstep) = wake_gA(1:tstep-1,1);
    wake_ygtA(1:tstep-1,tstep) = wake_gA(1:tstep-1,2);
    
    
    u = zeros(tstep,1);
    v = zeros(tstep,1);
    for i=1:tstep
        [u_g,v_g]    = VORT_VEL(w_str_A,wake_gB(i,1),wake_gB(i,2),wake_gA(:,1),wake_gA(:,2));
        [u_g1,v_g1]  = VORT_VEL(X(1:N_panels),wake_gB(i,1),wake_gB(i,2),vort_gA(:,1),vort_gA(:,2));
        [u_g2,v_g2]  = VORT_VEL(w_str_B,wake_gB(i,1),wake_gB(i,2),wake_gB(:,1),wake_gB(:,2));
        [u_g3,v_g3]  = VORT_VEL(X(N_panels+1:2*N_panels),wake_gB(i,1),wake_gB(i,2),vort_gB(:,1),vort_gB(:,2));
                u(i)   = sum([u_g;u_g1;u_g2;u_g3]);
                v(i)   = sum([v_g;v_g1;v_g2;v_g3]);
    end
    wake_gB(1:tstep-1,1) = wake_gB(1:tstep-1,1) + (u(1:tstep-1,1)+U(1))*dt;
    wake_gB(1:tstep-1,2) = wake_gB(1:tstep-1,2) + (v(1:tstep-1,1)+U(2))*dt;
    wake_xgtB(1:tstep-1,tstep) = wake_gB(1:tstep-1,1);
    wake_ygtB(1:tstep-1,tstep) = wake_gB(1:tstep-1,2);
%     waitbar(tstep / steps)
    gamma_pre_A = sum(X(1:N_panels,1));
    gamma_pre_B = sum(X(N_panels+1:2*N_panels,1));
end
toc
 disp("done:)")

%% plotting
figure
for i=2:length(t)
    xA = BODY_GLOB(coll_A,alpha_func_A(t(i)),X_0A);
    xB = BODY_GLOB(coll_B,alpha_func_B(t(i)),X_0B);
    plot(xA(:,1),xA(:,2),"r","LineWidth",3)
    hold on
    plot(xB(:,1),xB(:,2),"r","LineWidth",3)
    scatter( wake_xgtA(1:i-1,i), wake_ygtA(1:i-1,i),10,'b')
    hold on
    scatter(wake_xgtB(1:i-1,i),wake_ygtB(1:i-1,i),10,'b')
    hold off
    axis([-2,30,-7,7])
    pause(.001)
%     drawnow   
end