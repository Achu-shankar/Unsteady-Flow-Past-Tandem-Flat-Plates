% Code to solve flow past an oscill. flat plate using Lumped vortex method
% Unsteady Aerodynamics - Course Project
% Written by : Alan John Varghese and Achu Shankar
% Oct 2019
clearvars;
clc;


%% Importing user defined variables
USER_VAR

%% Private variables 
PRIVATE_VAR


%% time marching ( finding vortex strengths, update B, find new wake vortex and updates position of previous ones )
%progress bar
disp("computing...")
h = waitbar(0,'computing...');
steps = length(t);

for timestep = 1:length(t)
    
    % calculating all global coordinates
    [normal_gx,normal_gy] = BODY_GLOB(normal(1),normal(2),alpha_func(t(timestep)),0,0);
     normal_g             = [normal_gx,normal_gy];
    [xcoll_g,ycoll_g]     = BODY_GLOB(coll_b(:,1),coll_b(:,2),alpha_func(t(timestep)),X_0A,Y_0A);
    [xvort_g,yvort_g]     = BODY_GLOB(vort_b(:,1),vort_b(:,2),alpha_func(t(timestep)),X_0A,Y_0A);
    [n_wake_xg,n_wake_yg] = BODY_GLOB(n_wake(1),n_wake(2),alpha_func(t(timestep)),X_0A,Y_0A);
    [rotvel_gx,rotvel_gy] = BODY_GLOB(0,omega_alpha(timestep)*(coll_b(:,1)-xrot),alpha_func(t(timestep)),0,0);
     rotvel_g             = [rotvel_gx,rotvel_gy];
    % Influence matrix BODY_GLOB
    for i = 1:N_panels
        [uvelcoll_i,vvelcoll_i] = VORT_VEL(1,xcoll_g(i),ycoll_g(i),xvort_g,yvort_g);
        A(i,1:end-1) = normal_g*[uvelcoll_i,vvelcoll_i]';
    end
    [uvelcoll_w,vvelcoll_w] = VORT_VEL(1,xcoll_g,ycoll_g,n_wake_xg,n_wake_yg); %velocity induced at all coll points due to the new wake 
    A(1:end-1,N_panels+1)   = normal_g*[uvelcoll_w,vvelcoll_w]';


% Establishing the RHS B (body fixed)
    B(1:end-1,1) = -(U+v_h_func(timestep)+rotvel_g)*normal_g';

% including wake induced velocity in RHS
    C = B;
    for i=1:N_panels
        [uw,vw]                 = VORT_VEL(wakestrength,xcoll_g(i),ycoll_g(i),wakeglob_x,wakeglob_y);
        velwnorm                = [uw,vw]*normal_g';
        C(i,1)                  = C(i,1)-sum(velwnorm);
    end
    C(end,1) = gamma_prevt;
    
    %solving the linear equation AX = B
    X = mldivide(A,C);

    %updating the wake locations
%     [wakeglob_x(timestep),wakeglob_y(timestep)] = BODY_GLOB(newwake_locx,...
%                                                   newwake_locy,alpha_func(t(timestep)),h_func(t(timestep)),xrot,yrot);
    wakeglob_x(timestep) = n_wake_xg;
    wakeglob_y(timestep) = n_wake_yg;
    
    wakestrength(timestep) = X(end,1);
    uvel_w = zeros(timestep,1);
    vvel_w = zeros(timestep,1);
    for i=1:timestep
        [uvel_wsig,vvel_wsig] = VORT_VEL(wakestrength,wakeglob_x(i),wakeglob_y(i),wakeglob_x,wakeglob_y);
        [uvel_wig,vvel_wig]   = VORT_VEL(X(1:end-1),wakeglob_x(i),wakeglob_y(i),xvort_g,yvort_g);
        uvel_w(i) = sum([uvel_wsig;uvel_wig]);
        vvel_w(i) = sum([vvel_wsig;vvel_wig]);
    end
     wakeglob_x(1:timestep-1,1) =  wakeglob_x(1:timestep-1,1) + (uvel_w(1:timestep-1,1)+U(1))*dt ; 
     wakeglob_y(1:timestep-1,1) =  wakeglob_y(1:timestep-1,1) + (vvel_w(1:timestep-1,1)+U(2))*dt ;
     
     wakeglob_xt(1:timestep-1,timestep) = wakeglob_x(1:timestep-1,1);
     wakeglob_yt(1:timestep-1,timestep) = wakeglob_y(1:timestep-1,1);
     gamma_prevt = sum(X(1:end-1,1));
     
     waitbar(timestep / steps)
end
disp("done:)")

%% plotting
figure
for i=2:length(t)
    [x,y] = BODY_GLOB(coll_b(:,1),coll_b(:,2),alpha_func(t(i)),X_0A,Y_0A);
    plot(x,y,"b","LineWidth",3)
    hold on
    scatter(wakeglob_xt(1:i-1,i),wakeglob_yt(1:i-1,i),1,'r')
    hold off
    axis([-2,30,-7,7])
    pause(.001)
%     drawnow   
end