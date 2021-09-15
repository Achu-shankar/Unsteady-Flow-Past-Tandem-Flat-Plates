% Code to solve flow past an oscill. flat plate using Lumped vortex method
% Unsteady Aerodynamics - Course Project
% Written by : Alan John Varghese and Achu Shankar
% Oct 2019
clearvars;
clc;

%% Importing user defined variables
variables


%% Private variables 
% all vectors are defined as column vectors

% x and y coordinates of the panels(in body coordinates) 
x_0 = linspace(0,chord_length,N_panels+1)';
y_0 = zeros(N_panels+1,1);

%Vorticity points and collocation points of plate A(in body coordinates of A)
vort_bA(:,1) = x_0(1:end-1) + (x_0(2:end)-x_0(1:end-1))*0.25;
vort_bA(:,2) = zeros(N_panels,1);
coll_bA(:,1) = x_0(1:end-1) + (x_0(2:end)-x_0(1:end-1))*0.75;
coll_bA(:,2) = zeros(N_panels,1);
%Vorticity points and collocation points of plate B(in body coordinates of B)
vort_bB(:,1) = x_0(1:end-1) + (x_0(2:end)-x_0(1:end-1))*0.25;
vort_bB(:,2) = zeros(N_panels,1);
coll_bB(:,1) = x_0(1:end-1) + (x_0(2:end)-x_0(1:end-1))*0.75;
coll_bB(:,2) = zeros(N_panels,1);

% Wakestrength and global location of plate A
wakestrengthA = zeros(length(t),1); 
wakeglob_Ax   = zeros(length(t),1);
wakeglob_Ay   = zeros(length(t),1);
wakeglob_Axt   = zeros(length(t),length(t));
wakeglob_Ayt   = zeros(length(t),length(t));
% Wakestrength and global location of plate B
wakestrengthB = zeros(length(t),1); 
wakeglob_Bx   = zeros(length(t),1);
wakeglob_By   = zeros(length(t),1);
wakeglob_Bxt   = zeros(length(t),length(t));
wakeglob_Byt   = zeros(length(t),length(t));

%location of the newest wake of plate A(body fixed)
n_wakeA(1) = chord_length + 0.1*U(1)*dt; % should it include the component of U
n_wakeA(2) = 0;                           %vel in that case it has to be updated each time step.

%location of the newest wake of plate B (body fixed)
n_wakeB(1) = chord_length + 0.1*U(1)*dt; % should it include the component of U
n_wakeB(2) = 0;                           %vel in that case it has to be updated each time step.

normA  = [0,1]; % normal of each panel  (body fixed)
normB  = [0,1]; % normal of each panel  (body fixed)
gamma_prevtA = 0; % total gamma of due to the node vortices in the previous time step
gamma_prevtB = 0; % total gamma of due to the node vortices in the previous time step

A        = ones(N_panels+1,N_panels+1); %influence coefficient matrix
B        = zeros(N_panels+1,1);

% A        = zeros(2*N_panels+2,2*N_panels+2); %influence coefficient matrix
% B        = zeros(2*N_panels+2,1);
% 
% A(2*N_panels+1,1:N_panels)            = 1;
% A(2*N_panels+1,2*N_panels+1)          = 1;
% A(2*N_panels+2,N_panels+1:2*N_panels) = 1;
% A(2*N_panels+2,2*N_panels+2)          = 1;


%% time marching
disp("computing...")
% h = waitbar(0,'computing...');
steps = length(t);
%time marching ( find vortex strengths, update B, find new wake vortex and updates position of previous ones )
for timestep = 1:length(t)
    % calculating all global coordinates
    [norm_Agx,norm_Agy]   = BODY_GLOB(normA(1),normA(2),alpha_func(t(timestep)),0,0);
     norm_Ag            = [norm_Agx,norm_Agy];
    [norm_Bgx,norm_Bgy]   = BODY_GLOB(normB(1),normB(2),alpha_func(t(timestep)),0,0);
     norm_Bg            = [norm_Bgx,norm_Bgy];
     
    [xcoll_Ag,ycoll_Ag]     = BODY_GLOB(coll_bA(:,1),coll_bA(:,2),alpha_func(t(timestep)),X_0A,Y_0A);
    [xvort_Ag,yvort_Ag]     = BODY_GLOB(vort_bA(:,1),vort_bA(:,2),alpha_func(t(timestep)),X_0A,Y_0A);
    
    [xcoll_Bg,ycoll_Bg]     = BODY_GLOB(coll_bB(:,1),coll_bB(:,2),alpha_func(t(timestep)),X_0B,Y_0B);
    [xvort_Bg,yvort_Bg]     = BODY_GLOB(vort_bB(:,1),vort_bB(:,2),alpha_func(t(timestep)),X_0B,Y_0B);
    
    [n_wake_Axg,n_wake_Ayg]   = BODY_GLOB(n_wakeA(1),n_wakeA(2),alpha_func(t(timestep)),X_0A,Y_0A);
    [n_wake_Bxg,n_wake_Byg]   = BODY_GLOB(n_wakeB(1),n_wakeB(2),alpha_func(t(timestep)),X_0B,Y_0B);
    
    
    [rotvel_Agx,rotvel_Agy] = BODY_GLOB(0,omega_alpha(timestep)*(coll_bA(:,1)-xrot),alpha_func(t(timestep)),0,0);
     rotvel_Ag             = [rotvel_Agx,rotvel_Agy];
     
    [rotvel_Bgx,rotvel_Bgy] = BODY_GLOB(0,omega_alpha(timestep)*(coll_bB(:,1)-xrot),alpha_func(t(timestep)),0,0);
     rotvel_Bg             = [rotvel_Bgx,rotvel_Bgy];
     
%     Influence matrix  F1
    for i = 1:N_panels
        [uvelcoll_i,vvelcoll_i] = VORT_VEL(1,xcoll_Ag(i),ycoll_Ag(i),xvort_Ag,yvort_Ag);
        A(i,1:N_panels) = norm_Ag *[uvelcoll_i,vvelcoll_i]';
        
    end
    % Influence matrix  F2
    for i = 1:N_panels
        [uvelcoll_i,vvelcoll_i] = VORT_VEL(1,xcoll_Ag(i),ycoll_Ag(i),xvort_Bg,yvort_Bg);
        A(i, N_panels+1:2*N_panels) = norm_Ag *[uvelcoll_i,vvelcoll_i]';

    end
    
%     Influence matrix  F3
    [uvelcoll_w,vvelcoll_w] = VORT_VEL(1,xcoll_Ag,ycoll_Ag,n_wake_Axg,n_wake_Ayg); 
    A(1:N_panels,2*N_panels+1)   = norm_Ag*[uvelcoll_w,vvelcoll_w]';

    % Influence matrix  F4
    [uvelcoll_w,vvelcoll_w] = VORT_VEL(1,xcoll_Ag,ycoll_Ag,n_wake_Bxg,n_wake_Byg); 
    A(1:N_panels,2*N_panels+2)   = norm_Ag*[uvelcoll_w,vvelcoll_w]';

    % Influence matrix  G1
    for i = N_panels+1:2*N_panels
        [uvelcoll_i,vvelcoll_i] = VORT_VEL(1,xcoll_Bg(i- N_panels),ycoll_Bg(i- N_panels),xvort_Ag,yvort_Ag);
        A(i,1:N_panels) = norm_Bg *[uvelcoll_i,vvelcoll_i]';
    end
    
     % Influence matrix  G2
    for i = N_panels+1:2*N_panels
        [uvelcoll_i,vvelcoll_i] = VORT_VEL(1,xcoll_Bg(i- N_panels),ycoll_Bg(i- N_panels),xvort_Bg,yvort_Bg);
        A(i,N_panels+1:2*N_panels) = norm_Bg *[uvelcoll_i,vvelcoll_i]';
    end
    
     % Influence matrix  G3
    [uvelcoll_w,vvelcoll_w] = VORT_VEL(1,xcoll_Bg,ycoll_Bg,n_wake_Axg,n_wake_Ayg); 
    A(N_panels+1:2*N_panels,2*N_panels+1)   = norm_Bg*[uvelcoll_w,vvelcoll_w]';
    
      % Influence matrix  G4
    [uvelcoll_w,vvelcoll_w] = VORT_VEL(1,xcoll_Bg,ycoll_Bg,n_wake_Bxg,n_wake_Byg); 
    A(N_panels+1:2*N_panels,2*N_panels+2)   = norm_Bg*[uvelcoll_w,vvelcoll_w]';
    
% Establishing the RHS B (body fixed)
    B(1:N_panels,1)            = -(U+rotvel_Ag)*norm_Ag';
    B(N_panels+1:2*N_panels,1) = -(U+rotvel_Bg)*norm_Bg';

% including wake induced velocity in RHS
    C = B;
    for i=1:N_panels
        [uwA_A,vwA_A]  = VORT_VEL(wakestrengthA,xcoll_Ag(i),ycoll_Ag(i),wakeglob_Ax,wakeglob_Ay);
        [uwA_B,vwA_B]  = VORT_VEL(wakestrengthA,xcoll_Ag(i),ycoll_Ag(i),wakeglob_Ax,wakeglob_Ay);
        velwnorm       = [[uwA_A;uwA_B],[vwA_A;vwA_B]]*norm_Ag';
        C(i,1)         = C(i,1)-sum(velwnorm);
    end
    for i=1:N_panels
        [uwA_A,vwA_A]  = VORT_VEL(wakestrengthB,xcoll_Bg(i),ycoll_Bg(i),wakeglob_Bx,wakeglob_By);
        [uwA_B,vwA_B]  = VORT_VEL(wakestrengthA,xcoll_Bg(i),ycoll_Bg(i),wakeglob_Ax,wakeglob_Ay);
        velwnorm       = [[uwA_A;uwA_B],[vwA_A;vwA_B]]*norm_Bg';
        C(i+N_panels,1)         = C(i+N_panels,1)-sum(velwnorm);
    end
    C(2*N_panels+1,1) = gamma_prevtA;
    C(2*N_panels+2,1) = gamma_prevtB;
    
%     solving the linear equation AX = B
    X = mldivide(A,C);

%     updating the wake locations
    wakeglob_Ax(timestep) = n_wake_Axg;
    wakeglob_Ay(timestep) = n_wake_Ayg;
    wakeglob_Bx(timestep) = n_wake_Bxg;
    wakeglob_By(timestep) = n_wake_Byg;
    
    wakestrengthA(timestep) = X(end-1);
    wakestrengthB(timestep) = X(end);
    
    uvel_w = zeros(timestep,1);
    vvel_w = zeros(timestep,1);
    for i=1:timestep
        [uAg_Aw,vAg_Aw]           = VORT_VEL(wakestrengthA,wakeglob_Ax(i),wakeglob_Ay(i),wakeglob_Ax,wakeglob_Ay);
        [uAg_Awvort,vAg_Awvort]   = VORT_VEL(X(end-1),wakeglob_Ax(i),wakeglob_Ay(i),xvort_Ag,yvort_Ag);
        [uAg_Bw,vAg_Bw]           = VORT_VEL(wakestrengthB,wakeglob_Ax(i),wakeglob_Ay(i),wakeglob_Bx,wakeglob_By);
        [uAg_Bwvort,vAg_Bwvort]   = VORT_VEL(X(end),wakeglob_Ax(i),wakeglob_Ay(i),xvort_Bg,yvort_Bg);
        
        
        uvel_w(i) = sum([uAg_Aw;uAg_Awvort;uAg_Bw;uAg_Bwvort]);
        vvel_w(i) = sum([vAg_Aw;vAg_Awvort;vAg_Bw;vAg_Bwvort]);
    end
    
     wakeglob_Ax(1:timestep-1,1) =  wakeglob_Ax(1:timestep-1,1) + (uvel_w(1:timestep-1,1)+U(1))*dt ; 
     wakeglob_Ay(1:timestep-1,1) =  wakeglob_Ay(1:timestep-1,1) + (vvel_w(1:timestep-1,1)+U(2))*dt ;
     
     wakeglob_Axt(1:timestep-1,timestep) = wakeglob_Ax(1:timestep-1,1);
     wakeglob_Ayt(1:timestep-1,timestep) = wakeglob_Ay(1:timestep-1,1);
     
     
     for i=1:timestep
        [uAg_Aw,vAg_Aw]           = VORT_VEL(wakestrengthB,wakeglob_Bx(i),wakeglob_By(i),wakeglob_Bx,wakeglob_By);
        [uAg_Awvort,vAg_Awvort]   = VORT_VEL(X(end),wakeglob_Bx(i),wakeglob_By(i),xvort_Bg,yvort_Bg);
        [uAg_Bw,vAg_Bw]           = VORT_VEL(wakestrengthA,wakeglob_Bx(i),wakeglob_By(i),wakeglob_Ax,wakeglob_Ay);
        [uAg_Bwvort,vAg_Bwvort]   = VORT_VEL(X(end-1),wakeglob_Bx(i),wakeglob_By(i),xvort_Ag,yvort_Ag);
        
        uvel_w(i) = sum([uAg_Aw;uAg_Awvort;uAg_Bw;uAg_Bwvort]);
        vvel_w(i) = sum([vAg_Aw;vAg_Awvort;vAg_Bw;vAg_Bwvort]);
    end
    
     wakeglob_Bx(1:timestep-1,1) =  wakeglob_Bx(1:timestep-1,1) + (uvel_w(1:timestep-1,1)+U(1))*dt ; 
     wakeglob_By(1:timestep-1,1) =  wakeglob_By(1:timestep-1,1) + (vvel_w(1:timestep-1,1)+U(2))*dt ;
     
     wakeglob_Bxt(1:timestep-1,timestep) = wakeglob_Bx(1:timestep-1,1);
     wakeglob_Byt(1:timestep-1,timestep) = wakeglob_By(1:timestep-1,1);
     
     gamma_prevtA = sum(X(1:N_panels));
     gamma_prevtB = sum(X(N_panels+1:2*N_panels));
     
     waitbar(timestep / steps)
     if (timestep==50)
         1==1;
     end
end
% disp("done:)")
% 
% %% plotting
% figure
% for i=2:length(t)
%     [xA,yA] = BODY_GLOB(coll_bA(:,1),coll_bA(:,2),alpha_func(t(i)),X_0A,Y_0A);
%     [xB,yB] = BODY_GLOB(coll_bB(:,1),coll_bB(:,2),alpha_func(t(i)),X_0B,Y_0B);
%     plot(xA,yA,"r","LineWidth",3)
%     hold on
%     plot(xB,yB,"r","LineWidth",3)
%     scatter(wakeglob_Axt(1:i-1,i),wakeglob_Ayt(1:i-1,i),10,'b')
%     hold on
%     scatter(wakeglob_Bxt(1:i-1,i),wakeglob_Byt(1:i-1,i),10,'b')
%     hold off
%     axis([-2,30,-7,7])
%     pause(.001)
% %     drawnow   
% end