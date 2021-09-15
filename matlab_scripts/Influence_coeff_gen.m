%F1
for i=1:N_panels
    [u_g,v_g] = VORT_VEL(1,coll_gA(i,1),coll_gA(i,2),vort_gA(:,1),vort_gA(:,2));
    A(i,1:N_panels) = norm_gA*[u_g,v_g]';
end
%F2
for i=1:N_panels
    [u_g,v_g] = VORT_VEL(1,coll_gA(i,1),coll_gA(i,2),vort_gB(:,1),vort_gB(:,2));
    A(i,N_panels+1:2*N_panels) = norm_gA*[u_g,v_g]';
end

%F3
[u_g,v_g] = VORT_VEL(1,coll_gA(:,1),coll_gA(:,2),n_wake_gA(1),n_wake_gA(2));
A(1:N_panels,2*N_panels+1) = norm_gA*[u_g,v_g]';

%F4
[u_g,v_g] = VORT_VEL(1,coll_gA(:,1),coll_gA(:,2),n_wake_gB(1),n_wake_gB(2));
A(1:N_panels,2*N_panels+2) = norm_gA*[u_g,v_g]';

% G1
for i=1:N_panels
    [u_g,v_g] = VORT_VEL(1,coll_B(i,1),coll_B(i,2),vort_A(:,1),vort_A(:,2));
    A(i+N_panels,1:N_panels) = norm_gB*[u_g,v_g]';
end

%G2
for i=1:N_panels
    [u_g,v_g] = VORT_VEL(1,coll_B(i,1),coll_B(i,2),vort_B(:,1),vort_B(:,2));
    A(i+N_panels,N_panels+1:2*N_panels) = norm_gB*[u_g,v_g]';
end

%G3
[u_g,v_g] = VORT_VEL(1,coll_gB(:,1),coll_gB(:,2),n_wake_gA(1),n_wake_gA(2));
A(N_panels+1:2*N_panels,2*N_panels+1) = norm_gB*[u_g,v_g]';

%G4
[u_g,v_g] = VORT_VEL(1,coll_gB(:,1),coll_gB(:,2),n_wake_gB(1),n_wake_gB(2));
A(N_panels+1:2*N_panels,2*N_panels+2) = norm_gB*[u_g,v_g]';

