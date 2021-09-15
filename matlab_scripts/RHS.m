B(1:N_panels)            = -(U-rotv_gA)*norm_gA';
B(N_panels+1:2*N_panels) = -(U-rotv_gB)*norm_gB';

for i =1:N_panels
    [u_g,v_g]   = VORT_VEL(w_str_A,coll_A(i,1),coll_A(i,2),wake_gA(:,1),wake_gA(:,2));
    [u_g1,v_g1] = VORT_VEL(w_str_B,coll_A(i,1),coll_A(i,2),wake_gB(:,1),wake_gB(:,2));
    vel           = [[u_g;u_g1],[v_g;v_g1]]*norm_gA';
    B(i) = B(i) - sum(vel);
end
for i =1:N_panels
    [u_g,v_g]     = VORT_VEL(w_str_B,coll_B(i,1),coll_B(i,2),wake_gB(:,1),wake_gB(:,2));
    [u_g1,v_g1]   = VORT_VEL(w_str_A,coll_B(i,1),coll_B(i,2),wake_gA(:,1),wake_gA(:,2));
    vel         = [[u_g;u_g1],[v_g;v_g1]]*norm_gB';
    B(i+N_panels) = B(i+N_panels) - sum(vel);
end
B(2*N_panels+1) =gamma_pre_A;
B(2*N_panels+2) =gamma_pre_B;
