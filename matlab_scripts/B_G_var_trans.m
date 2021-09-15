%body to global transformations
norm_gA = BODY_GLOB(norm_A,alpha_func_A(tstep),[0,0]);
vort_gA = BODY_GLOB(vort_A,alpha_func_A(tstep),X_0A);
coll_gA = BODY_GLOB(coll_A,alpha_func_A(tstep),X_0A);
n_wake_gA = BODY_GLOB(n_wake_A,alpha_func_A(tstep),X_0A);
rotv_gA = BODY_GLOB([zeros(length(coll_A(:,1)),1),omega_alpha_A(tstep)...
          *(coll_A(:,1)-P_rotA(1))],alpha_func_A(t(tstep)),[0,0]);

norm_gB = BODY_GLOB(norm_B,alpha_func_B(tstep),[0,0]);
vort_gB = BODY_GLOB(vort_B,alpha_func_B(tstep),X_0B);
coll_gB = BODY_GLOB(coll_B,alpha_func_B(tstep),X_0B);
n_wake_gB = BODY_GLOB(n_wake_B,alpha_func_B(tstep),X_0B);
rotv_gB = BODY_GLOB([zeros(length(coll_B(:,1)),1),omega_alpha_B(tstep)...
         *(coll_A(:,1)-P_rotB(1))],alpha_func_B(t(tstep)),[0,0]);