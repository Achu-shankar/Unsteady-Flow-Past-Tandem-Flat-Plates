function X_g = BODY_GLOB(X_b,alpha,X_0)
    X_g = X_b*[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)] + X_0;
 end 