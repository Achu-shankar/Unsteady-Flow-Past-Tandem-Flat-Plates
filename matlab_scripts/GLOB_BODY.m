function X_b = GLOB_BODY(X_g,alpha,X_0)
    X_b = (X_g-X_0)*[cos(alpha) sin(alpha);-sin(alpha) cos(alpha)] ; 
end
