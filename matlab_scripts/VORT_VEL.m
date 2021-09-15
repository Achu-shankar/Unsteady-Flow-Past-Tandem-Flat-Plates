function [u,w]= VORT_VEL(gamma_i,x,y,x_i,y_i)
        sigma = 3/50^2;
        r = sqrt((x-x_i).^2+(y-y_i).^2);
        u = zeros(length(r),1);
        w = zeros(length(r),1);
        
        const = gamma_i./(2*pi*r);
        a = const.*(y-y_i) ;
        b = const.*(x-x_i) ;
        
        cond1 = r >=sigma;
        u(cond1) =  a(cond1)./r(cond1);
        w(cond1) = -b(cond1)./r(cond1);
        
        cond2 = r <sigma;
        u(cond2) =  a(cond2)./sigma;
        w(cond2) = -b(cond2)./sigma;
        
        u(isnan(u)) = 0;
        w(isnan(w)) = 0;        
        
end