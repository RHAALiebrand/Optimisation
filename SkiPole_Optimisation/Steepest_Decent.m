%% Steepest decent method
% We use fibonacci to perform the line searches 
function [r,c,f,N_iter]=Steepest_Decent(sf,x_start,conf_crit_SD,N_fibo,line_factor)
    % Set outputs
    r(1)=x_start(1);
    c(1)=x_start(2);
    f(1)=sf([r(1) c(1)]);
    % Looping process
    iterate=true;
    N_iter=1;
    while iterate==true 
        if N_iter~=1
            x0=[r_new;c_new];
        else 
            x0=x_start;
        end 
        % Directions
        [df_dx,df_dy]=differentiate(sf,x0(1),x0(2));
        d=-[df_dx; df_dy];
        d_unit=d/norm(d);
        x1=x0+d_unit.*x_start*line_factor; 
        
        % Fibo search
        [r_new,c_new]=Fibonacci(x0,x1,N_fibo,sf);
        r(N_iter+1)=r_new;
        c(N_iter+1)=c_new;
        f(N_iter+1)=sf([r_new c_new]);
        if abs(f(N_iter+1)-f(N_iter))<conf_crit_SD
            iterate=false; 
        end 
        x0_old=x0;
        N_iter=N_iter+1;
    end 
    N_iter=N_iter-1;
end 


