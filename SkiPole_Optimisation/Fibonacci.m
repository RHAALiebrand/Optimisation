%% Finbonacci line search
% Inputs
%      - a,b for bracketing
%      - func 
%      - search direction
%      - Number of iteration

function [r_new,c_new]=Fibonacci(x0,x1,N_iter,sf)
    %% Set-up the search direction so C=aR+b so that we search in R
    RC=(x1(2)-x0(2))/(x1(1)-x0(1));
    constant=x0(2)-RC*x0(1);

    %% Bracketing        
    
    a=x0(1);
    b=x1(1);
    
    if x0(1) > x1(1)
        b = x0(1);
        a = x1(1);
    end

    %% Setup Fibonacci reeks
    fibo_old=1;
    fibo_new=1;
    for i=1:N_iter
        if i==1 || i==2
        fibo(i)=1;
       continue;
        end
        fibo(i)=fibo_old+fibo_new;
        fibo_old=fibo_new;
        fibo_new=fibo(i);
    end

    % Now determine I2
    I2=(b-a)*fibo(N_iter-1)/fibo(N_iter);

    %% Start sectioning
    i=2; % Since we already did the bracktering
    while i<N_iter
        I1=(b-a);
        if I2>0.5*I1
            a_new=b-I2;
            b_new=a+I2;
        else if I2 <= 0.5*I1
                a_new=a+I2;
                b_new=b-I2;
        end
        end
        fun_eval1=sf([a_new (RC*a_new+constant)]);
        fun_eval2=sf([b_new (RC*b_new+constant)]);

        if fun_eval2>fun_eval1
            b=b_new;
            I2=fibo(N_iter-i)*I1/fibo(N_iter-i+2);
        else if fun_eval2<fun_eval1
            a=a_new;
            I2=fibo(N_iter-i)*I1/fibo(N_iter-(i-2));
        else if fun_eval2==fun_eval1
            b=b_new;
            I2=fibo(N_iter-i)*[b-a]/fibo(N_iter-(i-2));
            i=i+1;
        end
        end
        end
        i=i+1;
    end
    
r_new=0.5*(a+b);
c_new=0.5*(RC*a+constant+RC*b+constant);
end  
