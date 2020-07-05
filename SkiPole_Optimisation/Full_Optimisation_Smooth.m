%% Full optimisation SMOOTH
% We see quite a non-smooth behaviour, so idea:
% First make a smooth response surf using LHS+least squares and perform our
% own optimisation on that surface
clear all
close all

F=import_database();
%% Produce LHS space
N=2500;
R_min=0.3;
R_max=0.5;
H_min=0.2;
H_max=0.6;
C_min=5e-3;
C_max=4e-2;
t_min=0.8e-3;
t_max=2e-3;

R=linspace(R_min,R_max,N);
C=linspace(C_min,C_max,N);
[R_mesh,C_mesh]=meshgrid(R,C);

% LHS=lhsdesign(N,4);
% 
% R=LHS(:,1)*(R_max-R_min)+R_min;
% C=LHS(:,2)*(C_max-C_min)+C_min;
% H=LHS(:,3)*(H_max-H_min)+H_min;
% t=LHS(:,4)*(t_max-t_min)+t_min;
% for i=[1:N]
%         i
%         Fobj(i)=fobj(R(i),H(i),C(i),t(i),F);
%         g=g_i(R(i),H(i),C(i),t(i));
%         G1(i)=g(1);
%         G2(i)=g(2);
%         G3(i)=g(3);
% end

%% Reading the above iso looping
Fobj=dlmread(['Full_Optimisation/Fobj_',mat2str(N)]);
G1=dlmread(['Full_Optimisation/G1_',mat2str(N)]);
G2=dlmread(['Full_Optimisation/G2_',mat2str(N)]);
G3=dlmread(['Full_Optimisation/G3_',mat2str(N)]);
C=dlmread(['Full_Optimisation/C_',mat2str(N)]); % We also have to write+read these otherwise these change
R=dlmread(['Full_Optimisation/R_',mat2str(N)]); 
H=dlmread(['Full_Optimisation/H_',mat2str(N)]); 
t=dlmread(['Full_Optimisation/t_',mat2str(N)]); 

% dlmwrite(['Full_Optimisation/Fobj_',mat2str(N)],Fobj)
% dlmwrite(['Full_Optimisation/G1_',mat2str(N)],G1)
% dlmwrite(['Full_Optimisation/G2_',mat2str(N)],G2)
% dlmwrite(['Full_Optimisation/G3_',mat2str(N)],G3)
% dlmwrite(['Full_Optimisation/R_',mat2str(N)],R)
% dlmwrite(['Full_Optimisation/C_',mat2str(N)],C)
% dlmwrite(['Full_Optimisation/H_',mat2str(N)],H)
% dlmwrite(['Full_Optimisation/t_',mat2str(N)],t)
%% Fit functions
% Const_1=find(G1>0);sf
% Const_2=find(G2>0);
sf = MultiPolyRegress([R,H,C,t],Fobj',2);
sg1 = MultiPolyRegress([R,H,C,t],G1',2);
G2 = filloutliers(G2,'linear');
sg2 = MultiPolyRegress([R,H,C,t],G2',2);

bound_margin=0.9;
c_lower=C_min*(2-bound_margin);
c_upper=C_max*bound_margin;
rc_lower=R_min*(2-bound_margin);
rc_upper=R_max*bound_margin;
h_lower=H_min*(2-bound_margin);
h_upper=H_max*bound_margin;
t_lower=t_min*(2-bound_margin);
t_upper=t_max*bound_margin;
%% First order optimisation process
func=@(x) changed_nomenclature_fit(x,sf.PolynomialExpression);
cons=@(x) g_own_full(x,sg1.PolynomialExpression,sg2.PolynomialExpression);
x0=[0.36 0.5 0.02 1e-3];
maxiter=1000;
% SQP
options = optimoptions(@fmincon,'Display','Iter','MaxIter',maxiter,'Algorithm','sqp');
[x_sqp, fval_sqp, exitflag_sqp, output_sqp, lambda_sqp]=fmincon(func,x0,[],[],[],[],[rc_lower h_lower c_lower t_lower],[rc_upper h_upper c_upper t_upper],cons,options);

% IP
options = optimoptions(@fmincon,'Display','Iter','MaxIter',maxiter,'Algorithm','interior-point');
[x_ip, fval_ip, exitflag_ip, output_ip, lambda_ip]=fmincon(func,x0,[],[],[],[],[rc_lower h_lower c_lower t_lower],[rc_upper h_upper c_upper t_upper],cons,options);

% AS
options = optimoptions(@fmincon,'Display','Iter','MaxIter',maxiter,'Algorithm','active-set');
[x_as, fval_as, exitflag_as, output_as, lambda_as]=fmincon(func,x0,[],[],[],[],[rc_lower h_lower c_lower t_lower],[rc_upper h_upper c_upper t_upper],cons,options);

first_order=[fval_sqp,fval_ip,fval_as]

%% Zeroth order methods

% % Nelder-mead 
% options_nelder = optimset('Display','Iter');
% [x_nelder,fval_nelder,exitflag_nelder,output_nelder]=fminsearchcon(func,x0,[rc_lower h_lower c_lower t_lower],[rc_upper h_upper c_upper t_upper],[],[],cons)%,options);
% 
% % Genetic Algorithm
% rng default
% options_ga=optimoptions('ga','Display','Iter')
% [x_ga,fval_ga,exitflag_ga,output_ga]=ga(func,4,[],[],[],[],[rc_lower h_lower c_lower t_lower],[rc_upper h_upper c_upper t_upper],cons,options_ga);
% 
% % We can not really build a convergence path for this algo
% addpath('./Constrained_PS/');
% %options_ps = optimoptions('Display','Iter');
% [x_ps,fval_ps,exitflag_ps,output_ps,population_ps]=pso(func,4,[],[],[],[],[rc_lower h_lower c_lower t_lower],[rc_upper h_upper c_upper t_upper],cons)%,options_ps);
% iter_ps=size(population_ps(:,1));

% Read results iso running algo's all the time
% dlmwrite('Full_Optimisation\Full_Results_0thorder_2500',[x_nelder,fval_nelder,output_nelder.iterations,output_nelder.funcCount;x_ga,fval_ga,output_ga.generations,output_ga.funccount;x_ps,fval_ps,output_ps.generations,output_ps.generations*iter_ps(1)])
zeroth_results=dlmread('Full_Optimisation\Full_Results_0thorder_2500');
% Nomenclature r,h,c,t,fval,iter,funccount

%% Check evaluation speed
% N=100000
% tic
% for i=[1:N]
%     func(x_as)
% end 
% t=toc
% t_eval=t/N


function fval=changed_nomenclature_fit(x,polyexp)
    fval=polyexp(x(1),x(2),x(3),x(4));
end 


