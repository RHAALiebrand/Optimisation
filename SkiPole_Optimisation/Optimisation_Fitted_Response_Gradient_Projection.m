%% Self-implemented optimisation
% We see quite a non-smooth behaviour, so idea:
% First make a smooth response surf using LHS+least squares and perform our
% own optimisation on that surface
clear all
close all

%% Initialisation
ndv = 2;
ncons = 2;
dx = 1.0e-8; % Finite difference step

N=200;
R_min=0.3;
R_max=0.5;
C_min=5e-3;
C_max=4e-2;
xmin = [R_min, C_min]';
xmax = [R_max, C_max]';
h=0.25;
t=3e-3;
F=import_database();
N_its=25;  % How many total iterations?
N_fibo=100; 
N_bnd = 10;
line_factor=0.1;

% Define starting point
x0 = [0.48,0.015]'; % Assumed to be feasible
xtrace = x0;

%% Model import
% LHS design or data read
[Fobj,G1,G2,G3,C,R] = readData(N);
% [Fobj,G1,G2,G3,C,R] = createLHS(N,R_min,R_max,C_min,C_max,t,h)

% Remove outliers
G1 = filloutliers(G1,'linear');
G2 = filloutliers(G2,'linear');
Fobj = filloutliers(Fobj,'linear');

% Fit surface to the data and its derivatives
[sf , gof_f] = fit([R,C],Fobj','poly22');
[sg1, gof_g1] = fit([R,C],G1','poly22');
[sg2, gof_g2] = fit([R,C],G2','poly22');

%% Optimise
[xtrace,active,f,g1,g2,dfdx,dgdx]=Steepest_Decent_Bound(sf,x0,N_its,N_fibo,N_bnd,line_factor,sg1,sg2,xmin,xmax);

%% Check optimality
x = xtrace(:,end);
Check_Optimality(x,g1,g2,dfdx,dgdx,sf,sg1,sg2,dx,active)

plot(R_min,R_max,C_min,C_max,sf,sg1,sg2,xtrace)

function [Fobj,G1,G2,G3,C,R] = readData(N)
    % Reading the above iso looping
    Fobj=dlmread(['Own_Optimisation/Fobj_',mat2str(N)]);
    G1=dlmread(['Own_Optimisation/G1_',mat2str(N)]);
    G2=dlmread(['Own_Optimisation/G2_',mat2str(N)]);
    G3=dlmread(['Own_Optimisation/G3_',mat2str(N)]);
    C=dlmread(['Own_Optimisation/C_',mat2str(N)]); % We also have to write+read these otherwise these change
    R=dlmread(['Own_Optimisation/R_',mat2str(N)]); 
end

function [Fobj,G1,G2,G3,C,R] = createLHS(N,R_min,R_max,C_min,C_max,t,h)
    LHS=lhsdesign(N,2);

    R=LHS(:,1)*(R_max-R_min)+R_min;
    C=LHS(:,2)*(C_max-C_min)+C_min;

    Fobj=zeros(1,N);
    G1=zeros(1,N);
    G2=zeros(1,N);
    G3=zeros(1,N);

    for i=[1:length(R)]
        i
        Fobj(i)=fobj(R(i),h,C(i),t,F);
        g=g_i(R(i),h,C(i),t);
        G1(i)=g(1);
        G2(i)=g(2);
        G3(i)=g(3);
    end

    dlmwrite(['Own_Optimisation/Fobj_',mat2str(N)],Fobj)
    dlmwrite(['Own_Optimisation/G1_',mat2str(N)],G1)
    dlmwrite(['Own_Optimisation/G2_',mat2str(N)],G2)
    dlmwrite(['Own_Optimisation/G3_',mat2str(N)],G3)
    dlmwrite(['Own_Optimisation/R_',mat2str(N)],R)
    dlmwrite(['Own_Optimisation/C_',mat2str(N)],C)
end

function []= plot(R_min,R_max,C_min,C_max,sf,sg1,sg2,xtrace)
Nfine = 100;
Rfine=linspace(R_min,R_max,Nfine);
Cfine=linspace(C_min,C_max,Nfine);
[Rfine,Cfine]=meshgrid(Rfine,Cfine);

FobjSF = sf(Rfine,Cfine);
G1f = sg1(Rfine,Cfine);
G2f = sg2(Rfine,Cfine);

figure()
hold on
% contour(R,C,Fobj,50,'ShowText','on')
contour(Rfine,Cfine,FobjSF)
% contour(Rfine,Cfine,G2f,50,'ShowText','on')
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
hold on 
contour(Rfine,Cfine,G1f,[0.0 0.0],'b')
contour(Rfine,Cfine,G2f,[0.0 0.0],'r')
contour(Rfine,Cfine,G1f,[0.1 0.1],'--b')   
contour(Rfine,Cfine,G2f,[0.1 0.1],'--r')  
scatter(0.4422,0.0189,'cd','filled')
scatter(xtrace(1,:),xtrace(2,:)','k')
line(xtrace(1,:),xtrace(2,:))
%contour(R_mesh,C_mesh,G3,[0.0 0.0],'c')  INACTIVE
legendInfo={'Objective function','Constrain 1','Constrain 2','Active Const 1','Active Const 1','Estimated Global Minimum','Convergence Path'};
% start=length(legendInfo);
% for i=[1:maxiter] 
%     hold on
%     scatter(r_iter(i),c_iter(i))
%     legendInfo{start+i}=['Iter ',mat2str(i-1)];
% end
% legend(legendInfo)
caxis([min(FobjSF(:)) max(FobjSF(:))])
end






