%% Self-implemented optimisation
% We see quite a non-smooth behaviour, so idea:
% First make a smooth response surf using LHS+least squares and perform our
% own optimisation on that surface
clear all
close all

%% Steepest decent method for unconstrained problem
% Set avg parameters

t=3e-3;
r=0.4; % r   *c
h=0.25; % h   *c
c=0.012;

% Reading data
N=200;
Fobj=dlmread(['Own_Optimisation/Fobj_',mat2str(N)]);
C=dlmread(['Own_Optimisation/C_',mat2str(N)]); % We also have to write+read these otherwise these change
R=dlmread(['Own_Optimisation/R_',mat2str(N)]); 

% Produce fit
sf = fit([R, C],Fobj','poly22');
x_start=[0.35 ; 0.02];

% Set optimisation properties
conv_crit_SD=1e-6;  % How many steepest decent iterations?
N_fibo=50;      
line_factor=0.8;

[r,c,f,N_iter_SD]=Steepest_Decent(sf,x_start,conv_crit_SD,N_fibo,line_factor);

plot(sf,[R,C],Fobj')
hold on
scatter3(r(1),c(1),sf([r(1) c(1)]),'r','filled')
for i=[1:N_iter_SD]
    plot1=plot3([r(i); r(i+1)],[c(i); c(i+1)],[sf([r(i) c(i)]); sf([r(i+1) c(i+1)])],'r','LineWidth',2);
    M1="Convergence Path";
    scatter3(r(i+1),c(i+1),sf([r(i+1) c(i+1)]),'r','filled')
end

plot2=scatter3(0.4304,   0.0313,  6.0954,'cd','filled');
M2="Exact Solution";
legend([plot1, plot2],M1, M2)

xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')

figure()
plot([1:N_iter_SD+1],f)
hold on
plot([0 12],[6.0954 6.0954],'r--')
legend('Iterative solution','Exact solution')
xlabel('$Number \ of \ iterations$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$f$ [N]','fontsize',16,'Interpreter','LaTex')



