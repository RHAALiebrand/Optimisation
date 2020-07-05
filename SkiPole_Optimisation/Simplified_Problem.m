%% Simplified problem
% We only consider r/c and c so we exclude h/c and t for the reasons
% presented in the report
close all 
clear all 

%% Set up problem
N=30;
R=linspace(0.3,0.5,N);
C=linspace(5e-3,4e-2,N);
Fobj=zeros(length(C),length(R));
G1=zeros(length(C),length(R));
G2=zeros(length(C),length(R));
G3=zeros(length(C),length(R));
[R_mesh,C_mesh]=meshgrid(R,C);
% for i=[1:length(R)]
%     i
%     for j=[1:length(C)]
%         j
%         Fobj(j,i)=fobj(R(i),h,C(j),t,F);
%         g=g_i(R(i),h,C(j),t);
%         G1(j,i)=g(1);
%         G2(j,i)=g(2);
%         G3(j,i)=g(3);
%     end
% end

%% Read instead of building
Fobj=dlmread('Simplified_Problem/Fobj_30.dat');
G1=dlmread('Simplified_Problem/G1_30.dat');
G2=dlmread('Simplified_Problem/G2_30.dat');
G3=dlmread('Simplified_Problem/G3_30.dat');

bound_margin=0.9;
c_lower=0.005*(2-bound_margin);
c_upper=0.04*bound_margin;
rc_lower=0.3*(2-bound_margin);
rc_upper=0.5*bound_margin;

%% First order optimisation
% Optimisation
options = optimoptions(@fmincon,'Display','Iter','OutputFcn',@myoutput_simplified,'Algorithm','sqp');
func=@(x) fobj2(x);
cons=@(x) g_i2(x);
x0=[0.40 0.02];
%[x, fval, exitflag, output, lambda]=fmincon(func,x0,[],[],[],[],[rc_lower c_lower],[rc_upper c_upper],cons,options);
%maxiter=output.iterations; %Overwriting if converged

% Read results
r_iter=[];
c_iter=[];
f_iter=[];
maxiter=4;
for i=[0:maxiter] % Reading starts at 0
    results=dlmread(['Simplified_Problem/Iterative_Results/results_iteration_sqp_',mat2str(i)]);
    r_iter=[r_iter;results(1)];
    c_iter=[c_iter;results(2)];
    f_iter=[f_iter;results(3)];
end


%% Plotting
% 2-D Plot
figure()
hold on
contour(R_mesh,C_mesh,Fobj,30,'ShowText','on')
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
hold on 
contour(R_mesh,C_mesh,G1,[0.0 0.0],'b')
contour(R_mesh,C_mesh,G2,[0.0 0.0],'r')
contour(R_mesh,C_mesh,G1,[0.1 0.1],'--b')   
contour(R_mesh,C_mesh,G2,[0.1 0.1],'--r')  
scatter(0.463,0.0324,'cd','filled')
plot(r_iter,c_iter,'m')
%contour(R_mesh,C_mesh,G3,[0.0 0.0],'c')  INACTIVE
legendInfo={'Objective function','Constrain 1','Constrain 2','Active Const 1','Active Const 1','Estimated Global Minimum','Convergence Path'};
start=length(legendInfo);
for i=[1:maxiter] 
    hold on
    scatter(r_iter(i),c_iter(i))
    legendInfo{start+i}=['Iter ',mat2str(i-1)];
end
legend(legendInfo)
caxis([min(Fobj(:)) max(Fobj(:))])

% 3-D Plot
figure
surf(R_mesh,C_mesh,Fobj)
hold on
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
scatter3(0.463,0.0324,6.1,'cd','filled')
plot3(r_iter,c_iter,f_iter,'r')
legendInfo={'Objective function','Estimated Global Minimum','Convergence Path'};
start=length(legendInfo);
for i=[1:maxiter+1] 
    hold on
    scatter3(r_iter(i),c_iter(i),f_iter(i))
    legendInfo{start+i}=['Iter ',mat2str(i-1)];
end
legend(legendInfo)

%% Zeroth order optimization
% Particle Swarm 
max_iter_ps=20;
options_ps = optimoptions('particleswarm','Display','iter','MaxIterations',max_iter_ps,'SwarmSize',40,'OutputFcn',@pswplotranges);
[x_ps,fval_ps,exitflag_ps,output_ps]=particleswarm(func,2,[0.34 0.0122],[0.48 0.0325],options_ps);
maxiter=output_ps.iterations;

% Nelder-mead
max_iter_nelder=200;
options_nelder = optimset('Display','Iter','OutputFcn',@myoutput_nelder,'MaxIter',max_iter_nelder);
[x_nelder,fval_nelder,exitflag_nelder,output_nelder]=fminsearch(func,x0,options_nelder);

% GA
maxgen=10;
options_ga = optimoptions('ga','Display','iter','MaxGenerations',maxgen,'PopulationSize',10)%,'OutputFcn',@myoutputfcn_ga);;

rng default
[x_ga,fval_ga,exitflag_ga,output_ga]=ga(func,2,[],[],[],[],[0.34 0.0122],[0.48 0.0325],[],options_ga);

figure()
hold on
contour(R_mesh,C_mesh,Fobj,50,'ShowText','on')
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
hold on 
contour(R_mesh,C_mesh,G1,[0.0 0.0],'b')
contour(R_mesh,C_mesh,G2,[0.0 0.0],'r')
contour(R_mesh,C_mesh,G1,[0.1 0.1],'--b')   
contour(R_mesh,C_mesh,G2,[0.1 0.1],'--r')  
scatter(0.463,0.0324,'cd','filled')
%contour(R_mesh,C_mesh,G3,[0.0 0.0],'c')  INACTIVE
legendInfo={'Objective function','Constrain 1','Constrain 2','Active Const 1','Active Const 1','Estimated Global Minimum','Convergence Path Nelder','Found Optimum Nelder','Convergence Path Swarm','Found Optimum Swarm','Found Optimum GA'};
start=length(legendInfo);

for i=[0:output_nelder.iterations] 
    data=dlmread(['Simplified_Problem\nelder\',num2str(i)]);
    r_nelder(i+1)=data(1);
    c_nelder(i+1)=data(2);
end
plot(r_nelder,c_nelder,'g')
scatter(r_nelder(end),c_nelder(end),'gd','filled')

for i=[0:output_ps.iterations] 
    data=dlmread(['Simplified_Problem\ps\',num2str(i)]);
    r_ps(i+1)=data(1);
    c_ps(i+1)=data(2);
end
plot(r_ps,c_ps,'m')
scatter(r_ps(4),c_ps(4),'md','filled')
scatter(x_ga(1),x_ga(2),'rd','filled')
legend(legendInfo)
caxis([min(Fobj(:)) max(Fobj(:))])



function stop = pswplotranges(optimValues,state)
optimValues;
coor=reshape(optimValues.swarm,[1,numel(optimValues.swarm)]);
value=reshape(optimValues.swarmfvals,[1,numel(optimValues.swarmfvals)]);
results=[optimValues.bestx,optimValues.bestfval,coor,value];
dlmwrite(['Simplified_Problem/ps/x0_40_',mat2str(optimValues.iteration)],results)
stop = false; % This function does not stop the solver
switch state
    case 'init'
        nplot = size(optimValues.swarm,2); % Number of dimensions
        for i = 1:nplot % Set up axes for plot
            subplot(nplot,1,i);
            tag = sprintf('psoplotrange_var_%g',i); % Set a tag for the subplot
            semilogy(optimValues.iteration,0,'-k','Tag',tag); % Log-scaled plot
            ylabel(num2str(i))
        end
        xlabel('Iteration','interp','none'); % Iteration number at the bottom
        subplot(nplot,1,1) % Title at the top
        title('Log range of particles by component')
        setappdata(gcf,'t0',tic); % Set up a timer to plot only when needed
    case 'iter'
        nplot = size(optimValues.swarm,2); % Number of dimensions
        for i = 1:nplot
            subplot(nplot,1,i);
            % Calculate the range of the particles at dimension i
            irange = max(optimValues.swarm(:,i)) - min(optimValues.swarm(:,i));
            tag = sprintf('psoplotrange_var_%g',i);
            plotHandle = findobj(get(gca,'Children'),'Tag',tag); % Get the subplot
            xdata = plotHandle.XData; % Get the X data from the plot
            newX = [xdata optimValues.iteration]; % Add the new iteration
            plotHandle.XData = newX; % Put the X data into the plot
            ydata = plotHandle.YData; % Get the Y data from the plot
            newY = [ydata irange]; % Add the new value
            plotHandle.YData = newY; % Put the Y data into the plot
        end
        if toc(getappdata(gcf,'t0')) > 1/30 % If 1/30 s has passed
          drawnow % Show the plot
          setappdata(gcf,'t0',tic); % Reset the timer
        end
    case 'done'
        % No cleanup necessary
end
end 

function stop = myoutput_nelder(x,optimvalues,state);
 
    stop = false;
    if isequal(state,'iter')
      results=[x,optimvalues.fval];
      dlmwrite(['Simplified_Problem\nelder\x0_40_',mat2str(optimvalues.iteration)],results)
    end
end

function [state, options,optchanged] = myoutputfcn_ga(options,state,flag)
%GAOUTPUTFCNTEMPLATE Template to write custom OutputFcn for GA.
% [STATE, OPTIONS, OPTCHANGED] = GAOUTPUTFCNTEMPLATE(OPTIONS,STATE,FLAG)
% where OPTIONS is an options structure used by GA.
%
% STATE: A structure containing the following information about the state
% of the optimization:
% Population: Population in the current generation
% Score: Scores of the current population
% Generation: Current generation number
% StartTime: Time when GA started
% StopFlag: String containing the reason for stopping
% Selection: Indices of individuals selected for elite,
% crossover and mutation
% Expectation: Expectation for selection of individuals
% Best: Vector containing the best score in each generation
% LastImprovement: Generation at which the last improvement in
% fitness value occurred
% LastImprovementTime: Time at which last improvement occurred
%
% FLAG: Current state in which OutputFcn is called. Possible values are:
% init: initialization state
% iter: iteration state
% interrupt: intermediate state
% done: final state
%
% STATE: Structure containing information about the state of the
% optimization.
%
% OPTCHANGED: Boolean indicating if the options have changed.
%
%	See also PATTERNSEARCH, GA, GAOPTIMSET
% Copyright 2004-2006 The MathWorks, Inc.
% $Revision: 1.1.6.5 $ $Date: 2007/08/03 21:23:22 $
optchanged = false;
switch flag
case 'init'
disp('Starting the algorithm');
case {'iter','interrupt'}
disp('Iterating ...')
fname=[pwd,'\Simplified_Problem\ga\',num2str(state.Generation)];
save(fname,'state')
case 'done'
disp('Performing final task');
fname=[pwd,'\',num2str(state.Generation),'.mat'];
save(fname,'state')
end
end 