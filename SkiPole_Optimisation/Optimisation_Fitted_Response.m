%% Self-implemented optimisation
% We see quite a non-smooth behaviour, so idea:
% First make a smooth response surf using LHS+least squares and perform our
% own optimisation on that surface
clear all
close all

%% Produce LHS space
N=200;
R_min=0.3;
R_max=0.5;
C_min=5e-3;
C_max=4e-2;

R=linspace(R_min,R_max,N);
C=linspace(C_min,C_max,N);
[R_mesh,C_mesh]=meshgrid(R,C);

% LHS=lhsdesign(N,2);

% R=LHS(:,1)*(R_max-R_min)+R_min;
% C=LHS(:,2)*(C_max-C_min)+C_min;

% for i=[1:N]
%         Fobj(i)=fobj(R(i),h,C(i),t,F);
%         g=g_i(R(i),h,C(i),t);
%         G1(i)=g(1);
%         G2(i)=g(2);
%         G3(i)=g(3);
% end

%% Reading the above iso looping
Fobj=dlmread(['Own_Optimisation/Fobj_',mat2str(N)]);
G1=dlmread(['Own_Optimisation/G1_',mat2str(N)]);
G2=dlmread(['Own_Optimisation/G2_',mat2str(N)]);
G3=dlmread(['Own_Optimisation/G3_',mat2str(N)]);
C=dlmread(['Own_Optimisation/C_',mat2str(N)]); % We also have to write+read these otherwise these change
R=dlmread(['Own_Optimisation/R_',mat2str(N)]); 

% dlmwrite(['Own_Optimisation/Fobj_',mat2str(N)],Fobj)
% dlmwrite(['Own_Optimisation/G1_',mat2str(N)],G1)
% dlmwrite(['Own_Optimisation/G2_',mat2str(N)],G2)
% dlmwrite(['Own_Optimisation/G3_',mat2str(N)],G3)
% dlmwrite(['Own_Optimisation/R_',mat2str(N)],R)
% dlmwrite(['Own_Optimisation/C_',mat2str(N)],C)

Const_1=find(G1>0);
Const_2=find(G2>0);
G2 = filloutliers(Fobj,'linear');
[sf , gof_f] = fit([R, C],Fobj','poly22');
[sg1, gof_g1] = fit([R, C],G1','poly22');
G2 = filloutliers(G2,'linear');
[sg2, gof_g2] = fit([R, C],G2','poly22');

bound_margin=0.9;
c_lower=0.005*(2-bound_margin);
c_upper=0.04*bound_margin;
rc_lower=0.3*(2-bound_margin);
rc_upper=0.5*bound_margin;

%% First order optimization
maxiter=100;
options = optimoptions(@fmincon,'Display','Iter','MaxIter',maxiter,'OutputFcn',@myoutput_own,'Algorithm','sqp'); %myoutput is written to extract info during process, ugly but the only option after 2 hours of work
func=@(x) sf(x);

cons=@(x) g_own(x,sg1,sg2);
x0=[0.36 0.02];
[x, fval, exitflag, output, lambda]=fmincon(func,x0,[],[],[],[],[rc_lower c_lower],[rc_upper c_upper],cons,options);

% Read results
maxiter=output.iterations; %Overwriting if converged
r_iter=[];
c_iter=[];
f_iter=[];
for i=[0:maxiter-1] % Reading starts at 0
    results=dlmread(['Own_Optimisation/Iterative_Results/results_iteration_sqp',mat2str(i)]);
    r_iter=[r_iter;results(1)];
    c_iter=[c_iter;results(2)];
    f_iter=[f_iter;results(3)];
end

figure()
plot(sf,[R,C],Fobj')
hold on
scatter3(R([Const_1,Const_2]),C([Const_1,Const_2]),Fobj([Const_1,Const_2]),'r','filled')
hold on
plot3(r_iter,c_iter,f_iter,'r')
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
zlabel('$f(r/c,c)$','fontsize',16,'Interpreter','LaTex')
legendInfo={'Feasible LHS Point','Unfeasible (Constrained) LHS Point'};
start=length(legendInfo);
for i=[1:maxiter] 
    hold on
    scatter3(r_iter(i),c_iter(i),f_iter(i),'filled')
    legendInfo{start+i}=['Iter ',mat2str(i-1)];
end
hold on
legend(legendInfo)

%% Verification that we found the global minimum
A=[2*sf.p20 0.5*sf.p11 ; 0.5*sf.p11 2*sf.p02];
b=[sf.p10 ; sf.p01];
x=[c_iter(end); r_iter(end)];
x_star=-inv(A)*b;
grad_f=A*x+b;

%% Zeroth order methods
% Nelder-mead 
options_nelder = optimset('Display','Iter','OutputFcn',@myoutput_nelder);
[x_nelder,fval_nelder,exitflag_nelder,output_nelder]=fminsearch(func,x0,options_nelder);
x_nelder;
fval_nelder;

% Determine objective over mesh 
Fobj=sf(R_mesh,C_mesh);

% Genetic Algorithm
rng default
[x_ga,fval_ga,exitflag_ga,output_ga]=ga(func,2,[],[],[],[],[R_min C_min],[R_max C_max]);

% We can not really build a convergence path for this algo
options_ps = optimoptions('particleswarm','Display','Iter','OutputFcn',@pswplotranges);
[x_ps,fval_ps,exitflag_ps,output_ps]=particleswarm(func,2,[R_min C_min],[R_max C_max],options_ps);


figure()
hold on
contour(R_mesh,C_mesh,sf(R_mesh,C_mesh),30,'ShowText','on')
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
hold on 
contour(R_mesh,C_mesh,sg1(R_mesh,C_mesh),[0.0 0.0],'b')
contour(R_mesh,C_mesh,sg2(R_mesh,C_mesh),[0.0 0.0],'r')
contour(R_mesh,C_mesh,sg1(R_mesh,C_mesh),[0.1 0.1],'--b')   
contour(R_mesh,C_mesh,sg2(R_mesh,C_mesh),[0.1 0.1],'--r')  
scatter(0.4304,0.0313,'cd','filled')
legendInfo={'Objective function','Constrain 1','Constrain 2','Active Const 1','Active Const 1','Global Minimum','Convergence Path Nelder','Convergence Part Swarm'};
for i=[0:output_nelder.iterations] 
    data=load(['Own_Optimisation\Iterative_Results\results_iteration_nelder',num2str(i)]);
    r_iter(i+1)=data(1);
    c_iter(i+1)=data(2);
end
plot(r_iter,c_iter,'m')
for i=[0:output_ps.iterations]
    data=load(['Own_Optimisation\Iterative_Results\results_iteration_ps',num2str(i)]);
    r_iter(i+1)=data(1);
    c_iter(i+1)=data(2);
end 
plot(r_iter,c_iter,'g')
legend(legendInfo)
caxis([min(Fobj(:)) max(Fobj(:))])

%% Output function to write iterative data
function stop = myoutput_own(x,optimvalues,state);
    stop = false;
    if isequal(state,'iter')
      results=[x,optimvalues.fval];
      dlmwrite(['Own_Optimisation/Iterative_Results/results_iteration_sqp',mat2str(optimvalues.iteration)],results)
    end
end

function stop = myoutput_nelder(x,optimvalues,state);
    stop = false;
    if isequal(state,'iter')
      results=[x,optimvalues.fval];
      dlmwrite(['Own_Optimisation/Iterative_Results/results_iteration_nelder',mat2str(optimvalues.iteration)],results)
    end
end

function stop = pswplotranges(optimValues,state)
optimValues;
coor=reshape(optimValues.swarm,[1,numel(optimValues.swarm)]);
value=reshape(optimValues.swarmfvals,[1,numel(optimValues.swarmfvals)]);
results=[optimValues.bestx,optimValues.bestfval,coor,value];
dlmwrite(['Own_Optimisation/Iterative_Results/results_iteration_ps',mat2str(optimValues.iteration)],results)
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
