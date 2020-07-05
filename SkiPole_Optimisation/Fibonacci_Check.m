%% Check for the fibonacci algo
%  We check by plotting some lines with their optimium

% Set avg parameters
t=3e-3;
r=0.4; % r   *c
h=0.25; % h   *c
c=0.012;

% Reading the above iso looping
Fobj=dlmread('Own_Optimisation/Fobj_200');
G1=dlmread('Own_Optimisation/G1_200');
G2=dlmread('Own_Optimisation/G2_200');
G3=dlmread('Own_Optimisation/G3_200');
C=dlmread(['Own_Optimisation/C_',mat2str(N)]); % We also have to write+read these otherwise these change
R=dlmread(['Own_Optimisation/R_',mat2str(N)]); 
sf = fit([R, C],Fobj','poly22');

% CG Algo
x0=[0.35 ; 0.02];
[df_dx,df_dy]=differentiate(sf,x0(1),x0(2));
d_1=-[df_dx; df_dy];
d_1_unit=d_1/norm(d_1);
x1=x0+d_1_unit.*x0*0.8; 

% We want to perform a line search do we have to set up C=aR+b where a is
% the slope
RC=(x1(2)-x0(2))/(x1(1)-x0(1));
constant=x0(2)-RC*x0(1);

% Plot the line search
M=10;
R_disc=linspace(x0(1),x1(1),M);
C_disc=RC.*R_disc+constant;
plot(sf,[R,C],Fobj')
hold on
for i=[1:M]
    plot1=scatter3(R_disc(i),C_disc(i),sf([R_disc(i) C_disc(i)]),'p','filled');
    M1='First line search';
end

% Perform line search
N_iter=20;
[r_new,c_new]=Fibonacci(x0,x1,N_iter,sf);

%% Second iteration
l=sqrt((x1(1)-x0(1))^2+(x1(2)-x0(2))^2); % Scaling term for bracketing
x1=[r_new ; c_new];

% Second direction
[df_dx,df_dy]=differentiate(sf,x1(1),x1(2));
d_2=-[df_dx; df_dy];
d_2_unit=d_2/norm(d_2);
x2=x1+d_2_unit.*l*0.2;

% New C,R relation
RC=(x2(2)-x1(2))/(x2(1)-x1(1));
constant=x1(2)-RC*x1(1);

% Plotting
R_disc=linspace(x1(1),x2(1),M);
C_disc=RC.*R_disc+constant;

% Plot line search results
plot2=scatter3(r_new,c_new,sf([r_new c_new]),'rd','filled');
M2='First Fibonacci result';
for i=[2:M]
    plot3=scatter3(R_disc(i),C_disc(i),sf([R_disc(i) C_disc(i)]),'s','filled');
    M3='Second line search';
end
[r_new,c_new]=Fibonacci(x1,x2,N_iter,sf);

plot4=scatter3(r_new,c_new,sf([r_new c_new]),'cd','filled');
M4='Second Fibonacci result';
legend([plot1, plot2, plot3, plot4],M1, M2, M3, M4)

xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')

