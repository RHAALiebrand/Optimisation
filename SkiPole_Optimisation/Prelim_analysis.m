%% Preliminary analysis
% RUN MAIN TO RUN THIS

N=3;

% h/c   vs    r/c 
H=linspace(0.1,0.45,N);
R=linspace(0.3,0.5,N);
Fobj=zeros(length(R),length(H));
G1=zeros(length(R),length(H));
G2=zeros(length(R),length(H));
G3=zeros(length(R),length(H));
[H_mesh,R_mesh]=meshgrid(H,R);
for i=[1:length(H)]
    i
    for j=[1:length(R)]
        Fobj(j,i)=fobj(R(j),H(i),c,t,F);
        g=g_i(R(j),H(i),c,t);
        G1(j,i)=g(1);
        G2(j,i)=g(2);
        G3(j,i)=g(3);
    end
end

figure()
contour(H_mesh,R_mesh,Fobj,'ShowText','on')
xlabel('$h/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
hold on 
contour(H_mesh,R_mesh,G1,[0.0 0.0],'b')
contour(H_mesh,R_mesh,G2,[0.0 0.0],'r')
contour(H_mesh,R_mesh,G3,[0.0 0.0],'c')
contour(H_mesh,R_mesh,G1,[0.1 0.1],'--b')
contour(H_mesh,R_mesh,G2,[0.1 0.1],'--r')
legend('Objective function','Constrain 1')

% h/c   vs c
H=linspace(0.1,0.45,N);
C=linspace(5e-3,4e-2,N);
Fobj=zeros(length(C),length(H));
G1=zeros(length(C),length(H));
G2=zeros(length(C),length(H));
G3=zeros(length(C),length(H));
[H_mesh,C_mesh]=meshgrid(H,C);
for i=[1:length(H)]
    i
    for j=[1:length(C)]
        Fobj(j,i)=fobj(r,H(i),C(j),t,F);
        g=g_i(r,H(i),C(j),t);
        G1(j,i)=g(1);
        G2(j,i)=g(2);
        G3(j,i)=g(3);
    end
end
figure()
contour(H_mesh,C_mesh,Fobj,'ShowText','on')
xlabel('$h/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
hold on 
contour(H_mesh,C_mesh,G1,[0.0 0.0],'b')
contour(H_mesh,C_mesh,G2,[0.0 0.0],'r')
contour(H_mesh,C_mesh,G3,[0.0 0.0],'c')
contour(H_mesh,C_mesh,G1,[0.1 0.1],'--b')
contour(H_mesh,C_mesh,G2,[0.1 0.1],'--r')
legend('Objective function','Constrain 1','Constrain 2')


% r/c   vs  c
R=linspace(0.3,0.5,N);
C=linspace(5e-3,4e-2,N);
Fobj=zeros(length(C),length(R));
G1=zeros(length(C),length(R));
G2=zeros(length(C),length(R));
G3=zeros(length(C),length(R));
[R_mesh,C_mesh]=meshgrid(R,C);
for i=[1:length(R)]
    i
    for j=[1:length(C)]
        Fobj(j,i)=fobj(R(i),h,C(j),t,F);
        g=g_i(R(i),h,C(j),t);
        G1(j,i)=g(1);
        G2(j,i)=g(2);
        G3(j,i)=g(3);
    end
end

figure()
contour(R_mesh,C_mesh,Fobj,20,'ShowText','on')
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
hold on 
contour(R_mesh,C_mesh,G1,[0.0 0.0],'b')
contour(R_mesh,C_mesh,G2,[0.0 0.0],'r')
contour(R_mesh,C_mesh,G3,[0.0 0.0],'c')
contour(R_mesh,C_mesh,G1,[0.1 0.1],'--b')
contour(R_mesh,C_mesh,G2,[0.1 0.1],'--r')
legend('Objective function','Constrain 1','Constrain 2')
caxis([min(Fobj(:)) max(Fobj(:))])

