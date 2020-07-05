%% Problem investigation
% RUN MAIN TO RUN THIS
% We go and check the problem details
bound_margin=0.9;

Nc=3;
C_char_res=linspace(8e-3,bound_margin*4e-2,Nc)';
Fobj_C=zeros(Nc,1);

G1_C=zeros(Nc,1);
G2_C=zeros(Nc,1);
for i=[1:Nc]
    i
    [Fobj_C(i)]=fobj(r,h,C_char_res(i),t,F);
    g=g_i(r,h,C_char_res(i),t);
    G1_C(i)=g(1);
    G2_C(i)=g(2);
end

figure()
plot(C_char_res,Fobj_C)
grid
xlabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
ylabel('$f$ [N]','fontsize',16,'Interpreter','LaTex')
xlim([min(C_char_res) max(C_char_res)])

figure()
plot(C_char_res,G1_C)
hold on
plot([min(C_char_res) max(C_char_res)],[0 0],'--k')
grid
xlabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
ylabel('$g_1$ [-]','fontsize',16,'Interpreter','LaTex')
xlim([min(C_char_res) max(C_char_res)])

figure()
plot(C_char_res,G2_C)
hold on
plot([min(C_char_res) max(C_char_res)],[0 0],'--k')
grid
xlabel('$c$ [m]','fontsize',16,'Interpreter','LaTex')
ylabel('$g_2$ [-]','fontsize',16,'Interpreter','LaTex')
xlim([min(C_char_res) max(C_char_res)])

RC_char_res=linspace((2-bound_margin)*0.3,bound_margin*0.5,Nc)';
Fobj_RC=zeros(Nc,1);
G1_RC=zeros(Nc,1);
G2_RC=zeros(Nc,1);
for i=[1:Nc]
    i
    Fobj_RC(i)=fobj(RC_char_res(i),h,c,t,F);
    g=g_i(RC_char_res(i),h,c,t);
    G1_RC(i)=g(1);
    G2_RC(i)=g(2);
end

figure()
plot(RC_char_res,Fobj_RC)
grid
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$f$ [N]','fontsize',16,'Interpreter','LaTex')
xlim([min(RC_char_res) max(RC_char_res)])
figure()
plot(RC_char_res,G1_RC)
hold on
grid
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$g_1$ [-]','fontsize',16,'Interpreter','LaTex')
xlim([min(RC_char_res) max(RC_char_res)])

figure()
plot(RC_char_res,G2_RC)
hold on
grid
xlabel('$r/c$ [-]','fontsize',16,'Interpreter','LaTex')
ylabel('$g_2$ [-]','fontsize',16,'Interpreter','LaTex')
xlim([min(RC_char_res) max(RC_char_res)])
