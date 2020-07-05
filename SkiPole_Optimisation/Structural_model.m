%% Engineering optimization 
% Structural model
% Martin Janssens
% Rens Liebrand

function [m,g1,g2,g3]=Structural_model(r,h,c,t)
Params;

[A,I_x,I_y,x_bar,y_bar]=A_I_calc(r,h,c,t);

% Determine mass, which is an objective
m=A*L*rho_carbon;

% Determine moment arms
ex=abs(0.9*c-x_bar);
ey=abs(0.25*h*c-y_bar);

% Evaluate constains
g1=1-pi^2*E*min(I_x,I_y)/(780*K*L^2);      % Critical load
g2=ex/(18e-3)*(sec(sqrt(P_skier/(E*I_y))*L/2-1))-1;
g3=ey/(18e-3)*(sec(sqrt(P_skier/(E*I_x))*L/2-1))-1;
end 
