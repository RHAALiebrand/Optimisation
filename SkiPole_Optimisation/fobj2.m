%% Engineering optimization 
% Objective function in 2D so for simplified problem
% Martin Janssens
% Rens Liebrand


function [f]=fobj2(x) 
r=x(1);
c=x(2);
t=3e-3;
F=import_database();
h=0.25;
Params;
D=Aerodynamic_model(r,h,c,F);
[m,g1,g2,g3]=Structural_model(r,h,c,t);
f=D+m*nu_snow*G;
end 
