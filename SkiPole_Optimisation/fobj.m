%% Engineering optimization 
% Objective function
% Martin Janssens
% Rens Liebrand
function [f]=fobj(r,h,c,t,F)
Params;
D=Aerodynamic_model(r,h,c,F);
[m,g1,g2,g3]=Structural_model(r,h,c,t);
f=D+m*nu_snow*G;
end 

