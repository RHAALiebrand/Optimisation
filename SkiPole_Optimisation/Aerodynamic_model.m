%% Engineering optimization 
% Aerodynamics model
% Martin Janssens
% Rens Liebrand

% In/outputs
% Cd ~ h,c,r
% D=0.5rho*Uinf*h

function [D]=Aerodynamic_model(r,h,c,F)
Params;
Cd=F(r,h,c);
D=0.5*rho_air*Uinf^2*h*L*Cd;
end 

