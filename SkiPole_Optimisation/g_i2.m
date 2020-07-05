%% Engineering optimization 
% Constraints
% Martin Janssens
% Rens Liebrand

function [g,h]=g_i2(x)
r=x(1);
c=x(2);
t=3e-3;
h=0.25;
[m,g(1),g(2),g(3)]=Structural_model(r,h,c,t);
%g=[g1;g2;g3];
h=[];
end 