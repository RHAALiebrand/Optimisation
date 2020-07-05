%% Engineering optimization 
% Constraints
% Martin Janssens
% Rens Liebrand

function [g]=g_i(r,h,c,t)
[m,g1,g2,g3]=Structural_model(r,h,c,t);
g=[g1;g2;g3];
end 