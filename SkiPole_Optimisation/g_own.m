%% Engineering optimization 
% Constraints for simplified 2D problem
% Martin Janssens
% Rens Liebrand

function [g,h]=g_own(x,sg1,sg2)
    g(1)=sg1(x);
    g(2)=sg2(x);
    h=[]; 
end