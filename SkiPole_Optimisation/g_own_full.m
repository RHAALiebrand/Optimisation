%% Engineering optimization 
% Constraints for full optimisation since the fit polynominal needs another
% syntax 
% Martin Janssens
% Rens Liebrand

function [g,h]=g_own(x,sg1,sg2)
    g(1)=sg1(x(1),x(2),x(3),x(4));
    g(2)=sg2(x(1),x(2),x(3),x(4));
    h=[]; 
end