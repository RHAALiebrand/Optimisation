function [dfdx,dgdx] = gradients(x,h,f,g1,g2);
    % Returns the vector grad(f) and the matrix dgdx (of size ndv*ncons) at 
    % a given point x, by forward differences (I know I know...) Currently
    % hardcoded for 2 constraint functions

    % Input:
    % x  : design point (column vector) for which derivatives are computed.
    % h  : spacing between points
    % f  : function that is the objective function
    % g1 : function that returns the value of g1
    % g2 : ditto

    % Output:
    % gradf  : gradient (row)vector "[dfdx1 dfdx2]" of objective function.
    % dgdx   : matrix of dgdx 
    
    ndv = length(x);
    fx = f(x');
    g1x = g1(x');
    g2x = g2(x');
    
    dfdx = zeros(ndv,1);
    dgdx = zeros(ndv,2);

    % Compute gradients   
    for i=1:ndv
        if i == 1
            xi = [x(1)+h x(2:end)'];
        elseif i == length(x)
            xi = [x(1:end-1)' x(end)+h];
        else
            xi = [x(1:i-1)' x(i)+h x(i+1:end)'];
        end

        fxplush = f(xi);
        dfdx(i) = (fxplush - fx)/h;
        
        g1xplush = g1(xi);
        dgdx(1,i) = (g1xplush - g1x)/h;
        
        g2xplush = g2(xi);
        dgdx(2,i) = (g2xplush - g2x)/h;
    end    
end