function [H] = hessian(X,h,sf)

    H = zeros(length(X));

    % for each dimension of objective function
    for i=1:length(X)
        % derivative at first point (left)
        x1 = X;
        x1(i) = X(i) - h;
        [df1, dum] = gradients(x1,1e-8,sf,sf,sf);

        % derivative at second point (right)
        x2 = X;
        x2(i) = X(i) + h;
        [df2, dum] = gradients(x2,1e-8,sf,sf,sf);

        % differentiate between the two derivatives
        d2f = (df2-df1) / (2*h);

        % assign as row i of Hessian
        H(i,:) = d2f';
    end
end