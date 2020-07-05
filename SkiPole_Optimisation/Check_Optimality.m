function [] = Check_Optimality(x,g1,g2,dfdx,dgdx,sf,sg1,sg2,h,active)
% Check feasibility
disp(['g1 = ', num2str(g1)]);
disp(['g2 = ', num2str(g2)]);

if g1 > 1e-1 || g2 > 1e-1
    disp('The optimum does not satisfy the constraints to numerical tolerance')
else
    disp('The optimum satisfies the constraints')
end

% Compute and evaluate Lagrangian multipliers
if all(active == 0)
    % Complementary slackness - inactive multipliers are 0
    disp('Interior optimum - No Lagrangian multipliers needed for optimality')
else
    if g1 > 0
        gbar = g1;
        dgbardx = dgdx(1,:);
        Hg = hessian(x,h,sg1);
    elseif g2 > 0
        gbar = g2;
        dgbardx = dgdx(2,:);
        Hg = hessian(x,h,sg2);
    end
    mu = dgbardx'\-dfdx; % Least squares mu, since m < n
    if any(mu < 0)
        disp('Negative Lagrangian multiplier - could decrease f without violating constraints')
    else
        disp('The optimum is feasible')
    end
end

% Check whether it is an optimum
Hf = hessian(x,h,sf);
HL = Hf + mu'*Hg;
[~,p] = chol(HL);
if p == 0
    disp('Hessian is positive definite, this is the optimum!')
else
    disp('Hessian is not positive definite, this is not the optimum!')
end
end