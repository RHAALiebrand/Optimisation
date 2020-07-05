function [x,active,f,g1,g2,dfdx,dgdx]=Steepest_Decent_Bound(sf,x_start,N_its,N_fibo,N_bnd,line_factor,sg1,sg2,xmin,xmax)
    % Set outputs
    x = x_start;
    active = [0 0]; % Toggles for whether [g1 g2] are active - both 0 
                    % (inactive) if we assume a feasible x_start    
    % Looping process
    for i=[1:N_its]
        if i~=1
            x0=x_new;
        else 
            x0=x_start;
        end 
        
        % Directions
        [dfdx,dgdx] = gradients(x0,1e-8,sf,sg1,sg2);
        % If no constraint is active, steepest descent along gradf
        if all(active == 0)
            d = -[dfdx(1); dfdx(2)];
            d_unit = d/norm(d);
            x1 = x0 - d_unit.*x0*line_factor;
            x2 = x0 + d_unit.*x0*line_factor;
            
            % Ensure the design bounds are not violated along the line
            if any(x1>xmax)
                sub = find(x1>xmax);
                x1(sub) = xmax(sub);
            end
            if any(x2>xmax)
                sub = find(x2>xmax);
                x2(sub) = xmax(sub);
            end
            if any(x1<xmin)
                sub = find(x1<xmin);
                x1(sub) = xmin(sub);
            end
            if any(x2<xmin)
                sub = find(x2<xmin);
                x2(sub) = xmin(sub);
            end
            
            % Fibo search
            [r_new,c_new] = Fibonacci(x1,x2,N_fibo,sf);        
            x_new = [r_new c_new]';
            
            % Evaluate the constraints at the new point
            g1 = sg1(x_new');
            g2 = sg2(x_new');
            
            sd_it = 1; % Is this a steepest descent iteration?
        end
        
        % If a constraint is active, switch to gradient projection
        if g1 > 0 || g2 > 0
            % Which constraint is active? 
            if g1 > 0
                hx = g1;
                dhdx = dgdx(1,:);
                active(1) = 1;
            elseif g2 > 0   
                hx = g2;
                dhdx = dgdx(2,:);
                active(2) = 1;
            elseif g1 > 1 && g2 > 0
                "Both constraints are violated - don't know what to do!"
            end              
            % If you're on the boundary, project along it!
            % If you violated a constraint due to steepest descent, just
            % return to the boundary without projection
            if sd_it == 0
                P = 1 - dhdx'/(dhdx*dhdx')*dhdx;
                d = P*dfdx;
                d_unit = d/norm(d);
                x1 = x0 - d_unit.*x0*line_factor;
                x2 = x0 + d_unit.*x0*line_factor;
                
                % Again ensure the design bounds are not violated
                if any(x1>xmax)
                    sub = find(x1>xmax);
                    x1(sub) = xmax(sub);
                end
                if any(x2>xmax)
                    sub = find(x2>xmax);
                    x2(sub) = xmax(sub);
                end
                if any(x1<xmin)
                    sub = find(x1<xmin);
                    x1(sub) = xmin(sub);
                end
                if any(x2<xmin)
                    sub = find(x2<xmin);
                    x2(sub) = xmin(sub);
                end
                
                [r_new,c_new] = Fibonacci(x1,x2,N_fibo,sf);
                x_new = [r_new,c_new]';    
            end      
            x = [x x_new];
            
            % Since a constraint is violated, Newton iterations to get back 
            % to the constraint boundary      
            for i=[1:N_bnd]               
                % Update x with a Newton step towards boundary
                ck = -inv(dhdx*dhdx')*hx;
                x_new = x_new + dhdx'*ck;  
                
                % Compute h and dhdx at x_new
                [dfdx,dgdx] = gradients(x_new,1e-8,sf,sg1,sg2);
                if g1 > 0
                    hx = sg1(x_new');
                    dhdx = dgdx(1,:);
                elseif g2 > 0
                    hx = sg2(x_new');
                    dhdx = dgdx(2,:);
                end
%                 x = [x x_new];
            end
        else
            active = [0 0];
        end
        x = [x x_new];
        sd_it = 0;
    end 
    f = sf(x_new');
end 