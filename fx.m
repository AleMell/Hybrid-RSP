function cost_function_value_at_x_plus_noise = fx(x,Delta,Z, Delta_j)

% The cost/objective function

global bar_ns Mode lambda_s PhiMin_theoretical noise

% QUADRATIC

cost_function_value_at_x=1.*x(1).^2+5.*x(2).^2;

% ADJIMAN

%cost_function_value_at_x=(cos(x(1)) .* sin(x(2))) - (x(1) ./ ((x(2) .^ 2) + 1));

% DROP WAVE FUNCTION

%cost_function_value_at_x = - (1+cos(12*(x(1)^2+x(2)^2)^(1/2)))/(0.5*(x(1)^2+x(2)^2)+2);


% ROSENBROCK

% plot Rosenbrock 
%[x,y] = meshgrid(-3.5:.2:3.5);
%surf(x,y,(1-x).^2+100*(y-x.^2).^2)

%cost_function_value_at_x = (1-x(1))^2+10*(x(2)-x(1)^2)^2;

% ACKLEY

% plot Ackley 2D
% [x,y] = meshgrid(-10:.2:10);
% surf(x,y, -20*exp(-0.2*sqrt(1/n*(x.^2+y.^2)))-exp(1/n*(cos(2*pi*x)+cos(2*pi*y))))

%cost_function_value_at_x = -20*exp(-0.2*sqrt(1/n*sumsqr(x)))-exp(1/n*sum(cos(2*pi*x)));

% NOISE GENERATION

epsilon = 0.0001;

if (Mode==3)
    %random noise 
    
    noise = -bar_ns + (2*bar_ns).*rand;
    
elseif (Mode == 2)
    
    % exploding noise
    
    if(max(Delta_j)>=lambda_s*PhiMin_theoretical)
        noise = 0;
    else
        if (Z-noise>=cost_function_value_at_x - rho(Delta))
            
            if (cost_function_value_at_x - Z -rho(Delta) + epsilon <= bar_ns && cost_function_value_at_x - Z -rho(Delta) + epsilon >= -bar_ns)
                
                noise = cost_function_value_at_x - Z -rho(Delta)+epsilon;
            
            elseif (cost_function_value_at_x - Z -rho(Delta) + epsilon >= bar_ns)
                
                noise = bar_ns;
                
            else 
                noise = -bar_ns;
            end
            
        elseif (Z-noise<=cost_function_value_at_x - rho(Delta))   
            
            if (-cost_function_value_at_x + Z -rho(Delta) - epsilon <= bar_ns && -cost_function_value_at_x + Z -rho(Delta) - epsilon >= -bar_ns)
                
                noise = -cost_function_value_at_x + Z -rho(Delta)-epsilon;
            
            elseif (-cost_function_value_at_x + Z -rho(Delta) - epsilon >= bar_ns)
                
                noise = bar_ns;
                
            else 
                noise = -bar_ns;
            end
        
        end
    end
    
    %random noise 
    
    %noise = -bar_ns + (2*bar_ns).*rand;
    
else
    noise = 0;
end


cost_function_value_at_x_plus_noise = cost_function_value_at_x + noise;

end

