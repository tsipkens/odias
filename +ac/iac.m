
% IAC  Iterative-average charge algorithm.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

function [xbar, qbar0] = iac(x_star, qmodel)

qmax = 100;
qvec = 1:qmax;

% Initialize quantities.
xi = x_star;  % start with size equivalent to size-to-charge setpoint
xi1 = 0;
f_conv = zeros(size(xi));  % flag for converge entries, none at the start
qbar0 = f_conv;  % master list of average charge (qbar is subset being updated below)

% Vectorized form for vector m_star using Newton's method.
disp(' ');
tools.textheader('Running IAC');
disp(' ');
disp(' Starting iteration...');
disp(' ');
n = 50;
for ii=1:n
    
     % Physical model.
    qbar = qmodel(xi(~f_conv), qvec);
    
    % If qbar exceed the integration duration.
    % Error check, avoids erroneous truncation of
    % integer charge state during charge model (Fuch's) call.
    f_above = (qbar + sqrt(qbar)) > qmax;  % flag for cases that did not converge due to low zmax
    if any(f_above)
        qbar(f_above) = min(qbar(f_above), qmax);
        qmax = round(qmax * 1.5 / 10) * 10;  % increase zmax by 50%
        qvec = 1:qmax;
        
        disp([' Adjusted zmax to ', num2str(qmax)]);
    end
    qbar0(~f_conv) = qbar;  % copy over new values
    
    
    % Update size for unconverged.
    xi(~f_conv) = x_star(~f_conv) .* qbar;
    
    % Convergence check/components.
    if all(abs((xi1 - xi) ./ xi1) < 0.001)
        break;
    end
    
    f_conv = abs((xi1 - xi) ./ xi1) < 0.001;  % flag converged entries and skip next time
    xi1 = xi;  % temporarily store previous value (for convergence check)
    
end

disp(' Iteration complete.');
disp(' ');

% Check if algorithm did not converge, 
% which occurred if ii = n.
if ii==n
    warning('Algorithm did not converge!');
end
%}

xbar = xi;

tools.textheader();

end
