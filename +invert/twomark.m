
% TWOMARK  Performs inversion using the iterative Twomey-Markowski approach.
%  This add an intermediate smoothing step to the standard Twomey routine. 
%  
%  ------------------------------------------------------------------------
%  
%  INPUTS:
%   A           Model matrix
%   b           Data
%   n           Length of first dimension of solution, used in smoothing
%   xi          Initial guess
%
%  OUTPUTS:
%   x           Estimate
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR:  Timothy Sipkens, 2018-12-20

function x = twomark(A, b, n, xi)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('Sf','var')
    Sf = 1/300;  % matching Buckley et al.
elseif isempty(Sf)
    Sf = 1/300;
end

if ~exist('opt_smooth','var')
    opt_smooth = 'Buckley';
elseif isempty(opt_smooth)
    opt_smooth = 'Buckley';
end
%-------------------------------------------------------------------------%

iter = 40;  % maximum overall Twomey-Markowski iterations
iter_two = 150;  % max number of iterations in Twomey pass

x = xi;
x = invert.twomey(A, b, x, iter_two);  % initial Towmey procedure
R = roughness(x);  % roughness vector

tools.textbar([0, iter]);
iter_two = 150;  % max number of iterations in Twomey pass
for kk=1:iter  % iterate Twomey and smoothing procedure
    x_temp = x;  % store temporarily for the case that roughness increases
    x = invert.markowski(A, b, x, n, 1000, opt_smooth, Sf, 1);  % perform smoothing
    
    x = invert.twomey(A, b, x, iter_two);  % perform Twomey
    
    %-- Check roughness of solution ------------------%
    R(kk+1) = roughness(x);
    if R(kk+1) > R(kk)  % exit if roughness has stopped decreasing
        tools.textbar([iter, iter]);
        disp([' Exited Twomey-Markowski loop after ',num2str(kk),...
            ' iterations because roughness increased.']);
        x = x_temp;  % restore previous iteration
        break;
    end
    
    tools.textbar([kk, iter]);
end

disp(' Completed Twomey-Markowski procedure:');
disp(['  iter = ', num2str(kk)]);
disp(' ');

end
%=========================================================================%



%== ROUGHNESS ============================================================%
%   Computes an estimate of the roughness of the solution.
%   This function is used for convergence in the Twomey-Markowski loop and
%   is based on the average, absolute value of the second derivative.
function R = roughness(x)

R = abs(x(3:end) + ...
         x(1:(end-2),:) - ...
    2 .* x(2:(end-1)));
R = sum(R) ./ (length(x) - 2);

end


