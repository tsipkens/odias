
% TIKHONOV_LPR  Generates Tikhonov matrices/operators. 
%  
%  L = tikhonov_lpr(ORDER, X_LENGTH) generates the Tikhonov matrix, L, of
%  the order specified in ORDER and corresponding to a vector x of length
%  X_LENGTH, such that L*x = 0.
%  
%  L = tikhonov_lpr(ORDER, X_LENGTH, BC) add an input for the type of
%  boundary conidition to be applied at the edges.
%  
%  ------------------------------------------------------------------------
%
%  ORDER options:
%   0: 0th order Tikhonov promotes small solutions. 
%   1: 1st order Tikhonov minimizes the first derivative.
%   2: 2nd order Tikhonov minimized the second derivative.
%  
%  BC (boundary condition) options:
%   []: No explicit boundary conditions (uses that implied by order).
%    0: Zero BCs.
%    1: No slope BCs
%    2: No curvature BCs.
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR:   Timothy Sipkens, 2020-04-11

function L = tikhonov_lpr(order, x_length, bc)

% Parse type input.
if ~exist('bc', 'var'); bc = []; end
if isempty(bc); bc = order; end  % neglects boundary conditions

% Choose between order of Tikhonov operator to generate.
switch order
    case 0 % 0th order Tikhonov
    	L = speye(x_length);
        
    case 1 % 1st order Tikhonov
        L = -speye(x_length);
        L = spdiags(ones(x_length, 1), 1, L);
        
        if bc~=0
            L(end,:) = [];
        end
        
    case 2 % 2nd order Tikhonov
        L = -speye(x_length);
        L = spdiags(0.5 .* ones(x_length,2), [-1,1], L);
        
        if bc == 0  % zero
            L(1, 2) = 0;
            L(x_length, x_length - 1) = 0;
       elseif bc == 1  % no slope
            L(1, 2) = 1;
            L(x_length, x_length - 1) = 1;
        else
            L(1,:) = [];
            L(end,:) = [];
        end
        
    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end

end


