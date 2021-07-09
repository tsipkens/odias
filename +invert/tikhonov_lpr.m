
% TIKHONOV_LPR  Generates Tikhonov matrices/operators. 
%  
%  L = tikhonov_lpr(ORDER, X_LENGTH) generates the Tikhonov matrix, L, of
%  the order specified in ORDER and corresponding to a vector x of length
%  X_LENGTH, such that L*x = 0.
%  
%  ------------------------------------------------------------------------
%
%  ORDER options:
%   0: 0th order Tikhonov promotes small solutions. 
%   1: 1st order Tikhonov minimizes the first derivative.
%   2: 2nd order Tikhonov minimized the second derivative.
%  
%  ------------------------------------------------------------------------
% 
%  AUTHOR:   Timothy Sipkens, 2020-04-11

function L = tikhonov_lpr(order, x_length)

% Choose between order of Tikhonov operator to generate.
switch order
    case 0 % 0th order Tikhonov
    	L = speye(x_length);
        
    case 1 % 1st order Tikhonov
        L = -speye(x_length);
        L = spdiags(ones(x_length, 1), 1, L);
        L(end,:) = [];
        
    case 2 % 2nd order Tikhonov
        L = -speye(x_length);
        L = spdiags(0.5 .* ones(x_length,2), [-1,1], L);
        L(1,:) = [];
        L(end,:) = [];
        
    otherwise
        disp('The specified order of Tikhonov is not available.');
        disp(' ');
        return
end

end


