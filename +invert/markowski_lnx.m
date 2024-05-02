
% MARKOWSKI_LNX Funtion to perform Markowski smoothing.
%  Used in the Twomey-Markowski algorithm.
%  
%  AUTHOR: Timothy Sipkens, 2020-02-06
%  NOTE: n can be an integer or a Grid

function [x, G] = markowski_lnx(A, b, x, n, iter, ~, ~, ~)

xi = x;

%-- Perform smoothing over multiple iterations ---------------------------%
jj = 1;
while jj<=iter
    
    fl0 = x ~= 0;
    G = G_Markowski(sum(fl0));
    
    x(fl0) = exp(G * log(x(fl0))); % apply smoothing
    
    % Exit smoothing if mean square error exceeds unity.
    if tools.mean_sq_err(A, x, b) > 2
        return;
    end
    
    jj = jj+1;
end

disp(['SMOOTHING algorithm did not necessarily converge after ', ...
    num2str(iter), ' iteration(s).']);

end



%== G_MARKOWSKI ==========================================================%
%   Generates a smoothing matrix of the form based on original work of Markowski
function G = G_Markowski(n)

G = invert.tikhonov_lpr(2, n, 1);

G = abs(G) ./ 4;
G(1, 2) = 1/4;
G(1, 1) = 3/4;
G(end, end-1) = 1/4;
G(end, end) = 3/4;

end
%=========================================================================%


