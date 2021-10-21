
% MARKOWSKI Funtion to perform Markowski/Buckley-type smoothing.
%   Used in the Twomey-Markowski algorithm.
% Author: Timothy Sipkens, 2020-02-06
% Note: n can be an integer or a Grid
%=========================================================================%

function [x, G, SIGMA] = markowski(A, b, x, n, iter, opt_smooth, Sf, SIGMA_end)

xi = x;

G = G_Markowski(n);

%-- Perform smoothing over multiple iterations ---------------------------%
jj = 1;
while jj<=iter
    x = G * x; % apply smoothing
    
    % Exit smoothing if mean square error exceeds unity.
    if tools.mean_sq_err(A, x, b) > 2
        return;
    end
    
    jj = jj+1;
end

disp(['SMOOTHING algorithm did not necessarily converge after ',num2str(iter),' iteration(s).']);

end



%== G_MARKOWSKI ==========================================================%
%   Generates a smoothing matrix of the form based on original work of Markowski
function G = G_Markowski(n)

G = invert.tikhonov_lpr(2, n, 1);

G = abs(G) ./ 2;
G(1, 2) = 1/4;
G(1, 1) = 3/4;
G(end, end-1) = 1/4;
G(end, end) = 3/4;

end
%=========================================================================%


