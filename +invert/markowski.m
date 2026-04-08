
% MARKOWSKI Funtion to perform Markowski smoothing.
%  Used in the Twomey-Markowski algorithm.
%  
%  AUTHOR: Timothy Sipkens, 2020-02-06
%  NOTE: n can be an integer or a Grid

function [x, G, err] = markowski(A, b, x, n, iter, opt_smooth, ~, ~)

if ~exist('opt_smooth', 'var'); opt_smooth = []; end
if isempty(opt_smooth); opt_smooth = 3; end

err = [];

if opt_smooth == 3
    G = G3(n);
else
    G = G_Markowski(n);
end

xi = x;

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

err = ['SMOOTHING algorithm did not necessarily converge after ', ...
    num2str(iter), ' iteration(s).'];

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



%== G_MARKOWSKI ==========================================================%
%   Generates a smoothing matrix of the form based on original work of Markowski
function G = G3(n)

G = diag(ones(n-2,1), -2) + 2 .* diag(ones(n-1,1), -1) + ...
    3 .* diag(ones(n,1)) + 2 .* diag(ones(n-1,1), 1) + diag(ones(n-2,1), 2);
G = G ./ sum(G, 2);

end
%=========================================================================%


