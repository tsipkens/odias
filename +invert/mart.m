
% MART  Implementation of the MART algorithm.
% 
%--------------------------------------------------------------------------
% 
%  Inputs:
%   A       Model matrix
%   b       Data
%   xi      Initial guess       (Optional, default ones vector)
%   iter    MART iterations     (Optional, default 10)
%   sigma_fun     Function of x to be evaluated to determine convergence
%                   (Optional, default: ignore)
%   sigma         Value to which sigma_fun is compared to determine convergence
%                   (Optional, default: ignore)
%   f_bar         Boolean to determine whether or not to show textbar
%                   (Optional, default: 0)
%
%  Output:
%   x       MART estimate
% 
% -------------------------------------------------------------------------
% 
%  SEE: Gordon et al, J. Theor. Biol. 29(3), 471-481 (1970)
%       See also Katsekpor, Oper. Res. Forum (2023).
%  AUTHOR: Timothy Sipkens, 2024-03-21
%  Adapted from Samuel Grauer

function [x] = mart(A, b, xi, iter, f_bar)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('iter', 'var'); iter = []; end
if isempty(iter); iter = 1000; end % default of 10 iterations

if ~exist('xi', 'var'); xi = []; end
if ~isempty(xi)
    x = xi;
else
    x = ones(size(A,2), 1); % intiate vector of ones if xi is not specified
end

% whether to display text-based progress bar
if ~exist('f_bar', 'var'); f_bar = []; end
if isempty(f_bar); f_bar = 0; end
%-------------------------------------------------------------------------%



%-- Get ray weights -------------%
s = 1 ./ max(A, [], 2); % factor to scale data and model matrix
s_min = max(max(A,[],2)) * 1e-12;
s(max(A,[],2) < s_min) = 0;
A = bsxfun(@times, s, A); % scale model matrix
b = max(b.*s, 0); % remove any negative data and scale data
lam = 0.002; % relaxation parameter, order unity


%-{
%-- Note: Uncomment to switch to the normal implementation.
%-- MART iterations: Normal implementation -------------%
for kk=1:iter
    xm = x;  % x from previous step

    for ii = 1:length(b)
        if b(ii)~=0 % do not consider null data
            y = A(ii,:) * x; % scalar value
            if y~=0 % avoid division by zero, do not update those entries
                x = x .* exp(lam .* A(ii,:)' * log(b(ii) / y));
            end
        end
    end
    
    if f_bar; tools.textbar(kk/iter); end % outputs progress

    if max(x - xm) < max(x) * 1e-4
        if f_bar; tools.textbar(1); end
        break;
    end
end
%}


end
