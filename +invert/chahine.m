
% CHAHINE  An implementation of modified Chahine's method by Ferri et al.
%  
%  Differs from Chahine's original method in that the present method is a
%  vector method that does not perfer coarse discretization. 
%  
%  AUTHOR: Timothy Sipkens, 2024-04-03

function x = chahine(A, b, xi, iter, f_sigma, f_bar)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('iter', 'var'); iter = []; end
if isempty(iter); iter = 100; end % default of 100 iterations

% Intiate vector of ones if xi is not specified.
if ~exist('xi', 'var'); xi = []; end
if isempty(xi);  xi = ones(size(A,2), 1); end
x = xi;

% whether to display text-based progress bar
if ~exist('f_bar', 'var'); f_bar = []; end
if isempty(f_bar); f_bar = 0; end

% controls exit conditions
if ~exist('f_sigma', 'var'); f_sigma = []; end
if isempty(f_sigma); f_sigma = 1; end
%-------------------------------------------------------------------------%


%-- Start evaluation -----------------------------------------------------%
s = 1 ./ max(A, [], 2);  % factor to scale data and model matrix
s_min = max(max(A, [], 2)) * 1e-10;
s(max(A, [], 2) < s_min) = 0;
A = bsxfun(@times, s, A);  % scale model matrix
b = max(b .* s, 0);  % remove any negative data, scale data

if f_bar; disp('Chahine progress:'); tools.textbar([0, iter]); end
    % output textbar for progress


%-- Perform iterations ---------------------------------------------------%
for kk=1:iter
    y = A * x;  % modeled data
    w = A ./ sum(A, 1);  % weight from Ferri et al.
    
    % Avoid division by zero.
    fl1 = y ~= 0;
    fl2 = sum(A, 1) ~= 0;
    x(~fl2) = 0;  % set uninformed data to zero
    
    C = w(fl1, fl2) .* b(fl1) ./ y(fl1);  % calculate weight factor

    x(fl2) = x(fl2) .* sum(C, 1)';  % apply iteration
    
    if f_bar; tools.textbar([kk, iter]); end % outputs progresss

    % Exit conditions if sigma is specified.
    if f_sigma
        if tools.mean_sq_err(A, x, b) < 0.6
            if f_bar
                tools.textbar([iter, iter]);
                disp([' Algorithm converged at ', num2str(kk), ' iterations']);
            end
            break;
        end
    end
end

end
