
% EM  Performs expectation-maximization (EM) algorithm. 
%  Assumes Poisson-distributed data.
% 
%  X = invert.em(A, B) solved the inverse problem using EM algorithm, using
%  four iterations and an initial guess of ones. Returns reconstruction X.
%  
%  X = invert.em(A, B, XI) preferred usage adds an initial guess at which
%  the method starts. 
%  
%  X = invert.em(A, B, XI, ITER) explicitly specifies the number of
%  iterations as a scalar. 
%  
%  X = invert.em(A, B, XI, X0) adds a groudn truth, against which the
%  number of iterations is optimized to give the lowest Euvlidean error.
%  The variation X0 must be the same size as X and XI. 
%  
%  NOTE: Adding 
%  
%  AUTHOR:  Timothy Sipkens, 2024-03-22

function [x, ii] = em(A, b, xi, iter_x0)

%-- Parse inputs ---------------------------------------%
if ~exist('xi', 'var'); xi = []; end
if isempty(xi); xi = ones(size(A,2), 1); end

if ~exist('iter_x0', 'var'); iter_x0 = []; end
if isempty(iter_x0); iter_x0 = 4; end

if all(size(iter_x0) == size(xi))
    x0 = iter_x0;
    iter = 100;
else
    iter = iter_x0;
end
%-------------------------------------------------------%

% Initialize x.
x = xi;
x(x==0) = 1e-5; % remove initially zero values

% asj as "A sum j"
asj = sum(A)';
asj = asj + eps;

for ii=1:iter
    x1 = x;  % save copy of x
    
    ax = A * x;
    idx = ax ~= 0; % ignore empty rows of A
    
    % Expectation step.
    y = sum(b(idx) .* A(idx,:) ./ ax(idx), 1);
    
    % Maximization step.
    x = x ./ asj .* y';

    % If x0 is given, check exit condition.
    if exist('x0', 'var')
        if norm(x1 - x0) < norm(x - x0)
            x = x1;  % then revert to previous step
            break;  % exit loop and ultimately the function
        end
    end
end

end

