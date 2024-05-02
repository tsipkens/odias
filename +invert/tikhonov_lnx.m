
% TIKHONOV_LNX  Inversion minimizing the second derivative of lnx. 
%  This is a non-linear optimization problem.

function [x, Lpr0] = tikhonov_lnx(A, b, lambda, order, ~, xi)

% If initial guess is not specified.
if ~exist('xi', 'var'); xi = []; end
if isempty(xi); xi = ones(size(A,2), 1); end
xi = xi + eps;

if ~exist('order', 'var'); order = []; end  %  order to Tikhonov
if isempty(order); order = 3; end

% Generate Tikhonov derivative matrix.
Lpr0 = invert.tikhonov_lpr(order, length(xi), 0);

% Solve non-linear least-squares problem.
x = fminsearch(@(x) norm([A * x - b; lambda .* (Lpr0 * log(x))])^2, xi);

end
