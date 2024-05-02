
% TIKHONOV_OP  Finds optimal lambda for Tikhonov solver using known distribution, x.
% 
%  INPUTS:
%   A       Model matrix
%   b       Data
%   n_grid  Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   order   Order of regularization     (Optional, default is 1)
%   n       Number of lambda entries    (Optional, default is 70)
%
%  OUTPUTS:
%   x       Regularized estimate
%   lambda  Semi-optimal regularization parameter
%           (against exact solution if x_ex is specified or using Bayes factor)
%   output  Output structure with information for a range of the regularization parameter
%   Gpo_inv     Inverse of posterior covariance
%  
%  AUTHOR: Timothy Sipkens, 2019-07-17

function [x, lambda, Lpr, Gpo_inv] = tikhonov_op(A, b, x0, lam0, order_L, bc, xi, solver)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order_L', 'var'); order_L = []; end
if ~exist('bc', 'var'); bc = []; end
if ~exist('xi', 'var'); xi = []; end
if ~exist('solver', 'var'); solver = []; end

if ~exist('lam0', 'var'); lam0 = []; end
if isempty(lam0); lam0 = 1e3; end
%-------------------------------------------------------------------------%

lambda = fminsearch(@(lam) ...
    norm(x0 - invert.tikhonov(A, b, lam, order_L, bc, xi, solver))^2, lam0);

[x, ~, Lpr, Gpo_inv] = invert.tikhonov(A, b, lambda, order_L, bc, xi, solver);

end

