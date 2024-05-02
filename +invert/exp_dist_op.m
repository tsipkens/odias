
% EXP_DIST_OP  Regularization based on the exponential of the distance between elements/pixels.
%  
%  X = invert.exp_dist(A, B, LAMBDA, LPR) performs exponential distance
%  regularization of the basic problem AX = B, solving for X. The quantity
%  LAMBDA is a regularization parameter, and pre-computed LPR is the 
%  structure of the exponential distance prio matrix. 
%  See invert.exp_dist_lpr(...) for more information on LPR. 
%  
%  X = invert.exp_dist(A, B, LAMBDA, LD, VEC) replaces the pre-computed
%  matrix with a length scale, LD, and vector at which the size
%  distribution is to be reconstructed. LPR is then built in the function. 
%  
%  X = invert.exp_dist(A, B, LAMBDA, LPR, [], XI) adds an initial X for
%  iterative solvers. {LPR, []} can be replaced with {LD, VEC}. 
%  
%  X = invert.exp_dist(..., SOLVER) adds a string to specify the type of
%  solver to use. 
%  
%  ------------------------------------------------------------------------
%
%  OUTPUTS:
%   X           Regularized estimate
%   D           Inverse operator (x = D*[b;0])
%   LPR0        Cholesky factorization of prior covariance
%   GPO_INV     Inverse of the posterior covariance after inversion
%  
%  ------------------------------------------------------------------------
%  
%  Author: Timothy Sipkens, 2018-10-22

function [x, lambda, ld, Lpr, Gpo_inv] = exp_dist_op(A, b, x0, lam0, ld0, vec, xi, solver)

%-- Parse inputs -----------------------------------------------%
if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('solver','var'); solver = []; end % if computation method not specified
if ~exist('vec','var'); vec = []; end % if computation method not specified

% Optimize solution.
t = fminsearch(@(t) ...
    norm(x0 - invert.exp_dist(A, b, t(1), t(2), vec, xi, solver))^2, [lam0, ld0]);
lambda = t(1);
ld = t(2);

[x, ~, Lpr, Gpo_inv] = invert.exp_dist(A, b, lambda, ld, vec, xi, solver);

end
