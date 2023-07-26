
% EXP_DIST  Regularization based on the exponential of the distance between elements/pixels.
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

function [x,D,Lpr0,Gpo_inv] = exp_dist(A, b, lambda, ld, vec, xi, solver)

x_length = size(A,2);

%-- Parse inputs -----------------------------------------------%
if ~exist('xi','var'); xi = []; end % if no initial x is given
if ~exist('solver','var'); solver = []; end % if computation method not specified

if ~exist('vec', 'var'); vec = []; end
if and(numel(ld) == 1, ~isempty(vec))
    Lpr0 = invert.exp_dist_lpr(ld, vec); % use external function to evaluate prior covariance
else
    Lpr0 = ld;
end
%--------------------------------------------------------------%

% Scael Lpr by regularization parameter.    
Lpr = lambda .* Lpr0;


%-- Choose and execute solver --------------------------------------------%
[x,D] = invert.lsq(...
    [A;Lpr], [b;sparse(x_length,1)], xi, solver);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end
