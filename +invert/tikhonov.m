
% TIKHONOV  Performs inversion using various order Tikhonov regularization in 2D.
%
%  INPUTSs:
%   A           Model matrix
%   b           Data
%   order_L     Order of regularization -OR- 
%               pre-computed Tikhonov matrix structure
%                   (OPTIONAL, default is set by invert.tikhonov_lpr)
%   xi          Initial guess for solver
%                   (OPTIONAL, default is zeros)
%   f2          Flag of whether to use two-stage procedure of 
%               Huckle and Sedlacek.
%   solver      Solver (OPTIONAL, default is 'interior-point')
%
%  OUTPUTS:
%   x           Regularized estimate
%   D           Inverse operator (x = D*[b;0])
%   Lpr0        Tikhonov matrix structure
%   Gpo_inv     Inverse of posterior covariance
%  
%  AUTHOR: Timothy Sipkens, 2020-04-11

function [x,D,Lpr0,Gpo_inv] = tikhonov(A, b, lambda, order_L, xi, f2, solver)

n = size(A,2);  % length of x

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order_L', 'var'); order_L = []; end
    % if order not specified, use default of tikhonov_lpr

if ~exist('xi', 'var'); xi = []; end % if initial guess is not specified
if ~exist('solver', 'var'); solver = []; end

if ~exist('f2', 'var'); f2 = []; end
if isempty(f2); f2 = []; end
%-------------------------------------------------------------------------%


%-- Get Tikhonov smoothing matrix ----------------------------------------%
if all(size(order_L)==[1,1]) % if order is specified, build Lpr0
    Lpr0 = invert.tikhonov_lpr(...
        order_L, n);
else % is Lpr0 strucutre is provided directly
    Lpr0 = order_L;
end
Lpr0 = Lpr0;
Lpr = lambda .* Lpr0;


%-- Choose and execute solver --------------------------------------------%
pr_length = size(Lpr0,1);
[x,D] = invert.lsq(...
    [A;Lpr], [b;sparse(pr_length,1)], ...
    xi, solver);

% If flagged, apply two-step procedure.
if f2
    D2 = sparse(1:n, 1:n, 1 ./ max(abs(x), 0.01), n, n);
    [x,D] = invert.lsq(...
        [A;Lpr*D2], [b;sparse(pr_length,1)], ...
        xi, solver);
end


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end
%=========================================================================%

