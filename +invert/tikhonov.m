
% TIKHONOV  Performs inversion using various order Tikhonov regularization in 2D.
%
%  INPUTS:
%   A           Model matrix
%   b           Data
%   order_L     Order of regularization -OR- 
%               pre-computed Tikhonov matrix structure
%                   (OPTIONAL, default is set by invert.tikhonov_lpr)
%   xi          Initial guess for solver
%                   (OPTIONAL, default is zeros)
%   solver      Solver (OPTIONAL, default is 'interior-point')
%
%  OUTPUTS:
%   x           Regularized estimate
%   D           Inverse operator (x = D*[b;0])
%   Lpr0        Tikhonov matrix structure
%   Gpo_inv     Inverse of posterior covariance
%  
%  AUTHOR: Timothy Sipkens, 2020-04-11

function [x,D,Lpr0,Gpo_inv] = tikhonov(A, b, lambda, order_L, bc, xi, solver)

n = size(A,2);  % length of x

%-- Parse inputs ---------------------------------------------------------%
% If order not specified, will use default of tikhonov_lpr.
if ~exist('order_L', 'var'); order_L = []; end

if ~exist('xi', 'var'); xi = []; end % if initial guess is not specified
if ~exist('solver', 'var'); solver = []; end

% If boundary conditions not specified, will use default of tikhonov_lpr.
if ~exist('bc', 'var'); bc = []; end

% Determine number of stages.
if any(size(order_L) == [1,1])  % determine if orders were given
    n_stages = max(size(order_L));
    order = order_L;
else
    n_stages = 1;
    order = [];  % no orders were given, set to empty
end
n_stages = max([n_stages, ...
    max(size(lambda)), max(size(bc))]);  % check other inputs

% Fill out based on number of stages if inputs are scalars.
if numel(lambda) == 1; lambda = lambda .* ones(n_stages, 1); end
if numel(order) == 1; order = order .* ones(n_stages, 1); end
if numel(bc) == 1; bc = bc .* ones(n_stages, 1); end
%-------------------------------------------------------------------------%


%-- Get Tikhonov smoothing matrix ----------------------------------------%
% OPTION 1: if order is specified, build Lpr0.
if ~isempty(order)
    if isempty(bc); Lpr0 = invert.tikhonov_lpr(order(1), n);
    else; Lpr0 = invert.tikhonov_lpr(order(1), n, bc(1));
    end

% OPTION 2: Lpr0 strucutre is provided directly.
else; Lpr0 = order_L;
end

Lpr = lambda(1) .* Lpr0;  % incorporate regularization parameter


%-- Choose and execute solver --------------------------------------------%
pr_length = size(Lpr0,1);
[x,D] = invert.lsq(...
    [A; Lpr], [b; sparse(pr_length,1)], ...
    xi, solver);


%== Mutliple stage Tikhonov ==============================================%
% Depending on inputs, apply 2+ stage procedure.
% See also Huckle and Sedlacek for 2 stage procedure.
if n_stages > 1
    for ii=2:n_stages
        %-- Get Tikhonov smoothing matrix --------------------------------%
        % OPTION 1: if order is specified, build Lpr0.
        if ~isempty(order)
            if isempty(bc); Lpr0 = invert.tikhonov_lpr(order(ii), n);
            else; Lpr0 = invert.tikhonov_lpr(order(ii), n, bc(ii));
            end
        
        % OPTION 2: Lpr0 strucutre is provided directly.
        else; Lpr0 = order_L;
        end
        
        Lpr = lambda(ii) .* Lpr0;  % incorporate regularization parameter
        
        %-- Choose and execute solver ------------------------------------%
        D2 = sparse(1:n, 1:n, 1 ./ max(abs(x), 0.01), n, n);
        [x,D] = invert.lsq(...
            [A; Lpr * D2], [b; sparse(pr_length,1)], ...
            xi, solver);
    end
end
%=========================================================================%


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A + Lpr'*Lpr;
end

end

