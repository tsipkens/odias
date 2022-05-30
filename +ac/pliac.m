
% PLIAC  A modification to the IAC and PLAC methods using a power law/interpolation function.
%  Requires an iterative method to solve. 
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

function [xbar, qbar] = pliac(x_star, nu, q0, eta, c0, n)

%-- Parse inputs ----------------------%
if ~exist('eta', 'var'); eta = []; end
if isempty(eta); eta = 1; end

if isstruct(eta)  % if mass-mobility property struct
    prop = eta;
    c0 = (1 / (prop.k * 1e18)) ^ (1 / prop.zet);
    eta = 1 / prop.zet;
end

if ~exist('c0', 'var'); c0 = []; end
if isempty(c0); c0 = 1; end

if ~exist('n', 'var'); n = []; end
if isempty(n); n = 3.5; end  % controls interpolation function
%-------------------------------------%

pow = eta .* nu;  % combined power

qmodel = @(x, qvec) (1 + (c0 ^ nu * q0 .* x .^ pow) .^ n) .^ (1/n);  % average charge

[xbar, qbar] = ac.iac(x_star, qmodel);

end
