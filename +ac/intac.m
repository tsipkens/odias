
% INTAC  Interpolation-average charge method.
%  Requires an iterative method to solve. 
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

function [xbar, qbar] = intac(x_star, nu, q0, eta, c0, n)

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

%{
% Use IAC method.
qmodel = @(x, qvec) (1 + (c0 ^ nu * q0 .* x .^ pow) .^ n) .^ (1/n);  % average charge
[xbar, qbar] = ac.iac(x_star, qmodel);
%}

%-{
% Use specifically-derived Newton's method.
a = (c0 ^ nu * q0 .* x_star .^ pow) .^ n;
qmodel = @(x) -x + (a .* x .^ (pow*n) + 1) .^ (1/n);
dqmodel = @(x) -1 + (1/n) .* (a .* x .^ (pow*n) + 1) .^ (1/n - 1) .* ...
    (a .* (pow*n) .* x .^ (pow*n - 1));
updater = @(x) x - qmodel(x) ./ dqmodel(x);
y1 = c0 ^ nu * q0 .* x_star .^ (pow);
fl = ones(size(y1));
while any(fl)
    y0 = y1;
    y1 = updater(y0);
    
    fl = abs(y0 - y1) > 0.001;
end
qbar = y1;
xbar = qbar .* x_star;
%}

end
