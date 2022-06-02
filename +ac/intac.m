
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

%-{
% Use IAC method.
qmodel = @(x, qvec) (1 + (c0 ^ nu * q0 .* x .^ pow) .^ n) .^ (1/n);  % average charge
[xbar, qbar] = ac.iac(x_star, qmodel);
%}

%{
% Use specifically-derived Newton's method.
a = c0 ^ nu * q0 .* x_star .^ pow;
qmodel = @(qbar) qbar - (1 + (a .* qbar .^ pow) .^ n) .^ (1/n);
dqmodel = @(qbar) abs(1 - pow .* a .^ n .* qbar .^ (n*pow - 1) .* ...
    (1 + (a .* qbar .^ pow) .^ n) .^ (1/n - 1));
qbar1 = ones(size(x_star));
qbar0 = zeros(size(x_star));
while any((qbar1 - qbar0) > 0.01)
    qbar0 = qbar1;
    qbar1 = qbar0 - qmodel(qbar0) ./ abs(dqmodel(qbar0));
end

qbar = qbar1;
xbar = qbar .* x_star;
%}

end
