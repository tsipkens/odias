
% GET_POWER_LAW  Get power law parameters for high charge counts.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-30

function [nu, q0, p, X] = get_power_law(qbar0, d)

% Flag higher charge states. 
% Do not fit when average charge is less than 4. 
fl = qbar0 > 4;

A = [log(d(fl)), ones(size(d(fl)))]; % sensitivity matrix for least-squares
b = log(qbar0(fl));  % data = average charge to fit to

p = (A' * A) \ (A' * b);  % least-squares

% Extract power law parameters.
nu = p(1);
q0 = exp(p(2));

X = (1 / q0) ^ (1/nu);

end
