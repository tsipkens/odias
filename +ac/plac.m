
% PLAC  Charge correction using the power law-average charge method.
%  
%  XBAR = ac.plac(X_STAR, NU, Q0) computed the average transmitted particle
%  size at the setpoint, X_STAR, for a power law having an exponent NU
%  and pre-factor Q0. Assumes classification against a particle size
%  equiavlent to charge-equivalent diameter. 
%  
%  XBAR = ac.plac(X_STAR, NU, Q0, PROP) adds an input for a mass-mobility
%  struct, PROP. Computes the average transmitted particle mass. Assumes 
%  the mobility and charge-equivalent diameters are the same. 
%  
%  XBAR = ac.plac(X_STAR, NU, Q0, ETA, C0) allows for explicitly specifying
%  the power law exponent ETA and pre-factor C0 in relating the desired
%  particle size to the charge-equivalent diameter. 
%  
%  [XBAR, QBAR] = ac.plac(...) adds an output for the average transmitted
%  particle charge. 
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

function [xbar, qbar] = plac(x_star, nu, q0, eta, c0)

%-- Parse inputs ----------------------%
if ~exist('eta', 'var'); eta = []; end
if isempty(eta); eta = 1; end

if isstruct(eta)  % if mass-mobility property struct
    c0 = eta.m0;
    eta = eta.zet;
end

if ~exist('c0', 'var'); c0 = []; end
if isempty(c0); c0 = 1; end
%-------------------------------------%

p = eta .* nu;
qbar = (c0 .^ nu .* q0 .* x_star .^ p) ...
    .^ (1 ./ (1 - p));  % average charge

xbar = x_star .* qbar;  % average transmitted particle size

end
