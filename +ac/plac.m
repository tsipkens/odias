
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

function [xbar, qbar, qfun] = plac(x_star, nu, q0, eta, c0)

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
%-------------------------------------%

disp('Running PLAC...');  % add header to console

p = eta .* nu;  % combined power

qbar = (c0 .^ nu .* q0) .^ (1 ./ (1 - p)) .* ...
    x_star .^ (p ./ (1 - p));  % average charge
xbar = qbar .* x_star;  % average transmitted particle size

if p > 0.8
    warning('Power law average charge will be inaccurate. Returning NaN.');
    qbar = NaN(size(qbar));
    xbar = NaN(size(xbar));
end

qlim = 1.2;
if any(qbar < qlim)
    warning(['Average charges below ', num2str(qlim), ' truncated.']);
    qbar(qbar < qlim) = NaN;
end

% Additional output: a function handle for qbar.
qfun = @(x_star) (c0 .^ nu .* q0 .* x_star) ...
    .^ (1 ./ (1 - p)) ./ x_star;

tools.textdone();  % mark as complete
disp(' ');

end
