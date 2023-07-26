
% INTAC  Interpolation-average charge method.
%  Requires an iterative method to solve. 
%  
%  [XBAR, QBAR] = intac(X_STAR, NU, Q0, ...) computes the average
%  transmitted particle size/mass (XBAR) and particle charge (QBAR) at the
%  size/mass-to-charge setpoint X_STAR and where the charging model is
%  represented with a power law having a prefactor of Q0 and exponent of
%  NU. Also requires either a struct, provided as BET, or both BET and C0,
%  as per below. By default, uses an interpolating power of N = 2.5. 
%  
%  [XBAR, QBAR] = intac(X_STAR, NU, Q0, BET, C0) adds an exponent BET and
%  prefactor C0 relating the particle mobility diameter to the
%  single-particle mass. Not that BET is 1/ZET, where ZET is the
%  mass-mobility exponent (also denoted as Dm in some works).
%  
%  [XBAR, QBAR] = intac(X_STAR, NU, Q0, BET) uses the struct BET, which
%  must have fields BET.zet, which is the mass-mobility exponent  and is 
%  equal to 1/BET in the other variant, and BET.k, which is the 
%  mass-mobility prefactor. 
%  
%  [XBAR, QBAR] = intac(..., N) adds an input for the interpolating power.
%  By default, N = 2.5. 
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

function [xbar, qbar] = intac(x_star, nu, q0, bet, c0, n)

%-- Parse inputs ----------------------%
if ~exist('bet', 'var'); bet = []; end
if isempty(bet); bet = 1; end

if isstruct(bet)  % if mass-mobility property struct
    prop = bet;
    c0 = (1 / (prop.k * 1e18)) ^ (1 / prop.zet);
    bet = 1 / prop.zet;
end

if ~exist('c0', 'var'); c0 = []; end
if isempty(c0); c0 = 1; end

if ~exist('n', 'var'); n = []; end
if isempty(n); n = 2.5; end  % controls interpolation function
%-------------------------------------%

disp('Running INTAC...');  % add header to console

pow = bet .* nu;  % combined power

%{
% Use IAC method.
qmodel = @(x, qvec) (1 + (c0 ^ nu * q0 .* x .^ pow) .^ n) .^ (1/n);  % average charge
[xbar, qbar] = ac.iac(x_star, qmodel);
%}

%{
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

xbar = (x_star .^ n + (c0 .^ nu .* q0 .* x_star) .^ (n ./ (1 - pow))) .^ (1 ./ n);
qbar = xbar ./ x_star;

tools.textdone();  % mark as complete
disp(' ');

end
