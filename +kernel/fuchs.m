
% FUCHS  Evalutes the Fuchs unipolar charging model. 
%  
%  INPUTS: 
%   d    = particle diameter (m)
%   zmax = max. integer charge state
%   T    = temperature (K)
%   P    = pressure (bar)
%   nit  = ion·s/m3
%  
%  OUTPUTS:
%   zeroprob  Probability of no charge.
%   probs     Probability of increasin number of charge states.
%   n_mean    Mean number of charges.
%  
%  ------------------------------------------------------------------------
%  
%  ORIGINAL AUTHOR: Tyler Johnson
%  MODIFIED: Timothy Sipkens, 2021-10-25

function [zeroprob, probs, zmean] = fuchs(d, zmax, T, P, nit, eps)

b = fuchs_sub(zmax, T, d, P, eps);

init_dist = ecd(d, T, zmax); % Equilibrium charge dist (this function is below if you want it – doesn’t make much difference IIRC)

Z = birth_death(b,nit,init_dist);

posdist = Z(end, 3:end);
zeroprob = Z(end, 2);
partition = sum(posdist);
probs = posdist / partition;

% Removes negative probabilities or nonphysical
% results generated by computational rounding.
zmean = 0;
for i=1:1:length(probs)
    if probs(i)<0
        probs(i)=0;
    end
    zmean = zmean + i*probs(i);
end
    
end



%=========================================================================%
%   Subfunctions
%=========================================================================%

%== FUCHS_SUB ============================================================%
%   Determines the combination coefficient b
%   
%   INPUTS:
%    epsilon = Material dielectric constant
%    m = max number of charges to consider (takes longer / blows up with too many)
%    P = pressure (bar)
%    T is temperature in K
function [b] = fuchs_sub(nmax, T, dp, P, epsilon)

a = dp / 2;  % converts particle diameter to radius in meters

%  CONSTANTS
e = 1.6e-19;     % electron charge (C)
k = 1.38e-23;    % Boltzman's constant (J/K)
Na = 6.023e23;   % Avogadro's number (mol^-1)
KE = 9e9;        % prop. constant 
K = (epsilon - 1) / (epsilon + 1);

% ION PROPERTIES
mi = .109;   % ionic molecular weight (kg/mol)
mg = .029;   % air molecular weight (kg/mol)
Z = .00014;  % electrical mobility of ion (m2/Vs)
Zi = Z / P;  % electrical mobility at pressure P
D = k * T * Zi / e;  % diffusion coefficient (m2/s)
ci = sqrt(8 * k * T * Na / pi / mi);  % mean speed of the ions (m/s)
li = 1.329 * Zi / e * sqrt(k * T * mi * mg / (mi + mg) / Na);  % ionic mean free path of gaseous ions


% LIMITING SPHERE RADIUS (Fuchs, 1963, Eqn 5)
delta = (a^3) / (li^2) * ((1/5) * ((1 + li/a)^5) - ...
    (1/3) * (1 + (li^2) / (a^2)) * ((1 + li/a)^3) + ...
    (2/15) * ((1 + (li^2) / (a^2)) ^ (5/2)));


% CALCULATION OF PSI and COLISION PROBABILITY (see Fuchs, 1963)
n = 0:1:(nmax-1);  % charge state array from 0 to maximum considered charge (m)
r = delta:-a/1000:a;  % spacial array from limiting sphere radius to particle radius

% Perform numerical integration.
fun = @(x) int_fun(x, n', a, K, T);  % integration function
PSI = integral(@(x) fun(x), 0, a / delta, 'ArrayValued', true);  % (Fuchs, 1963, Eqn 11)

% Compute potential energies.
pot1 = KE*(e^2) .* ((n' ./ delta) - ...
    K .* ((a^3) / (2 * (delta^2) * ...
    (delta^2 - (a^2)))));  % of the ion at r=limited sphere radius (Fuchs, 1963, Eqn 4)
pot2 = KE .* (e^2) .* ((n' ./ r) - ...
    K .* ((a^3) ./ (2 * (r.^2) .* ...
    (r.^2 - (a^2)))));  % of the ion array calculated across spacial array (Fuchs, 1963, Eqn 4)

% Calculate impactor factor and gam.
impfac = (r.^2) .* (1 + (2/(3*k*T)) .* (pot1 - pot2));  % impactor factor squared (Fuchs, 1963, Eqn 7)
bmsqrd = min(impfac, [], 2);  % finds minimum of squared impactor factor
gam = bmsqrd ./ (delta^2);  % for b<b_ m, alpha=(bm/delta)^2 (Fuchs, 1963, Pg 188)


% Filters for gam.
% Only keeps impactor factors that are between 
% 0 and 1 (Physically possible situations). 
f1 = gam<=1;  % elements that are less than equal to one
f2 = gam>0;  % elements that are greater than equal to zero
f3 = and(f1==1, f2==1);  % combine above filters

% Potential energy of the ion at r=limited sphere radius 
% across charge state array (Fuchs, 1963, Eqn 4). 
pot_d = KE * (e^2) * ((n(f3)' / delta) - K*((a^3) / ...
    (2 * (delta^2) * (delta^2 - (a^2)))));

% Calculation of the combination coefficient.
b = zeros(size(gam));  % zero if f3 is false
b(f3) = ((4*pi*a*D) ./ ((((4*D*a) ./ ...
    (gam(f3) .* ci .* (delta^2))) .* ...
    (exp(pot_d ./ k ./ T))) + PSI(f3)));


end


%== INT_FUN ==============================================================%
% This function calculates the function of PSI integral at 
% the denominator in the limiting sphere model Fuchs (1963).
function [out] = int_fun(x,n,a,K,T)

k = 1.38e-23; % Boltzman's Constant (J/K)
e = 1.6e-19;  % Electron Charge (C)
KE = 9e9;     % Prop. Constant

r = a ./ x;  % normalize distance wrt distance away from particle surface
pot = KE .* (e.^2) * ((n./r) - K .* ((a.^3) ./ ...
    (2 .* (r.^2) .* (r.^2 - (a.^2)))));  % potential energy of the ion at normalized radius (Fuch1963 Eqn 4)
out = exp(pot ./ k ./ T);  % term to be integrated in Fuchs (1963), Eqn 11

% Remove inf values in favour of very large values.
out(isinf(out)) = realmax();

end


%== ECD ==================================================================%
%   Equilibrium charge distribution function.
function f = ecd(d, T, nmax)

k = 1.38e-23; % Boltzmann Constant
e = 1.6e-19;  % electron charge
K = 1 / (4 * pi * 8.85e-12);  % Coloumb law constant
partition = sqrt(pi) / sqrt((K * e ^ 2)/(d * k * T));

n = (1:1:nmax);
f = exp((-K*(e^2) .* (n.^2)) ./ (d*k*T)) ./ partition;

end


%== BIRTH_DEATH ==========================================================%
function [Z] = birth_death(b, nit_total, dist)

% Total charging time (sec).
nit = 0:1e11:nit_total;  % defines nit vector
newm = size(b, 1);

fun = @(t,y) charge(t, y, b, newm);  % defined dy/dt array for ODE
[~, Y] = ode45(@(t,y) fun(t,y), nit, dist');

Z = [nit', Y];

end


%== CHARGE ===============================================================%
%   This function calculates the fraction of charged particles according 
%   (Boisdron and Brock, 1970, Eqn 18) birth-and-death model. The model 
%   can be used with Combination Coefficients (b) determined by any models 
%   available in the literature (cont., fm, and trans. regime). 
function [dy] = charge(~, y, b, zmax)

dy = [-b(1) * y(1); ...
    b(1:(zmax-1)) .* y(1:(zmax-1)) - b(2:end) .* y(2:end)];

end


