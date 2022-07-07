
% GET_TRANSITION  Gets the transition size-to-charge setpoint.
%  
%  AUTHOR: Timothy Sipkens, 2022-07-17

function [xt, X] = get_transition(nu, q0, eta, c0, x_star)

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

pow = nu .* eta;

xt = (c0 .^ nu .* q0) .^ (-1 ./ pow);

if nargout > 1
    X = xt ./ x_star;
end

end
