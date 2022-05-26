
% IAC  Applies the iterative-average-charge algorithm of Corbin/Sipkens et l.
%  
%  AUTHOR: Timothy Sipkens, 2021-10-06

function [mi, qbar0, di0] = iac(m_star, prop, f_deq, model, opt)

if ~exist('f_deq', 'var'); f_deq = []; end
if isempty(f_deq); f_deq = 0; end

if ~exist('model', 'var'); model = []; end
if isempty(model); model = 'Fuchs'; end  % assume unipolar Fuch's model

zmax = 100;
zvec = 1:zmax;


% Initialize quantities. 
mi = m_star;
mi1 = 0;
f_conv = zeros(size(mi));  % flag for converge entries, none at the start
qbar0 = f_conv;
di0 = f_conv;


%{
% Alternate form for scalar m_star.
% Generally requires preset and high zmax above
% and makes use of subfunction defined below. 
fun = @(x) x ./ physfun(x, prop, f_deq, zvec, model, opt) - m_star;
mi = fzero(fun, m_star);
[qbar0, di] = physfun(mi, prop, f_deq, zvec, model, opt);  % final calc. for output
%}


%-{
% Vectorized form for vector m_star using Newton's method.
tools.textheader('Running IAC');
n = 50;
for ii=1:n
    
    %== PHYSICAL MODEL ===================================================%
    %   The physical model to predict the average charge.
    
    % First compute the mobility diameter. 
    di = (mi(~f_conv) .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % use mass-mobility relation

    % OPTIONAL: Convert to charging equivalent diameter.
    if f_deq  % check input flag
        di = working.dm2deq(di);
    end
    
    [~, qbar] = kernel.tfer_charge(di * 1e-9, zvec, 298, model, opt);
    %=====================================================================%
    
    
    % If qbar exceed the integration duration.
    % Error check, avoids erroneous truncation of
    % integer charge state during charge model (Fuch's) call.
    f_above = (qbar + sqrt(qbar)) > zmax;  % flag for cases that did not converge due to low zmax
    if any(f_above)
        qbar(f_above) = min(qbar(f_above), zmax);
        zmax = zmax * 1.5;  % increase zmax by 50%
        zvec = 1:zmax;
        
        disp([' Adjusted zmax to ', num2str(zmax)]);
    end
    qbar0(~f_conv) = qbar;  % copy over new values
    di0(~f_conv) = di;
    
    
    % Update mass for unconverged.
    mi(~f_conv) = m_star(~f_conv) .* qbar;
    
    
    % Convergence check/components.
    if all(abs((mi1 - mi) ./ mi1) < 0.001)
        break;
    end
    
    f_conv = abs((mi1 - mi) ./ mi1) < 0.001;  % flag converged entries and skip next time
    mi1 = mi;  % temporarily store previous value (for convergence check)
    
end

%{
% Validation plot, checking charge distribution.
figure(gcf);
plot(1:zmax, pj);
%}

% Check if algorithm did not converge, 
% which occurred if ii = n.
if ii==n
    warning('Algorithm did not converge!');
end
%}

tools.textheader();

end


%{
%== PHYSFUN ==============================================================%
%   Defines the physical model converting mp to qbar.
function [qbar, di] = physfun(mi, prop, f_deq, zvec, model, opt)

mi = max(mi, 0.00001);

% Copy out mass-mobility parameters.
m0 = prop.m0;
Dm = prop.Dm;

di = (mi .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % use mass-mobility relation
    
% OPTIONAL: Convert to charging equivalent diameter.
if f_deq
    di = working.dm2deq(di);
end

[~, qbar] = kernel.tfer_charge(di * 1e-9, zvec, 298, model, opt);

end
%}


