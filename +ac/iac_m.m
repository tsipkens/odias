
% IAC_M  Applies the iterative-average-charge algorithm of Corbin/Sipkens et al.
%  Wrapper function for the mass case, which involves a multiple step 
%  particle size coversion.
%  
%  AUTHOR: Timothy Sipkens, 2021-10-06

function [mbar, qbar] = iac_m(m_star, prop, f_deq, cmodel, opt)

if ~exist('f_deq', 'var'); f_deq = []; end
if isempty(f_deq); f_deq = 0; end

if ~exist('cmodel', 'var'); cmodel = []; end
if isempty(cmodel); cmodel = 'Fuchs'; end  % assume unipolar Fuch's model

% The physical model to predict the average charge.
qmodel = @(mi, qvec) physfun(mi, prop, f_deq, qvec, cmodel, opt);

% Call parent IAC function.
[mbar, qbar] = ac.iac(m_star, qmodel);

end


%== PHYSFUN ==============================================================%
%   The physical model to convert mass to average charge. 
function [qbar, di] = physfun(mi, prop, f_deq, qvec, cmodel, opt)

mi = max(mi, 0.00001);

di = (mi .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % use mass-mobility relation
    
% OPTIONAL: Convert to charging equivalent diameter.
if f_deq
    di = working.dm2deq(di);
end

[~, qbar] = kernel.tfer_charge(di * 1e-9, qvec, 298, cmodel, opt);

end


