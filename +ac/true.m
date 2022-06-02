
% TRUE  Computes the transfer function-weighted average transmitted charge.
%  This is not an average charge method as defined by Sipkens et al. in
%  that the average transmitted particle size is computed directly. 
%  
%  NOTE: First argument is a placeholder for X_STA, which is only included 
%  to match the pattern of the other average charge function. However, the
%  value is not used in this method.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-30

function [mbar, mbar_geo] = true(~, K, x, p)

x = x';

% If distribution is ommitted, assume ones.
if ~exist('p', 'var'); p = []; end
if isempty(p); p = ones(size(x));
else; p = p'; end

mbar = sum(p .* K .* x, 2) ./ ...
    sum(p .* K, 2);

mbar_geo = exp(sum(p .* K .* log(x), 2) ./ ...
    sum(p .* K, 2));

end
