
% FCF  Full charge fraction algorithm for computing average charge.
%  
%  NOTE: Not a true "average charge" function as this method only ignores the
%  particle size distribution and transfer function. 
%  
%  AUTHOR: Timothy Sipkens, 2022-07-05

function [xbar, qbar, fq_star] = fcf(x_star, fq, x, q)

% Use second dimension of A for charges. 
if ~exist('q', 'var'); q = []; end
if isempty(q); q = 0:(size(fq, 2) - 1); end


disp('Computing average charge using full charge fraction (FCF)...');  % add header to console

% Compute charge fractions at qx*.
fq_star = interp2(x, q, fq, ...
    x_star .* q, ones(size(x_star)) .* q);

% Get average charge. 
% nansum's ignore fq where extrapolation would be required above.  
qbar = nansum(fq_star .* q ./ ...
    nansum(fq_star, 2), 2);

qbar(nansum(fq_star, 2) < 1) = NaN;  % interpolation failed at these points

xbar = qbar .* x_star;

tools.textdone();  % mark as complete
disp(' ');

end
