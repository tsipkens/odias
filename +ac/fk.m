
% FK  Full transfer function-average charge algorithm.
%  Not a true average charge function as this method only ignores the
%  particle size distribution. 
%  
%  NOTE: If input Kq has been multiplied by the size distribution, the true
%  result is output.
%  
%  AUTHOR: Timothy Sipkens, 2022-02-25

function [xbar, qbar] = fk(~, Kq, x, q)

% Use second dimension of A for charges. 
if ~exist('q', 'var'); q = []; end
if isempty(q); q = 0:(size(Kq, 2) - 1); end


disp('Running FK...');  % add header to console

qbar = nansum(nansum(Kq, 3) .* q ./ ...
    nansum(sum(Kq, 3), 2), 2);
xbar = nansum(squeeze(nansum(Kq, 2)) .* x' ./ ...
    nansum(sum(Kq, 3), 2), 2);

% xbar = qbar .* x_star;

tools.textdone();  % mark as complete

end
