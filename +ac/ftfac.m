
% FTFAC  Full transfer function-average charge algorithm.
%  
%  AUTHOR: Timothy Sipkens, 2022-02-25

function [xbar, qbar] = ftfac(x_star, Kq, q)

% Use second dimension of A for charges. 
if ~exist('q', 'var'); q = []; end
if isempty(q); q = 0:(size(Kq, 2) - 1); end


disp('Running FKAC...');  % add header to console

qbar = sum(sum(Kq, 3) .* q ./ ...
    sum(sum(Kq, 3), 2), 2);

xbar = qbar .* x_star;

tools.textdone();  % mark as complete

end
