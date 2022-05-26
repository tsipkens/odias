
% FKAC  Full kernel-average charge algorithm.
%  
%  AUTHOR: Timothy Sipkens, 2022-02-25

function [xbar, qbar] = fkac(x_star, A)

disp('Running FKAC...');  % add header to console

z = 0:(size(A, 2) - 1);

qbar = sum(A .* z ./ ...
    sum(A, 2), 2);

xbar = qbar' .* x_star;

tools.textdone();  % mark as complete

end
