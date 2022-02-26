
% FKAC  Full kernel-average charge algorithm.
%  
%  AUTHOR: Timothy Sipkens, 2022-02-25

function [xbar, zbar] = fkac(x_star, kernel)


tools.textheader('Running FKAC');  % add header to console

z = 0:size(kernel, 2);

zbar = exp(sum(kernel .* z' ./ ...
    sum(kernel, 2), 2));

xbar = zbar' .* x_star;

tools.textheader();  % mark as complete

end
