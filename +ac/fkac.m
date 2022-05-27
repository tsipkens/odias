
% FKAC  Full kernel-average charge algorithm.
%  
%  AUTHOR: Timothy Sipkens, 2022-02-25

function [xbar, qbar] = fkac(x_star, Aq, z)

% Use second dimension of A for charges. 
if ~exist('z', 'var'); z = []; end
if isempty(z); z = 0:(size(Aq, 2) - 1); end


disp('Running FKAC...');  % add header to console

qbar = sum(sum(Aq, 3) .* z ./ ...
    sum(sum(Aq, 3), 2), 2);

xbar = qbar .* x_star;

tools.textdone();  % mark as complete

end
