
% GEN_ELPI Compute kernel/transfer function for ELPI+.
%  Uses the fits from Järvinen et al. (2014).
%  
%  K = kernel.gen_elpi(DA) computes the kernel, K, for the aerodynamic
%  dimaeters DA in nanometers. 
%  
%  AUTHOR: Timothy Sipkens, 2021-04-05

function [K, d50] = gen_elpi(da)

% Values from Järvinen et al. (2014).
d50 = [15.7, 30.4, 54.1, 94.3, 154., ...
       254., 380., 600., 943., 1620, ...
       2460, 3640, 5340];
s50 = [3.32, 3.65, 3.89, 3.05, 3.62, ...
       6.30, 8.43, 7.16, 6.21, 5.32, ...
       5.33, 4.14, 3.66];

En = @(da, d50, s50) 1 ./ (1 + (d50 ./ da) .^ (2 .* s50));


% Generate kernel.
K = zeros(length(d50), length(da));
B = zeros(1, length(da));
for ii=length(d50):-1:1
    t0 = En(da, d50(ii), s50(ii));
    
    K(ii, :) = t0 - B;
    
    B = 1 - (1 - B) .* (1 - t0);
end


end

