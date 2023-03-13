
% WHITE  Evaluate White's (1951) charging model. 
%  
%  AUTHOR: Timothy Sipkens, 2023-02-03

function [qbar1, qbar_pl, A, B] = white(deq, nit)

deq = deq .* 1e-9;  % convert to m

T = 300;

mg = 2 * 2.3258671e-26;

% Use CGS units, as per original work.
e = 4.8e-10;
kB = 1.3807e-16;
nit1 = nit / 1e6;
mg1 = mg * 1e3;
c1 = sqrt(3 * kB * T / mg1);

% Constants in White's equation.
A = 1 ./ 2 .* (kB * T) ./ (e^2) .* 1e2;  % 1e2 accounts for CGS units
B = pi ./ 2 .* c1 .* e^2 .* nit1 ./ (kB * T) .* 1e2;

qbar1 = deq .* A .* log(1 + deq .* B);


qbar_pl = [];  % if output commented below
%{
% For Taylor series expansion of White's. 
deqs = 1e-7;  % i.e., 100 nms
p = B .* deqs ./ ((B .* deqs + 1) .* log(B .* deqs + 1));

k = deqs .* A .* log(1 + deqs .* B);

qbar_pl = k .* (deq ./ (deqs / 1e2)) .^ (1 + p);  % power law approximation
%}

end
