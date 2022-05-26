

clear;
close all;
addpath cmap tfer_pma;


% Defaults.
Rm0 = 3;
nit = 4e13;
eps0 = 13.5;

Dm0 = 2.48;
rho100_0 = 510;
Q = 0.9183;



% Set charging model parameters. 
opt0.nit = nit;
opt0.eps = eps0;


% Set transfer function evaluation grid.
nx = 150;
d = logspace(0, 3, nx)';  % mobility diameters (as surrogate for deq)


zmax = 100;
zvec = 1:zmax;
[fq, qbar] = kernel.tfer_charge(d .* 1e-9, zvec, 298, 'Fuchs', opt0);

% [~, qbar] = kernel.tfer_charge(di * 1e-9, zvec, 298, model, opt);


figure(1);
plot(d, qbar);
set(gca, 'XScale', 'log', 'YScale', 'log')
