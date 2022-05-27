

clear;
close all;
addpath cmap;


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

% Get charge distribution.
zmax = 100;
zvec = 0:zmax;
[fq, qbar0] = kernel.tfer_charge(d .* 1e-9, zvec, 298, 'Fuchs', opt0);

% Get power law fit.
fl = qbar0 > 4;  % flag higher charge states
A = [log(d(fl)), ones(size(d(fl)))];
b = log(qbar0(fl));
p = (A' * A) \ (A' * b);
nu = p(1)
q0 = exp(p(2))
% qbarh = exp(polyval(p, log(d)));
% qbarh = exp(nu .* log(d) + log(q0));
qbarh = q0 .* d .^ nu;

qbarl = ones(size(d));

n = 3;
qbart = (qbarh .^ n + 1) .^ (1/n);




% Plot charge information.
figure(1);
plot(d, qbarh, 'b');
hold on;
plot(d, qbart, 'b');
plot(d, qbar0, 'k');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
ylim([0.5, inf]);


%{
% Plot distribution of charges.
figure(3);
h = pcolor(d, zvec, fq);
set(gca, 'XScale', 'log', 'YScale', 'log');
set(h, 'EdgeColor', 'none');
colormap(flipud(ocean));

hold on;  % add qbar0 computed above
plot(d, qbar0, 'k');
hold off;
%}


%%
%== Basic average charge ==============================%
%   i.e., perfect classification using charge-equivalent diameter.
qmodel = @(d, zvec) out2(@kernel.tfer_charge, d .* 1e-9, zvec, 298, 'Fuchs', opt0);

[~, q_fkac] = ac.fkac(d, fq');
[~, q_iac] = ac.iac(d, qmodel);
[~, q_plac] = ac.plac(d, nu, q0);

figure(1);
hold on;
plot(d, q_fkac, 'r--');
plot(d, q_iac, 'g:');
plot(d, q_plac, 'c');
hold off;



function o = out2(fun, varargin)
    [~, o] = fun(varargin{:});
end


