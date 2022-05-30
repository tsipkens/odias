

clear;
close all;
addpath cmap;


% Set charging model parameters. 
opt.nit = 4e13;
opt.eps = 13.5;


% Set transfer function evaluation grid.
nx = 150;
d = logspace(0.5, 2.25, nx)';  % mobility diameters (as surrogate for deq)

% Get charge distribution.
zmax = 100;
zvec = 0:zmax;
[fq, qbar0] = kernel.tfer_charge(d .* 1e-9, zvec, 298, 'Fuchs', opt);


% Get power law fit.
[nu, q0, p] = ac.get_power_law(qbar0, d);

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
c0 = 1;
eta = 1/3;  % power for deq (e.g., ~1/3 for spheres and x = mass)
x = d .^ (1/eta);

% Assume delta transfer function equally efficient at all sizes.
qmodel = @(x, zvec) out2(@kernel.tfer_charge, x .^ eta .* 1e-9, zvec, 298, 'Fuchs', opt);

[~, q_iac] = ac.iac(x, qmodel);
[~, q_plac] = ac.plac(x, nu, q0, eta, c0);
[~, q_pliac] = ac.pliac(x, nu, q0, eta, c0);

figure(2);
plot(d, qbar0, 'k');  % direct, charge at size x
hold on;
plot(d, qbarh, 'k--');
plot(d, q_iac, 'Color', [0.9, 0.4, 0.4]);  % red, ave. charge for setpoint x*
plot(d, q_plac, 'Color', [0.2, 0.2, 0.8]);  % blue, ave. charge for setpoint x*
plot(d, q_plac, '--', 'Color', [0.1, 0.7, 0.1]);  % green, ave. charge for setpoint x*
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
ylim([0.5, inf]);

function o = out2(fun, varargin)
    [~, o] = fun(varargin{:});
end


