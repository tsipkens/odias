
% MAIN_AC_DMA  A function evaluating the average charge algorithms for a DMA.
%  For unipolar-DMA-CPC.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

clear;
close all;
addpath cmap autils;



% Defaults.
Rm = 3;
nit = 4e13; % 4e13; % 1e12
eps = 13.5;

charge_type = 'f';
% charge_type = 'w';  % bipolar


% Set charging model parameters. 
opt.nit = nit;
opt.eps = eps;
opt.n = 2.5;


% Set transfer function evaluation grid.
nx = 1300;

d = logspace(0, log10(4e3), nx)';  % reconstruction points (4e13)
d_star = logspace(0.5, 1.23, 70)';  % mass-to-charge setpoints (4e13)
d_star = logspace(0.5, 1.8, 70)';  % mass-to-charge setpoints (4e13)

% d = logspace(0, log10(2e4), nx)';  % reconstruction points (1e12)
% d_star = logspace(0.5, 1.8, 70)';  % mass-to-charge setpoints (1e12)
% d_star = logspace(0.5, 2, 70)';  % mass-to-charge setpoints (1e12)

a = 3 * pi * 1.82e-5;
ddr = 1 ./ dm2zp(d .* 1e-9, 298, 1) ./ (1e-9 .* a);  % mechanical mobility
ddr_star = 1 ./ dm2zp(d_star .* 1e-9, 298, 1) ./ (1e-9 .* a);  % mechanical mobility
fl_fit = d < 100;
p = polyfit(log(d(fl_fit)), log(ddr(fl_fit)), 1)
ddr_fit = exp(p(2)) .* d .^ 2;


figure(10);
plot(d, ddr);
hold on;
plot(d, d, '--k');
plot(d, ddr_fit);
plot(d, d ./ Cc(d), ':');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');


nb = length(d_star);

z = 0:300;


% Get properties and then update.
prop = kernel.prop_dma;  % use default properties

tools.textheader('Computing kernel')
[K, fq, Kq, qbar0] = kernel.gen_dma(d_star, d, z, [], {charge_type, opt});
tools.textheader();


% Get power law fit.
[nu, q0] = ac.get_power_law(qbar0, d);


figure(2);
h = pcolor(d, z, fq);
set(h, 'EdgeColor', 'none');
set(gca, 'XScale', 'log');



%%
prop0 = struct();
prop0.zet = 1;
prop0.k = 1;

% "True" average transmitted particle size.
d_bar_t = ac.true(d_star, K, d);

[d_bar_fcf, q_bar_fcf, fq_star] = ac.fcf(d_star, fq, d, z);
[d_bar_fk, q_bar_fk] = ac.fk(d_star, Kq, z);

% Run the PLAC algorithm with default settings. 
[d_bar_plac, q_bar_plac] = ...
    ac.plac(d_star, nu, q0, prop0);
[d_bar_intac, q_bar_intac] = ...
    ac.intac(d_star, nu, q0, prop0, [], opt.n);  % instead solved with IAC

% Assymptotic approach.
a1 = 0.2;
a2 = 2;
a3 = 1 ./ q0;

fl_fit = d_star < (a3 .* 0.8);
d_bar_fun = @(a, fl) d_star(fl) .* ...
    (a(1) .* d_star(fl) .^ a(2) ./ (a3 - d_star(fl)) .^ a(2) + 1);
a0 = lsqnonlin(@(a) log(d_bar_fun(a, fl_fit)) - log(d_bar_fcf(fl_fit)), [a1; a2]);
d_bar_vintac = d_bar_fun(a0, true(size(fl_fit)));


figure(3);
imagesc(d_star, z, fq_star');
set(gca, 'YDir', 'normal');
set(gca, 'XScale', 'log');

figure(4);
imagesc(d_star, z, nansum(Kq, 3)');
set(gca, 'YDir', 'normal');
set(gca, 'XScale', 'log');


% Plot errors at default.
figure(1); clf;
plot(d_star, d_bar_intac, 'k-');
hold on;
plot(d_star, d_bar_plac, '-');
plot(d_star, d_bar_fcf, 'ro-');
plot(d_star, d_bar_fk, 'bo-');
plot(d_star, d_bar_vintac, 'g-');
plot(d_star, d_bar_t, 'ko-');
plot(d_star, d_star, 'k--');
xline(1 ./ q0);
xline(1 ./ q0 ./ 2, ':');
xline(0.75 ./ q0, '--');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');

ylim([d_star(1), 800]);
xlim([-inf, 100]);





%%
% Add classification with respect to drag force.

bet = 1 / 2;
c0 = sqrt(1 / exp(p(2)));

% Run the PLAC algorithm with default settings. 
[Dr_bar_plac, q_bar_plac] = ...
    ac.plac(ddr_star, nu, q0, bet, c0);
[Dr_bar_intac, q_bar_intac] = ...
    ac.intac(ddr_star, nu, q0, bet, c0, opt.n);  % instead solved with IAC

[Dr_bar_fcf, q_bar_fcf, fq_star] = ac.fcf(ddr_star, fq, ddr, z);
[Dr_bar_fk, q_bar_fk] = ac.fk(ddr_star, Kq, z);


figure(11); clf;
plot(ddr_star, Dr_bar_intac, 'k-');
hold on;
plot(ddr_star, Dr_bar_plac, '-');
plot(ddr_star, Dr_bar_fcf, 'ro-');
plot(ddr_star, Dr_bar_fk, 'bo');
% plot(Dr_star, Dr_bar_t, 'ko-');
% plot(Dr_star, Dr_star, 'k--');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');


figure(1);
alt = (sqrt(q0) .* d_star) .^ (1 ./ (1 - nu/2));
hold on;
plot(d_star, alt, '-');
plot(d_star, (d_star .^ opt.n + alt .^ opt.n) .^ (1/opt.n), '-');
hold off;

