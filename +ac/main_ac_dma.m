
% MAIN_AC_DMA  A function evaluating the average charge algorithms for a DMA.
%  For unipolar-DMA-CPC.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

clear;
close all;
addpath cmap autils;



% Defaults.
Rm = 3;
nit = 4e13;
eps = 13.5;

charge_type = 'Fuchs';


% Set charging model parameters. 
opt.nit = nit;
opt.eps = eps;
opt.n = 2.5;


% Set transfer function evaluation grid.
nx = 1200;
d = logspace(0, log10(4e3), nx)';  % reconstruction points
d_star = logspace(0.5, 2, 10)';  % mass-to-charge setpoints
nb = length(d_star);

z = 1:300;


% Get properties and then update.
prop = kernel.prop_dma;  % use default properties

tools.textheader('Computing kernel')
[K, fq, Kq, qbar0] = kernel.gen_dma(d_star, d, z, [], {'Fuchs', opt});
tools.textheader();


% Get power law fit.
[nu, q0] = ac.get_power_law(qbar0, d);



%%
prop0 = struct();
prop0.zet = 1;
prop0.k = 1;

% "True" average transmitted particle size.
d_bar_t = ac.true(d_star, K, d);
d_bar_t(find(d_bar_t > 0.4 * d(end)):end) = NaN;

[d_bar_fcf, q_bar_fcf] = ac.fcf(d_star, fq, d, z);
[d_bar_ftf, q_bar_ftf] = ac.ftf(d_star, Kq, z);

% Run the PLAC algorithm with default settings. 
[d_bar_plac, q_bar_plac] = ...
    ac.plac(d_star, nu, q0, prop0);
[d_bar_intac, q_bar_intac] = ...
    ac.intac(d_star, nu, q0, prop0, opt.n);  % instead solved with IAC


% Plot errors at default.
figure(1); clf;
plot(d_star, d_bar_intac, 'o-');
hold on;
plot(d_star, d_bar_plac, 'o-');
plot(d_star, d_bar_fcf, 'ro-');
plot(d_star, d_bar_ftf, 'bo-');
plot(d_star, d_bar_t, 'ko-');
plot(d_star, d_star, 'k--');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
ylim([d_star(1), 1000]);



