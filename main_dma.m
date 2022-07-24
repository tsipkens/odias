
% MAIN  Inversion of mobility distributions.
% AUTHOR: Timothy Sipkens, 2020-04-11
%=========================================================================%

clear;
close all;
addpath cmap;

d = logspace(log10(10), log10(1e3), 500)';  % reconstruction points
d_star = logspace(log10(13.1), log10(763.5), 114)';  % mobility setpoints

prop = kernel.prop_dma;


A = kernel.gen_dma(d_star, d, [], prop);


mu_d = [200, 200/3];
s_d = [1.55, 1.55 / 1.15];
w_d = [1, 0.5];

%{
mu_d = 200;
s_d = 1.2;
w_d = 1;
%}

[b, Lb, x0] = tools.gen_data(A, d, mu_d, s_d, w_d, d_star);


%%

disp(' ');

%-- Least-squares ---------%
disp('Running least-squares ...');
x_lsq = invert.lsq(Lb * A, Lb * b);
tools.textdone();
disp(' ');


%-- Twomey ----------------%
disp('Running Twomey:');
xi = invert.get_init(Lb * A, Lb * b, d, d_star);
x_two = invert.twomey(Lb * A, Lb * b, xi, [], 1, 1);
disp(' ');


%-- Twomey-Markowski ------%
disp('Running Twomey-Markowski:');
xi = invert.get_init(Lb * A, Lb * b, d, d_star);
x_twomark = invert.twomark(Lb * A, Lb * b, length(xi), xi);
disp(' ');


%-- 1st order Tikhonov ----%
disp('Running Tikhonov (1st) ...');
lambda_tk1 = 3.8e1;
[x_tk1, ~, ~, Gpo_inv_tk1] = ...
    invert.tikhonov(Lb * A, Lb * b, lambda_tk1, 1);
Gpo_tk1 = inv(Gpo_inv_tk1);
e.tk1 = (x_tk1 - x0)' * Gpo_inv_tk1 * (x_tk1 - x0);
tools.textdone();
disp(' ');


%-- 2nd order Tikhonov ----%
disp('Running Tikhonov (2nd) ...');
lambda_tk2 = 1e3;
[x_tk2, ~, ~, Gpo_inv_tk2] = ...
    invert.tikhonov(Lb * A, Lb * b, lambda_tk2, 2);
Gpo_tk2 = inv(Gpo_inv_tk2);
e.tk2 = (x_tk2 - x0)' * Gpo_inv_tk2 * (x_tk2 - x0);
tools.textdone();
disp(' ');


%-- Two-step 2nd order Tikhonov --%
disp('Running Tikhonov (2nd, two-step) ...');
lambda_tk2 = 3e3;
[x_tk22, ~, ~, Gpo_inv_tk22] = ...
    invert.tikhonov(Lb * A, Lb * b, lambda_tk2, 2, [], 1);
Gpo_tk22 = inv(Gpo_inv_tk22);
e.tk22 = (x_tk22 - x0)' * Gpo_inv_tk2 * (x_tk22 - x0);
tools.textdone();
disp(' ');


%-- Exponential distance --%
disp('Running exponential distance ...');
lambda_ed = 6e0;
ld = log10(s_d(1));
[x_ed, ~, ~, Gpo_inv_ed] = ...
    invert.exp_dist(Lb * A, Lb * b, lambda_ed, ld, d);
Gpo_ed = inv(Gpo_inv_ed);
e.ed = (x_ed - x0)' * Gpo_inv_ed * (x_ed - x0);
tools.textdone();
disp(' ');
disp(' ');


e



%%
figure(2);

x_tk = x_tk22;
Gpo_tk = Gpo_tk22;

subplot(2, 2, 1);
tools.plotci(d, x_tk, Gpo_tk, x0);
title('Tikhonov (2nd, two-step)');

subplot(2, 2, 2);
tools.plotci(d, x_ed, Gpo_ed, x0);
title('Exponential distance');

subplot(2, 2, 3);
tools.plotci(d, x_two, [], x0);
title('Twomey');
hold on;
plot(d, xi, 'b--');
hold off;

subplot(2, 2, 4);
tools.plotci(d, x_twomark, [], x0);
title('Twomey-Markowski');
hold on;
plot(d, xi, 'b--');
hold off;


%%
% Optimize Tikhonov + show Bayes factor.
%{
[a0, a1, a2] = invert.tikhonov_op(Lb*A,Lb*b,[1e-1,1e3],2);

figure(3);
semilogx([a2.lambda], -[a2.B]);
%}

