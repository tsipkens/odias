
% MAIN  Inversion of mobility distributions.
% AUTHOR: Timothy Sipkens, 2020-04-11
%=========================================================================%

clear;
close all;
addpath cmap tfer_pma;

m = logspace(-3, 2, 500)';  % reconstruction points
m_star = logspace(-3, 2, 80)';  % mass-to-charge setpoints

prop = kernel.prop_pma;
d = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % get mobility diameters

sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', 3);
A = kernel.gen_pma(sp, m, d, 1:100, prop, [], 'Fuchs');
Ab = kernel.gen_pma(sp, m, d, 1:100, prop);

mu = [1, 0.1];
s = [2.5, 1.9];
w = [1, 0.5];

%-{
mu = 1;
s = 2.5;
w = 1;
%}


[b, Lb, x0] = tools.gen_data(A, m, mu, s, w, m_star);

hold on;
tools.gen_data(Ab, m, mu, s, w, m_star);
hold off;

%%

disp(' ');

%-- Least-squares ---------%
disp('Running least-squares ...');
x_lsq = invert.lsq(Lb * A, Lb * b);
tools.textdone();
disp(' ');


%-- Twomey ----------------%
disp('Running Twomey:');
xi = invert.get_init(Lb * A, Lb * b, m, m_star);
x_two = invert.twomey(Lb * A, Lb * b, xi, [], 1, 1);
disp(' ');


%-- Twomey-Markowski ------%
disp('Running Twomey-Markowski:');
xi = invert.get_init(Lb * A, Lb * b, m, m_star);
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
lambda_ed = 0.9e1;
ld = log10(s(1));
[x_ed, ~, ~, Gpo_inv_ed] = ...
    invert.exp_dist(Lb * A, Lb * b, lambda_ed, ld, m);
Gpo_ed = inv(Gpo_inv_ed);
e.ed = (x_ed - x0)' * Gpo_inv_ed * (x_ed - x0);
tools.textdone();
disp(' ');
disp(' ');


e



%%
figure(2);

x_tk = x_tk2;
Gpo_tk = Gpo_tk2;

subplot(2, 2, 1);
tools.plot_ci(m, x_tk, Gpo_tk, x0);
title('Tikhonov (2nd)');

subplot(2, 2, 2);
tools.plot_ci(m, x_ed, Gpo_ed, x0);
title('Exponential distance');

subplot(2, 2, 3);
tools.plot_ci(m, x_two, [], x0);
title('Twomey');
hold on;
plot(m, xi, 'b--');
hold off;

subplot(2, 2, 4);
tools.plot_ci(m, x_twomark, [], x0);
title('Twomey-Markowski');
hold on;
plot(m, xi, 'b--');
hold off;


%%
% Optimize Tikhonov + show Bayes factor.
%{
[a0, a1, a2] = invert.tikhonov_op(Lb*A,Lb*b,[1e-1,1e3],2);

figure(3);
semilogx([a2.lambda], -[a2.B]);
%}

