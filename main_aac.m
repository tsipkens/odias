
clear;
close all;
addpath cmap tfer;

da = logspace(log10(10), log10(1e3), 500)';  % reconstruction points
da_star = logspace(log10(13.1), log10(700), 114);  % mobility setpoints
% 763.5

prop = prop_aac();
prop = massmob.add(prop, 'soot')


A = tfer_aac(da_star, da, prop)';

% Generate distribution.
mu_d = 120;
s_d = 1.7;
x0 = normpdf(log(da), log(mu_d), log(s_d));

b0 = A * x0;
[b, Lb] = tools.get_noise(b0, 1e2, 1e-6);


figure(1);
semilogx(da_star, b, '.');
hold on;
semilogx(da_star, b0, 'Color', [0.5,0.5,0.5]);
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
xi = invert.get_init(Lb * A, Lb * b, da, da_star);
x_two = invert.twomey(Lb * A, Lb * b, xi, [], 1, 1);
disp(' ');


%-- Twomey-Markowski ------%
disp('Running Twomey-Markowski:');
xi = invert.get_init(Lb * A, Lb * b, da, da_star);
x_twomark = invert.twomark(Lb * A, Lb * b, length(xi), xi);
disp(' ');


%-- 1st order Tikhonov ----%
disp('Running Tikhonov (1st) ...');
lambda_tk1 = 3.8e1;
[x_tk1, ~, ~, Gpo_inv_tk1] = ...
    invert.tikhonov(Lb * A, Lb * b, lambda_tk1, 1, 0);
Gpo_tk1 = inv(Gpo_inv_tk1);
e.tk1 = (x_tk1 - x0)' * Gpo_inv_tk1 * (x_tk1 - x0);
tools.textdone();
disp(' ');


%-- 2nd order Tikhonov ----%
disp('Running Tikhonov (2nd) ...');
lambda_tk2 = 1e3;
[x_tk2, ~, ~, Gpo_inv_tk2] = ...
    invert.tikhonov(Lb * A, Lb * b, lambda_tk2, 2, 0);
Gpo_tk2 = inv(Gpo_inv_tk2);
e.tk2 = (x_tk2 - x0)' * Gpo_inv_tk2 * (x_tk2 - x0);
tools.textdone();
disp(' ');


%-- Two-step 2nd order Tikhonov --%
disp('Running Tikhonov (2nd, two-step) ...');
lambda_tk2 = 3e3;
[x_tk22, ~, ~, Gpo_inv_tk22] = ...
    invert.tikhonov(Lb * A, Lb * b, lambda_tk2, 2, 0);
Gpo_tk22 = inv(Gpo_inv_tk22);
e.tk22 = (x_tk22 - x0)' * Gpo_inv_tk2 * (x_tk22 - x0);
tools.textdone();
disp(' ');


%-- Exponential distance --%
disp('Running exponential distance ...');
lambda_ed = 1e1;
ld = log10(s_d(1));
[x_ed, ~, ~, Gpo_inv_ed] = ...
    invert.exp_dist(Lb * A, Lb * b, lambda_ed, ld, da);
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

subplot(2, 3, 1);
tools.plotci(da, x_tk1, Gpo_tk1, x0);
title('Tikhonov (1st)');

subplot(2, 3, 2);
tools.plotci(da, x_tk2, Gpo_tk2, x0);
title('Tikhonov (2nd)');

subplot(2, 3, 3);
tools.plotci(da, x_tk, Gpo_tk, x0);
title('Tikhonov (2nd, two-step)');

subplot(2, 3, 4);
tools.plotci(da, x_ed, Gpo_ed, x0);
title('Exponential distance');

subplot(2, 3, 5);
tools.plotci(da, x_two, [], x0);
title('Twomey');
hold on;
plot(da, xi, 'b--');
hold off;

subplot(2, 3, 6);
tools.plotci(da, x_twomark, [], x0);
title('Twomey-Markowski');
hold on;
plot(da, xi, 'b--');
hold off;


%%
% Optimize Tikhonov + show Bayes factor.
%{
[a0, a1, a2] = invert.tikhonov_op(Lb*A,Lb*b,[1e-1,1e3],2);

figure(3);
semilogx([a2.lambda], -[a2.B]);
%}

