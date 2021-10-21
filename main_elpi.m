
close all;
addpath cmap;

d = logspace(log10(10), log10(1e5), 700)';  % reconstruction points

[A, d_star] = kernel.gen_elpi(d');  % kernel for ELPI+



mu_d = 1000;
s_d = 1.8;

% Bimodal.
x0 = normpdf(log10(d), log10(mu_d), log10(s_d)) + ...
    0.5 .* normpdf(log10(d), log10(mu_d / 5), log10(s_d / 1.15));

b0 = A * x0;

[b, Lb] = tools.get_noise(b0, 1e2, 1e-6);


figure(1);
semilogx(d_star, b, '.');
hold on;
semilogx(d_star, b0, 'Color', [0.5,0.5,0.5]);
hold off;




%%
disp(' ');


%-- 1st order Tikhonov ----%
disp('Running Tikhonov (1st) ...');
lambda_tk1 = 2e1;
[x_tk1, ~, ~, Gpo_inv_tk1] = ...
    invert.tikhonov(Lb*A, Lb*b, lambda_tk1, 1);
Gpo_tk1 = inv(Gpo_inv_tk1);
e.tk1 = (x_tk1 - x0)' * Gpo_inv_tk1 * (x_tk1 - x0);
tools.textdone();
disp(' ');



%-- 1st order Tikhonov ----%
disp('Running Tikhonov (2nd) ...');
lambda_tk2 = 4e1;
[x_tk2, ~, ~, Gpo_inv_tk2] = ...
    invert.tikhonov(Lb*A, Lb*b, lambda_tk2, 1);
Gpo_tk2 = inv(Gpo_inv_tk2);
e.tk2 = (x_tk2 - x0)' * Gpo_inv_tk2 * (x_tk2 - x0);
tools.textdone();
disp(' ');



figure(2);
semilogx(d, x0, 'k--');
hold on;

x_tk = x_tk2;
Gpo_tk = Gpo_tk2;
semilogx(d, x_tk, 'c');
semilogx(d, x_tk + 2 .* sqrt(diag(Gpo_tk)),'c--');
semilogx(d, max(x_tk - 2 .* sqrt(diag(Gpo_tk)),0), 'c--');
hold off;


