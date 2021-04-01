
% MAIN  Inversion of mobility distributions.
% AUTHOR: Timothy Sipkens, 2020-04-11
%=========================================================================%

addpath cmap;

d = logspace(log10(10),log10(1e3),400)';  % reconstruction points
d_star = logspace(log10(10),log10(1e3),150)';  % mobility setpoints

prop_dma = kernel.prop_dma;


A = kernel.gen_dma(d_star, d, prop_dma);


mu_d = 200;
s_d = 1.4;

% Unimodal.
% x0 = normpdf(log10(d), log10(mu_d), log10(s_d));

% Bimodal.
x0 = normpdf(log10(d), log10(mu_d), log10(s_d)) + ...
    0.5 .* normpdf(log10(d), log10(mu_d / 3), log10(s_d / 1.15));

b0 = A * x0;

[b, Lb] = tools.get_noise(b0, 1e3, 1e-6);


figure(1);
semilogx(d_star, b, '.');
hold on;
semilogx(d_star, b0, 'Color', [0.5,0.5,0.5]);
hold off;


%%

disp(' ');

%-- Least-squares ---------%
disp('Running least-squares ...');
x_lsq = invert.lsq(Lb*A, Lb*b);
tools.textdone();
disp(' ');


%-- 1st order Tikhonov ----%
disp('Running Tikhonov (1st) ...');
lambda_tk1 = 3.3e1;
[x_tk1, ~, ~, Gpo_inv_tk1] = ...
    invert.tikhonov(Lb*A, Lb*b, lambda_tk1, 1);
Gpo_tk1 = inv(Gpo_inv_tk1);
e.tk1 = (x_tk1 - x0)' * Gpo_inv_tk1 * (x_tk1 - x0);
tools.textdone();
disp(' ');


%-- 2nd order Tikhonov ----%
disp('Running Tikhonov (2nd) ...');
lambda_tk2 = 5e2;
[x_tk2, ~, ~, Gpo_inv_tk2] = ...
    invert.tikhonov(Lb*A, Lb*b, lambda_tk2, 2);
Gpo_tk2 = inv(Gpo_inv_tk2);
e.tk2 = (x_tk2 - x0)' * Gpo_inv_tk2 * (x_tk2 - x0);
tools.textdone();
disp(' ');


%-- Exponential distance --%
disp('Running exponential distance ...');
lambda_ed = 5e0;
ld = 1.2 .* log10(s_d);
[x_ed, ~, ~, Gpo_inv_ed] = ...
    invert.exp_dist(Lb*A, Lb*b, lambda_ed, ld, d);
Gpo_ed = inv(Gpo_inv_ed);
e.ed = (x_ed - x0)' * Gpo_inv_ed * (x_ed - x0);
tools.textdone();
disp(' ');
disp(' ');


e



%%
figure(2);
semilogx(d, x0, 'k--');
hold on;

x_tk = x_tk2;
Gpo_tk = Gpo_tk2;
semilogx(d, x_tk, 'c');
semilogx(d, x_tk + 2 .* sqrt(diag(Gpo_tk)),'c--');
semilogx(d, max(x_tk - 2 .* sqrt(diag(Gpo_tk)),0), 'c--');

semilogx(d, x_ed, 'Color', [0.7,0.7,0.7]);
semilogx(d, x_ed + 2 .* sqrt(diag(Gpo_ed)), '--', 'Color', [0.7,0.7,0.7]);
semilogx(d, max(x_ed - 2 .* sqrt(diag(Gpo_ed)),0), '--', 'Color', [0.7,0.7,0.7]);

hold off;


