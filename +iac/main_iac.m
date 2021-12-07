
% MAIN_IAC  A function evaluating the iterative-average-charge algorithm.

clear;
close all
clc;


clear;
close all;
addpath cmap tfer_pma;

% Set charging model parameters. 
opt.nit = 4e13;
opt.eps = 13.5;

nx = 700;  nb = 400;
m = logspace(-4, 3, nx)';  % reconstruction points
m_star = logspace(-4, 1, nb)';  % mass-to-charge setpoints
Rm = 3;

% Get properties and then update.
prop = kernel.prop_pma;
Dm0 = 2.48;
rho100_0 = 510;
Q = 0.9183;
prop = prop_update_flow(prop, Q .* 1.66667e-5);
prop = working.prop_update_massmob(prop, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation


d = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % get mobility diameters

sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', Rm);
A = kernel.gen_pma(sp, m, d, 1:200, prop, [], 'Fuchs', opt);



%%
mu = [1, 0.1];
s = [2.5, 1.9];
w = [1, 0.5];

%-{
mu = 0.1;
s = 5; %2.5;
w = 1;
%}


figure(1);
[b, Lb, x0] = tools.gen_data(A, m, mu, s, w, m_star, 1e3);



%== FIG. 2 ==============================================%
%   Largely independent code to generate FIG. 2. 
disp('For <strong>FIG. 2</strong>...');
zmax = 200;
m0 = 6;  % mass for FIG. 2
d0 = (m0 .* (1:zmax) .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);
fz = kernel.tfer_charge(d0 .* 1e-9, 1:400, [], 'Fuchs', opt);
[~, qbar0, ~] = working.iac(m0, prop, [], 'Fuchs', opt);

figure(2);
cmap_sweep(zmax, tempo(zmax));
plot(fz, '-');
xlim([0, 80]);

hold on;
plot(1:zmax, diag(fz,0)', 'ko', ...
    'LineWidth', 1.2, 'MarkerFaceColor', [1,1,1]);
hold off;

limy = ylim();
hold on;
plot([1,1] .* qbar0, limy, 'k--');
hold off;
%========================================================%


[m_star_iac, qbar, d_star] = working.iac(m_star, prop, [], 'Fuchs', opt);






%%

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



e



%%


am = sum(A, 2);  am = am ./ max(am);

figure(3);
subplot(7,7,[1,12]);
plot(m, sum(A));
set(gca, 'XScale', 'log');

subplot(7,7,[15,47]);
h = pcolor(m, m_star, A);
set(gca, 'XScale', 'log', 'YScale', 'log');
set(h, 'EdgeColor', 'none');
colormap matter;

hold on;
for ii=0:10
    plot((2 ^ ii) .* m_star, m_star, 'w:');
end
hold off;

xlabel('x / m');
ylabel('b / m*');

subplot(7,7,[20,49]);
plot(sum(A, 2), m_star);
set(gca, 'YScale', 'log');



figure(4);
plot(m, x_tk2 ./ max(x_tk2));
hold on;
plot(m, x0 ./ max(x0), 'k');
plot(m_star_iac, b ./ am ./ max(b ./ am));
plot(m_star, b ./ max(b), 'k-');
hold off;
set(gca, 'XScale', 'log');

disp('Output:')
mtot1 = exp(sum(x_tk2 ./ sum(x_tk2) .* log(m)))
mtot2 = full(exp(sum(b ./ am ./ sum(b ./ am) .* log(m_star_iac))))
disp(' ')




%%
%{
% Full charging model. VERY SLOW!
[m_star_fac, qbar1, fz1] = working.fac(m_star, prop, [], [], opt);
%}



figure(5);
loglog(m_star, qbar);
hold on;

fun = @(x) (x(1) .* m_star .^ x(2) + x(3));
% x1 = lsqnonlin(@(x) log(qbar) - log(fun(x)), [26,0.67,1])
% loglog(m_star, fun(x1));
loglog(m_star, fun([6.4168, 0.4204, 0.4587]));  % previous fit

linidx = (length(qbar)-100):length(qbar);
x2 = lsqnonlin(@(x) ...
    log(qbar(linidx)) - x(1) .* log(m_star(linidx)) - x(2), [1,0])
loglog(m_star, exp(x2(1) .* log(m_star) + x2(2)), 'k--')
loglog(m_star([1,end]), [1,1], 'k--');

plot(m_star, qbar, 'k.');
hold off;



figure(6);
loglog(d_star, qbar);
hold on;

linidx = (length(qbar)-100):length(qbar);
x3 = lsqnonlin(@(x) ...
    log(qbar(linidx)) - x(1) .* log(d_star(linidx)) - x(2), [1,0])
% loglog(d_star, exp(x3(1) .* log(d_star) + x3(2)), 'k--')
loglog(d_star([1,end]), [1,1], 'k--');
loglog(d_star, 0.055 .* d_star, 'k--')

if exist('qbar1', 'var')
    loglog(d_star, qbar1);
end
hold off;

% xlim([8, 300]);
% ylim([0.9, 20]);



