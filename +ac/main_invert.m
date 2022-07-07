
% MAIN_INVERT  Inversion of mobility distributions.
%  
%  AUTHOR: Timothy Sipkens, 2022-07-07

clear;
close all;
addpath cmap tfer_pma autils;

m = logspace(-3, 2, 500)';  % reconstruction points
m_star = logspace(-2.5, 0.5, 80)';  % mass-to-charge setpoints

prop = kernel.prop_pma;
prop = massmob.add(prop, 'soot');
d = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % get mobility diameters

sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', 3);
[A, ~, ~, ~, qbar0] = kernel.gen_pma(sp, m, d, 0:100, prop, [], 'Fuchs');


%%
mu = [1, 0.1];
s = [2.2, 1.5];
w = [1, 1];

%-{
mu = [1, 0.1];
s = [2.2, 1.06];
w = [1, 0.1];
%}

%{
mu = 1;
s = 2.5;
w = 1;
%}

figure(1);
clf;
[b, Lb, x0] = tools.gen_data(A, m, mu, s, w, m_star);


%%
[nu, q0] = ac.get_power_law(qbar0, d)
n = 2.5;

figure(2);
plot(d, qbar0, 'k', 'LineWidth', 2);
hold on;
plot(d, q0 .* d .^ nu, 'r-');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
ylim([0.7, inf]);  xlim([min(d), max(d)]);

eta = 1 / prop.zet;
c0 = (1 / prop.k) ^ eta;
x_intac = ac.invert_intac(m_star, b, m, nu, q0, prop, [], n);

figure(1);
subplot(2, 1, 1);
m_star2 = exp(log(m_star) - (log(m_star(2)) - log(m_star(1))) ./ 2);
stairs(m_star2, b, 'b');  % more representative of TSI display
hold on;
plot(m_star, b, 'b.');
hold off;
set(gca, 'XScale', 'log');
xlim([min(m), max(m)]);

subplot(2, 1, 2);
plot(m, x0, 'k');
hold on;
plot(m, x_intac, 'r');
hold off;
set(gca, 'XScale', 'log');
xlim([min(m), max(m)]);

