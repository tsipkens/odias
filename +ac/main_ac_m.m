
% MAIN_IAC_M  A function evaluating the iterative-average charge algorithm.
%  Allows for a quick check of methods for specific m_star.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

clear;
close all;
addpath cmap tfer_pma autils;



% Defaults.
Rm = 3;
nit = 4e13;
eps = 13.5;

Dm0 = 2.48;
rho100_0 = 510;
Q = 0.9183;

charge_type = 'Fuchs';


% Set charging model parameters. 
opt.nit = nit;
opt.eps = eps;


% Set transfer function evaluation grid.
nx = 700;
m = logspace(-3.5, 3.5, nx)';  % reconstruction points
m_star = [5e-3, 0.00965967, 0.0443217, 0.101746, 1, 3]';  % mass-to-charge setpoints
nb = length(m_star);

z = 1:300;


% Get properties and then update.
prop = kernel.prop_pma;
prop = prop_update_flow(prop, Q .* 1.66667e-5);
prop = massmob.add(prop, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation

d = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % use mass-mobility relation
d_star = (m_star .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % use mass-mobility relation

sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', Rm);

tools.textheader('Computing kernel')
[K, ~, fq, Kq, qbar0] = kernel.gen_pma(sp, m, d, z, prop, [], 'Fuchs', opt);  % get kernel
tools.textheader();

% Get power law fit.
fl = qbar0 > 4;  % flag higher charge states
A = [log(d(fl)), ones(size(d(fl)))];
b = log(qbar0(fl));
p = (A' * A) \ (A' * b);
nu = p(1)
q0 = exp(p(2))
qbarh = q0 .* d .^ nu;

qbarl = ones(size(d));

n = 3;
qbart = (qbarh .^ n + 1) .^ (1/n);



%%
figure(3);
zmid = exp((log(z(2:end)) + log(z(1:(end-1)))) / 2);
zmid = [0.7071, zmid];
h = pcolor(d, zmid, fq);
hold on;
plot(d, qbar0, 'k');
plot(d, qbarh);
plot(d, qbarl);
plot(d, qbart);
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
set(h, 'EdgeColor', 'none');
cm = ocean;
colormap(flipud(cm(50:end-2,:)));
colorbar;



%%
%== AC ALGORITHMS ========================================================%
% Compute "true" average.
pr = ones(size(m));
gsmd = 1.8;  % GSD for mobility
gsm = exp(prop.zet * log(gsmd));  % GSD with respect to mass
pr = normpdf(log(m), log(0.3801), log(gsm));  % alternate size distribution
[m_t, m_t_geo] = ac.true(m_star, K, m, pr);
table(m_star, m_t, m_t_geo)

% Run IAC algorithm. 
[m_iac, q_iac] = ...
    ac.iac_m(m_star, prop, [], charge_type, opt);
table(m_star, m_iac, q_iac)

% Run FTFAC algorithm.
[m_ftfac, q_ftfac] = ac.ftfac(m_star, Kq, z);
table(m_star, m_ftfac, q_ftfac)

% Run PLAC algorithm.
[m_plac, q_plac, qfun_plac] = ac.plac(m_star, nu, q0, prop);
table(m_star, m_plac, q_plac)


figure(2);
plot(m_star, q_iac, 'or');
hold on;
plot(m_star, q_ftfac, 'o');
plot(m_star, q_plac, 'ok');
plot(m, qfun_plac(m), 'k--');
plot(m, qbar0);  % map mass to charge directly (as opposed to transmitted)
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log');


