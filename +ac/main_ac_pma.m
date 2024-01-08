
% MAIN_AC_M  A function evaluating the average charge algorithms for PMAs at limited m_star.
%  Allows for a quick check of methods for specific m_star.
%  
%  Includes Fig. 1 in current draft.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-26

clear;
close all;
addpath cmap autils tfer;



% Defaults.
Rm = 3;
nit = 1e12;
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
m = logspace(-4, 3.5, nx)';  % reconstruction points
m_star = [1e-3, 5e-3, 0.00965967, 0.0443217, 0.101746, 1, 3]';  % mass-to-charge setpoints
nb = length(m_star);

z = 0:300;


% Get properties and then update.
prop = prop_pma;
prop = prop_update_flow(prop, Q .* 1.66667e-5);
prop = massmob.add(prop, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation

d = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % use mass-mobility relation
d_star = (m_star .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % use mass-mobility relation

sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', Rm);

tools.textheader('Computing kernel')
[K, ~, fq, Kq, qbar0] = kernel.gen_pma(sp, m, d, z, prop, [], 'Fuchs', opt);  % get kernel
tools.textheader();

Kq_nn = Kq(:,2:end,:);  % no neutrals
K_nn = squeeze(sum(Kq_nn, 2));
z_nn = z(2:end);

% Get power law fit.
[nu, q0] = get_power_law(qbar0, d);

qbarh = q0 .* d .^ nu;
qbarl = ones(size(d));
    
n = 3;
qbart = (qbarh .^ n + 1) .^ (1/n);

xt = ac.get_transition(nu, q0, prop)


%%
% Transfer function plot.
figure(3);
h = pcolor(m_star, m, K');
set(h, 'EdgeColor', 'none');
set(gca, 'XScale', 'log', 'YScale', 'log');

hold on;
plot(m_star, m_star, 'r-o');
hold off;



%%
%== AC ALGORITHMS ========================================================%
% Compute "true" average.

pr = ones(size(m));  % no size distribution

%-{
% Alternatively, add a lognormal size distribution.
gsmd = 1.8;  % GSD for mobility
gsm = exp(prop.zet * log(gsmd));  % GSD with respect to mass
gmd = 0.3;
pr = normpdf(log(m), log(gmd), log(gsm));  % alternate size distribution
%}

[m_t, m_t_geo] = ac.true(m_star, K, m, pr);
table(m_star, m_t, m_t_geo)

% Run IAC algorithm. 
[m_iac, q_iac] = ...
    ac.iac_m(m_star, prop, [], charge_type, opt);
table(m_star, m_iac, q_iac)

% Run FK algorithm.
[m_fk, q_fk] = ac.fk(m_star, Kq, z, pr);  % including neutrals (incorrect approach)
table(m_star, m_fk, q_fk)

[m_fk_nn, q_fk_nn] = ac.fk(m_star, Kq_nn, z_nn, pr);
table(m_star, m_fk, q_fk)

% Run FCFAC algorithm.
[m_fcf, q_fcf] = ac.fcf(m_star, fq, m, z);
table(m_star, m_fcf, q_fcf)

% Run PLAC algorithm.
[m_plac, q_plac, qfun_plac] = ac.plac(m_star, nu, q0, prop);
table(m_star, m_plac, q_plac)

[m_intac, q_intac] = ac.intac(m_star, nu, q0, prop);
table(m_star, m_intac, q_intac)


figure(2);
plot(m_star, q_iac, 'sr');
hold on;
plot(m_star, q_fk_nn, 'o');
plot(m_star, q_fcf, '^');
plot(m_star, q_plac, 'ok');
plot(m, qfun_plac(m), 'k--');
plot(m_star, q_intac, '<g');
plot(m, qbar0);  % map mass to charge directly (as opposed to transmitted)
yline(1, 'k--');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');

