
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
m = logspace(-4, 4, nx)';  % reconstruction points
m_star = [5e-4, 0.00965967, 0.0443217, 0.101746, 1, 3]';  % mass-to-charge setpoints
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
[~, ~, ~, Aq] = kernel.gen_pma(sp, m, d, z, prop, [], 'Fuchs', opt);  % get kernel
tools.textheader();

tools.textheader('Computing charge fractions')
[fq, qbar0] = kernel.tfer_charge(d .* 1e-9, z, 298, 'Fuchs', opt);
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
%== AC ALGORITHMS ========================================================%

% Run IAC algorithm. 
[m_iac, q_iac] = ...
    ac.iac_m2(m_star, prop, [], charge_type, opt);
table(m_star, m_iac, q_iac)

% Run FKAC algorithm.
[m_fkac, q_fkac] = ac.fkac(m_star, Aq, z);
table(m_star, m_fkac, q_fkac)

% Run PLAC algorithm.
[m_plac, q_plac, qfun_plac] = ac.plac(m_star, nu, q0, prop);
table(m_star, m_plac, q_plac)


figure(2);
plot(m_star, q_iac, 'or');
hold on;
plot(m_star, q_fkac, 'o');
plot(m_star, q_plac, 'ok');
plot(m, qfun_plac(m), 'k--');
plot(m, qbar0);  % map mass to charge directly (as opposed to transmitted)
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log');




