
clear;
close all;
addpath cmap tfer_pma autils;


% Defaults.
Rm = 3;
nit = 4e12;
eps = 13.5;

Dm0 = 2.48;
rho100_0 = 510;
Q = 0.9183;

charge_type = 'Fuchs';


% Set charging model parameters. 
opt.nit = nit;
opt.eps = eps;


% Set evaluation grid.
nx = 1000;
d = logspace(log10(1), log10(4e3), nx)';

z = 0:300;

% Compute kernel.
tools.textheader('Computing kernel')
[fq, qbar0] = kernel.tfer_charge(d .* 1e-9, z, [], 'Fuchs', opt);
tools.textheader();

% Get power law fit.
[nu, q0] = ac.get_power_law(qbar0, d);

qbarh = q0 .* d .^ nu;
qbarl = ones(size(d));


%%

% White's model.
qbar_white = working.white(d, nit);


% [qbar_ranga, d_ranga] = working.read_meancharge(...
%     '..\..\Program - Aerosol inversion\ranga-charging-code\v6');


load('+working/li_v4_collkernel.mat');
qbar_ranga = [];
for ii=1:length(dvec)
    qbar_ranga(ii) = kernel.collkernel2charge(collkernel0{ii}, nit);
end
d_ranga = dvec;
fl_r = d_ranga <= 1e3;
d_ranga = d_ranga(fl_r);
qbar_ranga = qbar_ranga(fl_r);


% Plot Fuch's model
figure(3);
subplot(5, 1, 1:3);
zmid = exp((log(z(2:end)) + log(z(1:(end-1)))) / 2);
zmid = [0.05, 0.7071, zmid(2:end)];
h = pcolor(d, zmid, fq);
set(h, 'EdgeColor', 'none');
hold on;
plot(d, qbar0, 'k');
plot(d, qbarh);
plot(d, qbarl);
plot(d, qbar_white, 'm');
plot(d_ranga, qbar_ranga, 'g-');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
cm = ocean;
colormap(flipud(cm(50:end-2,:)));
ylim([0.1, 300]);
% colorbar;


% Plot exponent as a function of size.
p_fuchs = working.logdiff(d, qbar0);
p_white = working.logdiff(d, qbar_white);
p_ranga = working.logdiff(d_ranga, qbar_ranga');

d_mid = exp((log(d(2:end)) + log(d(1:end-1))) ./ 2);
d_mid_ranga = exp((log(d_ranga(2:end)) + log(d_ranga(1:end-1))) ./ 2);


subplot(5, 1, 4:5);
plot(d_mid, p_fuchs, 'k');
hold on;
plot(d_mid, p_white, 'm');
plot(d_mid, nu .* ones(size(d_mid)), 'r');
plot(d_mid_ranga, p_ranga, 'g-');
hold off;
set(gca, 'XScale', 'log');
xlim([d(1), d(end)]);
ylim([0.7, 2]);

