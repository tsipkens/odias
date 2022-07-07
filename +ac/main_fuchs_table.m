
% MAIN_FUCHS_TABLE  Generate a table of power laws for Fuchs model.
%  
%  AUTHOR: Timothy Sipkens, 2022-07-07

clear;
close all;
clc;

% Set vectors of model parameters.
eps_vec = [3, 10, 50];
nit_vec = [0.05, 0.1, 0.5, 1, 4] .* 1e13;

%-{
eps_vec = logspace(log10(3), log10(50), 30);
nit_vec = logspace(log10(0.05), log10(5), 40) .* 1e13;
%}

n_eps = length(eps_vec);
n_nit = length(nit_vec);

% Set charge fraction evaluation vector.
nx = 500;
d = logspace(0.5, 3.5, nx)';  % charge equivalent diameter

% Set integer charge vector for evaluation.
zmax = 150;
zvec = 0:zmax;

nu_vec = [];
q0_vec = [];

disp('Overall progress:');
tools.textbar([0, n_nit, 0, n_eps]);
for ii=1:n_eps
    opt.eps = eps_vec(ii);

    for jj=1:n_nit
        opt.nit = nit_vec(jj);

        % Get charge distribution.
        [fq, qbar0] = kernel.tfer_charge(d .* 1e-9, zvec, 298, 'Fuchs', opt);

        % Get power law fit.
        [nu, q0, p] = ac.get_power_law(qbar0, d);
        nu_vec(ii, jj) = nu;
        q0_vec(ii, jj) = q0;

        figure(1);
        plot(d, qbar0, 'k', 'LineWidth', 2);
        hold on;
        plot(d, q0 .* d .^ nu, 'r-');
        hold off;
        set(gca, 'XScale', 'log', 'YScale', 'log');
        ylim([0.7, inf]);  xlim([min(d), max(d)]);
        
        disp('Overall progress:');
        tools.textbar([0, n_eps, 0, n_nit]);
        tools.textbar([jj, n_nit, ii, n_eps]);
    end
end

%{
table(nu_vec)
table(q0_vec)
%}

figure(1);
imagesc(nu_vec);
colorbar;

figure(2);
imagesc(q0_vec);
colorbar;


