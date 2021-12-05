
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


% Get properties and then update.
prop = kernel.prop_pma;
Dm0 = 2.48;
rho100_0 = 510;
Q = 0.9183;
prop = prop_update_flow(prop, Q .* 1.66667e-5);
prop = working.prop_update_massmob(prop, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation




%%
[m_star_iac_bin] = working.iac(m_star, prop, [], 'Fuchs', opt);


%%
type_bin = '1C_diff';
prop_bin = working.prop_update_massmob(prop, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation

Rm = 3;
scan_vec = [2,3,5,10];

% span_vec = [0.5, 0.75, 1, 1.25, 1.5] .* rho100_0;

A_bar = {};
for kk=1:length(scan_vec)
    Rm = scan_vec(kk);
    
    %{
    rho100_bin = 1;
    prop_bin = working.prop_update_massmob(prop, ...
        'Dm', Dm0, 'rho100', rho100_bin);  % universal relation
    %}
    
    d_bar = (m .* 1e-18 ./ prop_bin.m0) .^ (1 / prop_bin.Dm);  % get mobility diameters
    sp_bar = get_setpoint(prop_bin, 'm_star', m_star .* 1e-18, 'Rm', Rm);
    
    A_bar{kk} = kernel.gen_pma(sp_bar, m, d_bar, (1:200)', prop_bin, type_bin, 'Fuchs', opt);
end

%
x_bar = ones(size(m))';
% scan_vec = [1/4, 1/2, 1, 2, 4, 8, 16];  % for mod. to GMD
% scan_vec = [1.2, 1.5, 2, 4, 8, 30, inf];  % for GSD


figure(82);
clf;
set(gca, 'XScale', 'log');
m_bar_ftf = [];
for kk=1:length(scan_vec)
    
    %-{
    m_bar_ftf(kk,:) = exp(sum(x_bar .* A_bar{kk} .* log(m)' ./ ...
        sum(x_bar .* A_bar{kk},2), 2));
    %}
    
    %{
    muk = 0.1;
    sigk = 1.5;
    
    % muk = scan_vec(kk) .* 0.1;
    sigk = scan_vec(kk);
    if isinf(sigk)
        x_bar = ones(size(m))';
    else
        x_bar = normpdf(log(m), log(muk), log(sigk))';
    end
    x_bar = x_bar + max(max(x_bar)) .* (1e-8);
    % x_flag = x_bar < max(max(x_bar)) .* (1e-6);
    % x_bar(x_flag) = 1e1 * length(m) * eps;
    
    hold on; plot(m, x_bar); hold off;
    
    m_bar_ftf(kk,:) = exp(sum(x_bar .* A_bar{1} .* log(m)' ./ ...
        sum(x_bar .* A_bar{1},2), 2));
    % m_bar_ftf(kk,log(m_star) > log(muk) + log(6 .* sigk)) = NaN;
    % m_bar_ftf(kk,log(m_star) < log(muk) - log(6 .* sigk)) = NaN;
    x_flag = sum(x_bar .* A_bar{2},2) < (1e8 * length(m) * eps);
    m_bar_ftf(kk,x_flag) = NaN;
    %}
    
end

figure(80);
loglog(m_star, m_bar_ftf);
hold on;
loglog(m_star, m_star_iac_bin);
hold off;

figure(81);
cm = ocean;  cm = cm(1:(end-40), :);
cmap_sweep(size(m_bar_ftf,1), cm);
semilogx(m_star, (m_bar_ftf ./ m_star_iac_bin') - 1);
hold on; semilogx(m_star, (m_bar_ftf(1,:) ./ m_star_iac_bin') - 1, 'r--'); hold off;
xlim([m_star(70), m_star(end-25)]);

mr = m_star_iac_bin ./ m_star;
m1 = find(mr > 1.000001, 1);
m2 = find(mr > 2, 1);
m3 = find(mr > 3, 1);
m4 = find(mr > 4, 1);
m8 = find(mr > 8, 1);
m16 = find(mr > 16, 1);
hold on;
xline(m_star(m1), 'r');
xline(m_star(m2), 'r');
xline(m_star(m3), 'r');
xline(m_star(m4), 'r');
xline(m_star(m8), 'r');
xline(m_star(m16), 'r');
yline(0, 'k');
hold off;

ylim([-0.3, 0.3]);
ylim([-1, 2]);


figure(35);
imagesc(log10(m_star), log10(m), (x_bar .* A_bar{2})');
set(gca, 'YDir', 'normal');
colormap(flipud(ocean));

