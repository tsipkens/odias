
% MAIN_IAC  A function evaluating the iterative-average-charge algorithm.


clear;
close all;
addpath cmap tfer_pma;



% Defaults.
Rm = 3;
nit = 4e13;
eps0 = 13.5;

Dm0 = 2.48;
rho100_0 = 510;
Q = 0.9183;



% Set charging model parameters. 
opt.nit = nit;
opt.eps = eps0;


% Set transfer function evaluation grid.
nx = 900;  nb = 450;
m = logspace(-4.5, 3, nx)';  % reconstruction points
m_star = logspace(-4, 1, nb)';  % mass-to-charge setpoints


% Get properties and then update.
prop0 = kernel.prop_pma;
prop0 = prop_update_flow(prop0, Q .* 1.66667e-5);
prop0 = working.prop_update_massmob(prop0, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation



% Run the IAC algorithm. 
[m_star_iac] = working.iac(m_star, prop0, [], 'Fuchs', opt);

%-{
% Full charging model. VERY SLOW for all of the points!
sel = 1:25:length(m_star);
m_star_fac = working.fac(m_star(sel), prop0, [], [], opt);
%}


%%

cfg = tools.load_config('+iac/config/v1.sig.json');
% cfg = tools.load_config('+iac/config/v1.Rm.json');


type = cfg.tfer_type;
prop = working.prop_update_massmob(prop0, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation


scan_vec = cfg.scan;

% span_vec = [0.5, 0.75, 1, 1.25, 1.5] .* rho100_0;

A_bar = {};
if ~strcmp(cfg.perturb, 'distr')
    for kk=1:length(scan_vec)
        % For changing resolution.
        if strcmp(cfg.perturb, 'Rm')
            Rm = scan_vec(kk);
        end

        % For changing density.
        if strcmp(cfg.perturb, 'rho')
            rho100 = scan_vec(kk);
            prop = working.prop_update_massmob(prop0, ...
                'Dm', Dm0, 'rho100', rho100);  % universal relation
        end

        d_bar = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % get mobility diameters
        sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', Rm);

        A_bar{kk} = kernel.gen_pma(sp, m, d_bar, (1:300)', prop, type, 'Fuchs', opt);
        % A_bar{kk} = kernel.gen_pma(sp, m, d_bar, 0', prop, type, [], opt);
    end
else
    d_bar = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % get mobility diameters
    sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', Rm);
    A_bar{1} = kernel.gen_pma(sp, m, d_bar, (1:300)', prop, type, 'Fuchs', opt);
end

%
x_bar = ones(size(m))';
% scan_vec = [1/4, 1/2, 1, 2, 4, 8, 16];  % for mod. to GMD
% scan_vec = [1.2, 1.5, 2, 4, 8, 30, inf];  % for GSD



m_bar_ftf = [];
for kk=1:length(scan_vec)
    
    if ~strcmp(cfg.perturb, 'distr')
        m_bar_ftf(kk,:) = exp(sum(x_bar .* A_bar{kk} .* log(m)' ./ ...
            sum(x_bar .* A_bar{kk},2), 2));
    
    else
        muk = 0.1;
        sigk = 1.5;

        % muk = scan_vec(kk) .* 0.1;
        sigk = scan_vec(kk);
        if isinf(sigk)
            x_bar = ones(size(m))';
        else
            x_bar = normpdf(log(m), log(muk), log(sigk))';
        end
        % x_bar = x_bar + max(max(x_bar)) .* (1e-8);
        % x_flag = x_bar < max(max(x_bar)) .* (1e-6);
        % x_bar(x_flag) = 1e1 * length(m) * eps;

        m_bar_ftf(kk,:) = exp(sum(x_bar .* A_bar{1} .* log(m)' ./ ...
            sum(x_bar .* A_bar{1},2), 2));
        m_bar_ftf(kk,:) = sum(x_bar .* A_bar{1} .* m' ./ ...
            sum(x_bar .* A_bar{1},2), 2);
        % m_bar_ftf(kk,log(m_star) > log(muk) + log(6 .* sigk)) = NaN;
        % m_bar_ftf(kk,log(m_star) < log(muk) - log(6 .* sigk)) = NaN;
        
        t0 = sum(x_bar .* A_bar{1},2);
        x_flag = t0 < (0.01 .* max(t0));
        m_bar_ftf(kk,x_flag) = NaN;
    end
    
end



figure(80);
loglog(m_star, m_bar_ftf);
hold on;
loglog(m_star, m_star_iac, 'k', 'LineWidth', 2);
loglog(m_star(sel), m_star_fac, 'ko', 'MarkerSize', 5);
hold off;


figure(81);
cm = ocean;  cm = cm(1:(end-40), :);
cmap_sweep(size(m_bar_ftf,1), cm);
semilogx(m_star, (m_star_iac' ./ m_bar_ftf) - 1);
hold on; semilogx(m_star, (m_star_iac' ./ m_bar_ftf(1,:)) - 1, 'r--'); hold off;

mr = m_star_iac ./ m_star;
m1 = find(mr > 1.000001, 1) - 1;
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

hold on;
semilogx(m_star(sel), (m_star_fac' ./ m_bar_ftf(end, sel)) - 1, 'ko', 'MarkerSize', 5);
semilogx(m_star(sel), (m_star_fac' ./ m_bar_ftf(1, sel)) - 1, 'ko', 'MarkerSize', 5);
hold off;

ylim([-0.3, 0.3]);
if strcmp(cfg.perturb, 'distr')
    ylim([-0.5, 1]);
end
xlim([10^-4, m_star(end-25)]);


figure(35);
imagesc(log10(m_star), log10(m), (x_bar .* A_bar{end})');
set(gca, 'YDir', 'normal');
colormap(flipud(ocean));
xlim([-4, log10(m_star(end-15))]);
ylim([-4, log10(m(end))]);

hold on;
plot(log10(m_star), log10(m_star_iac'), 'r-');
plot(log10(m_star(sel)), log10(m_star_fac'), 'ro', 'MarkerSize', 5);
hold off;




figure(81);  % show Fig. 81 first
