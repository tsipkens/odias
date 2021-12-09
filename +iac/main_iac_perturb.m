
% MAIN_IAC  A function evaluating the iterative-average-charge algorithm.


clear;
close all;
addpath cmap tfer_pma;



% Defaults.
Rm0 = 3;
nit = 4e13;
eps0 = 13.5;

Dm0 = 2.48;
rho100_0 = 510;
Q = 0.9183;



% Set charging model parameters. 
opt0.nit = nit;
opt0.eps = eps0;


% Set transfer function evaluation grid.
nx = 900;  nb = 400;
m = logspace(-5, 3, nx)';  % reconstruction points
m_star = logspace(-4, 1, nb)';  % mass-to-charge setpoints


% Get properties and then update.
prop0 = kernel.prop_pma;
prop0 = prop_update_flow(prop0, Q .* 1.66667e-5);
prop0 = working.prop_update_massmob(prop0, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation



% Run the IAC algorithm. 
m_star_iac0 = working.iac(m_star, prop0, [], 'Fuchs', opt0);

%-{
% Full charging model. VERY SLOW for all of the points!
sel = 1:25:length(m_star);
m_star_fac = working.fac(m_star(sel), prop0, [], [], opt0);
%}


%%

% cfg = tools.load_config('+iac/config/v1.default.json');

% cfg = tools.load_config('+iac/config/v1.sig.json');
% cfg = tools.load_config('+iac/config/v1.Rm.json');
% cfg = tools.load_config('+iac/config/v1.mm.rho.json');
% cfg = tools.load_config('+iac/config/v1.mm.Dm.json');
% cfg = tools.load_config('+iac/config/v1.nit.json');
cfg = tools.load_config('+iac/config/v1.nit.b.json');

% If not a field, indicate the the full transfer function
% and IAC methods are not being perturbed (only the full
% trasnfer function approach). 
if ~isfield(cfg, 'both')
    cfg.both = 0;
end
cfg  % show configuration

type = cfg.tfer_type;
prop = working.prop_update_massmob(prop0, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation


scan_vec = cfg.scan;

% span_vec = [0.5, 0.75, 1, 1.25, 1.5] .* rho100_0;

m_star_iac = [];
A_bar = {};
if ~strcmp(cfg.perturb, 'distr')
    for kk=1:length(scan_vec)
        % For changing resolution.
        if strcmp(cfg.perturb, 'Rm')
            Rm = scan_vec(kk);
        else
            Rm = Rm0;
        end
        
        if strcmp(cfg.perturb, 'nit')
            opt = opt0;
            opt.nit = scan_vec(kk) .* 1e13;
        else
            opt = opt0;
        end
        
        % For changing density.
        if contains(cfg.perturb, 'mm')
            if strcmp(cfg.perturb, 'mm.rho')
                rho100 = scan_vec(kk) .* rho100_0;
                prop = working.prop_update_massmob(prop0, ...
                    'Dm', Dm0, 'rho100', rho100);  % universal relation
            elseif strcmp(cfg.perturb, 'mm.Dm')
                Dm = scan_vec(kk) .* Dm0;
                prop = working.prop_update_massmob(prop0, ...
                    'Dm', Dm, 'rho100', rho100_0);  % universal relation
            end
        end
    
        % If changing both the IAC and full transfer 
        % function outputs. 
        if cfg.both == 0
            m_star_iac(kk,:) = m_star_iac0;
        else
            m_star_iac(kk,:) = working.iac(m_star, prop, [], 'Fuchs', opt);
        end
        
        
        d_bar = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % get mobility diameters
        sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', Rm);

        A_bar{kk} = kernel.gen_pma(sp, m, d_bar, (1:300)', prop, type, 'Fuchs', opt);
        % A_bar{kk} = kernel.gen_pma(sp, m, d_bar, 0', prop, type, [], opt);
    end
else
    for kk=1:length(scan_vec)
        m_star_iac(kk,:) = m_star_iac0;
    end
    
    d_bar = (m .* 1e-18 ./ prop0.m0) .^ (1 / prop0.Dm);  % get mobility diameters
    sp = get_setpoint(prop0, 'm_star', m_star .* 1e-18, 'Rm', Rm);
    A_bar{1} = kernel.gen_pma(sp, m, d_bar, (1:300)', prop0, type, 'Fuchs', opt0);
end



m_bar_ftf = [];
m_bar_jtf = [];
for kk=1:length(scan_vec)
    
    if ~strcmp(cfg.perturb, 'distr')
        m_bar_jtf(kk,:) = exp(sum(A_bar{kk} .* log(m)' ./ ...
            sum(A_bar{kk},2), 2));
        m_bar_ftf(kk,:) = m_bar_jtf(kk,:);
    
    else
        m_bar_jtf(kk,:) = exp(sum(A_bar{1} .* log(m)' ./ ...
            sum(A_bar{1},2), 2));
        
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
        % m_bar_ftf(kk,:) = sum(x_bar .* A_bar{1} .* m' ./ ...
        %     sum(x_bar .* A_bar{1},2), 2);  % arithmatic mean
        % m_bar_ftf(kk,log(m_star) > log(muk) + log(6 .* sigk)) = NaN;
        % m_bar_ftf(kk,log(m_star) < log(muk) - log(6 .* sigk)) = NaN;
        
        t0 = sum(x_bar .* A_bar{1},2);
        x_flag = t0 < (0.01 .* max(t0));
        m_bar_ftf(kk,x_flag) = NaN;
    end
    
end

% Check if transfer function not fully resolved
% and NaN if truncated.
x_flag2 = logical([]);
for kk=1:length(scan_vec)
    if length(A_bar) == 1; kk2 = 1;
    else; kk2 = kk; end
    x_flag2(kk,:) = (A_bar{kk2}(:,1) > (0.001 .* max(A_bar{kk2},[],2)));
    m_bar_ftf(kk,x_flag2(kk,:)) = NaN;
end


% Fit to IAC. 
p = [];
m_bar_co1 = [];
for kk=1:length(scan_vec)
    xflag3 = m_star > 0.3;
    xflag3(end-20:end) = 0;
    
    p(kk,:) = polyfit(log(m_star(xflag3)), log(m_star_iac(kk,xflag3)), 1);
    m_bar_co1(kk, :) = exp(polyval(p(kk,:), log(m_star)));
    
    m_bar_co2(kk, :) = m_star;
end


figure(80);
loglog(m_star, m_bar_ftf);
hold on;
loglog(m_star, m_star_iac, 'k', 'LineWidth', 1);
% loglog(m_star(sel), m_star_fac, 'ko', 'MarkerSize', 5);
loglog(m_star, m_bar_co1, 'r:');
hold off;


figure(81);
cm = ocean;  cm = cm(1:(end-40), :);
cmap_sweep(size(m_bar_ftf,1), cm);
semilogx(m_star, (m_star_iac ./ m_bar_ftf) - 1);
hold on; semilogx(m_star, (m_star_iac(1,:) ./ m_bar_ftf(1,:)) - 1, 'r--'); hold off;

mr = m_star_iac(1,:)' ./ m_star;
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

%{
hold on;
semilogx(m_star(sel), (m_star_fac' ./ m_bar_ftf(end, sel)) - 1, 'ko', 'MarkerSize', 5);
semilogx(m_star(sel), (m_star_fac' ./ m_bar_ftf(1, sel)) - 1, 'ko', 'MarkerSize', 5);
hold off;
%}

%-{
hold on;
semilogx(m_star, (m_bar_co1(end,:) ./ m_bar_ftf(end,:)) - 1, 'r', 'MarkerSize', 5);
semilogx(m_star, (m_bar_co2(end,:) ./ m_bar_ftf(end,:)) - 1, 'r', 'MarkerSize', 5);
hold off;
%}

%{
if strcmp(cfg.perturb, 'distr')
    hold on;
    semilogx(m_star(:), (m_bar_jtf(end, :) ./ m_bar_ftf(end, :)) - 1, 'r-');
    hold off;
end
%}

ylim([-0.3, 0.3]);
if or(strcmp(cfg.perturb, 'distr'), strcmp(cfg.perturb, 'mm.Dm'))
    ylim([-0.5, 0.5]);
end
xlim([10^-4, 1]);


figure(35);
imagesc(log10(m_star), log10(m), (x_bar .* A_bar{end})');
set(gca, 'YDir', 'normal');
colormap(flipud(ocean));
xlim([-4, log10(m_star(end-15))]);
ylim([-4, log10(m(end))]);

hold on;
plot(log10(m_star), log10(m_star_iac'), 'r-');
plot(log10(m_star(sel)), log10(m_star_fac'), 'ro', 'MarkerSize', 5);
plot(log10(m_star), log10(m_bar_co1'), 'r-');
hold off;




figure(81);  % show Fig. 81 first


