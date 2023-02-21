
% MAIN_AC_PERTURB  A function perturbing the AC inputs for a PMA. 
%  
%  Includes Fig. 2-4 in current draft.
%  
%  AUTHOR: Timothy Sipkens, 2022-06

clear;
close all;
addpath cmap tfer_pma autils;



% Defaults.
Rm0 = 3;
nit = 4e12;
eps0 = 13.5;

Dm0 = 2.48;
rho100_0 = 510;
Q = 0.9183;

charge_type = 'Fuchs';



% Set charging model parameters. 
opt0.nit = nit;
opt0.eps = eps0;
opt0.geo = 1;  % compute geometric mean in IAC
opt0.n = 2.5;


% Set transfer function evaluation grid.
nx = 2e3;  nb = 700;
m = logspace(-5, 2, nx)';  % reconstruction points
m_star = logspace(-4, 0, nb)';  % mass-to-charge setpoints

zvec = 0:300;

% Get properties and then update.
prop0 = kernel.prop_pma;
prop0 = prop_update_flow(prop0, Q .* 1.66667e-5);
prop0 = massmob.add(prop0, ...
	'Dm', Dm0, 'rho100', rho100_0);  % universal relation

d0 = (m .* 1e-18 ./ prop0.m0) .^ (1 / prop0.Dm);  % use mass-mobility relation
d_star = (m_star .* 1e-18 ./ prop0.m0) .^ (1 / prop0.Dm);  % use mass-mobility relation

% Get Fuchs power law.
disp('Getting power law for default...');
sp = get_setpoint(prop0, 'm_star', m_star .* 1e-18, 'Rm', Rm0);
[K0, ~, fq0, Kq0, qbar0] = kernel.gen_pma(sp, m, d0, zvec, prop0, [], 'Fuchs', opt0);  % get kernel
[nu0, qp0] = ac.get_power_law(qbar0, d0);



%%
% Kernel ignoring neutrals. 
Kq0_nn = Kq0(:,2:end,:);  % no neutrals
K0_nn = squeeze(sum(Kq0_nn, 2));
zvec_nn = zvec(2:end);

% "True" average transmitted particle size.
m_bar_t0 = ac.true([], K0_nn, m);
m_bar_g0 = exp(ac.true([], K0_nn, log(m)));

% Run the IAC algorithm with default settings. 
[m_bar_iac0, q_bar_iac0] = ...
    ac.iac_m(m_star, prop0, [], charge_type, opt0);

% Run the FK algorithm with default settings. 
[m_bar_fkac0, q_bar_fkac0] = ...
    ac.fk(m_star, Kq0_nn, m, zvec_nn);
[m_bar_fcfac0, q_bar_fcfac0] = ...
    ac.fcf(m_star, fq0, m, zvec);

% Run the PLAC algorithm with default settings. 
[m_bar_plac0, q_bar_plac0] = ...
    ac.plac(m_star, nu0, qp0, prop0);
[m_bar_intac0, q_bar_intac0] = ...
    ac.intac(m_star, nu0, qp0, prop0, opt0.n);  % instead solved with IAC


% Plot errors at default.
figure(2); clf;

subplot(5, 1, 1:2);
plot(m_star, m_bar_iac0);
hold on;
plot(m_star, m_bar_fkac0);
plot(m_star, m_bar_fcfac0);
plot(m_star, m_bar_intac0);
plot(m_star, m_bar_plac0);
plot(m_star, m_bar_t0, 'k');
plot(m_star, m_star, 'k--');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
ylim([1e-4, 1e2]);

subplot(5, 1, 3:5);
plot(m_star, (m_bar_iac0 - m_bar_t0) ./ m_bar_t0);
hold on;
plot(m_star, (m_bar_fkac0 - m_bar_t0) ./ m_bar_t0);
plot(m_star, (m_bar_fcfac0 - m_bar_t0) ./ m_bar_t0);
plot(m_star, (m_bar_intac0 - m_bar_t0) ./ m_bar_t0);
plot(m_star, (m_bar_plac0 - m_bar_t0) ./ m_bar_t0);
plot(m_star, (m_bar_g0 - m_bar_t0) ./ m_bar_t0, 'k');
plot(m_star, (m_star - m_bar_t0) ./ m_bar_t0, 'k');
yline(0, 'k--');
hold off;
ylim([-0.15, 0.15]);  % reset y-axis
set(gca, 'XScale', 'log');



% Transfer function plot.
figure(3);
h = pcolor(m_star, m, K0');
set(h, 'EdgeColor', 'none');
set(gca, 'XScale', 'log', 'YScale', 'log');

hold on;
plot(m_star, m_star, 'r-');
hold off;

cm = ocean;
colormap(flipud(cm(50:end-2,:)));



%%
%== START: Perturbation analysis =========================================%

% Range of configuration file specifying perturbation scenarios. 
cfg = tools.load_config('+ac/config/v3.sig1.json');  % fast
% cfg = tools.load_config('+ac/config/v1.mu.json');  % fast
% cfg = tools.load_config('+ac/config/v1.Rm.json');
% cfg = tools.load_config('+ac/config/v1.mm.b.json');
% cfg = tools.load_config('+ac/config/v1.nit.b.json');
% cfg = tools.load_config('+ac/config/v1.eps.b.json');

% cfg = tools.load_config('+ac/config/v2.nit.json');
% cfg = tools.load_config('+ac/config/v1.mm.rho.json');
% cfg = tools.load_config('+ac/config/v1.mm.Dm.json');

% If not a field, both = 0. 
% This indicate the the AC methods are not being perturbed. 
% In other words, the perturbed quantities do not change those methods.
if ~isfield(cfg, 'both')
    cfg.both = 0;
end
cfg  % show configuration

type = cfg.tfer_type;
prop = prop0;  % copy with universal relation


scan_vec = cfg.scan;
sz = max(size(scan_vec));


m_bar_iac = [];
m_bar_t = [];
m_bar_fkac = [];
m_bar_fcfac = [];
m_bar_plac = [];
m_bar_intac = [];
for kk=1:sz
    disp(['Running condition no. ', num2str(kk), '.']);

    if ~strcmp(cfg.perturb, 'distr')

        % For changing PMA resolution.
        if strcmp(cfg.perturb, 'Rm')
            Rm = scan_vec(kk);
        else
            Rm = Rm0;
        end
        
        % For changing charge model parameters.
        opt = opt0;
        if strcmp(cfg.perturb, 'nit')
            opt.nit = scan_vec(kk) .* 1e13;
        end
        if strcmp(cfg.perturb, 'eps0')
            opt.eps = scan_vec(kk);
        end
        
        % For changing mass-mobility parameters.
        if contains(cfg.perturb, 'mm')
            if strcmp(cfg.perturb, 'mm.rho')
                rho100 = scan_vec(kk) .* rho100_0;
                prop = massmob.add(prop0, ...
                    'Dm', Dm0, 'rho100', rho100);  % universal relation
            elseif strcmp(cfg.perturb, 'mm.Dm')
                Dm = scan_vec(kk) .* Dm0;
                prop = massmob.add(prop0, ...
                    'Dm', Dm, 'rho100', rho100_0);  % universal relation
            elseif strcmp(cfg.perturb, 'mm')
                Dm = scan_vec(1, kk);
                rho100 = scan_vec(2, kk);
                prop = massmob.add(prop0, ...
                    'Dm', Dm, 'rho100', rho100);  % universal relation
            end
        end
        
        % If changing both the IAC and full transfer function outputs. 
        %{
        if cfg.both == 0
            m_bar_iac(kk,:) = m_bar_iac0;
        else
            m_bar_iac(kk,:) = ac.iac_m(m_star, prop, [], charge_type, opt);
        end
        %}
        
        sp = get_setpoint(prop, 'm_star', m_star .* 1e-18, 'Rm', Rm);
        d = (m .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);  % use mass-mobility relation
        [K, ~, fq, Kq, qbar] = kernel.gen_pma(sp, m, d, zvec, prop, [], 'Fuchs', opt);  % get kernel
        [nu, qp] = ac.get_power_law(qbar, d);

        Kq_nn = Kq(:,2:end,:);  % no neutrals
        K_nn = squeeze(sum(Kq_nn, 2));
        zvec_nn = zvec(2:end);
        
        % "True" average transmitted particle size.
        m_bar_t(kk,:) = ac.true([], K_nn, m);
        
        if cfg.both == 1
            % Run the FK algorithm. 
            [m_bar_fkac(kk,:), q_bar_fkac] = ...
                ac.fk(m_star, Kq_nn, m, zvec_nn);
            
            % Run the FCFAC algorithm. 
            [m_bar_fcfac(kk,:), q_fcfac] = ...
                ac.fcf(m_star, fq, m, zvec);
            
            % Run the PLAC algorithm with default settings. 
            [m_bar_plac(kk,:), q_bar_plac] = ...
                ac.plac(m_star, nu, qp, prop);
            [m_bar_intac(kk,:), q_bar_intac] = ...
                ac.intac(m_star, nu, qp, prop);  % instead solved with IAC
        
        % ELSE: Perturb only the true average mass, thereby assessing what
        % happens when the assumed model parameters are incorrect.
        else
            m_bar_fkac(kk,:) = m_bar_fkac0;
            m_bar_fcfac(kk,:) = m_bar_fcfac0;
            m_bar_plac(kk,:) = m_bar_plac0;
            m_bar_intac(kk,:) = m_bar_intac0;
        end
    else
        % Only perturbing distribution, so AC methods do not change.
        % No kernel reevaluation required.

        if isfield(cfg, 'mu')
            muk = cfg.mu;
            sigk = cfg.scan(kk);
        end

        if isfield(cfg, 'sig')
            sigk = cfg.sig;
            muk = cfg.scan(kk);
        end

        sigk = exp(prop0.zet .* log(sigk));  % convert to mass
        
        % m_bar_iac(kk,:) = m_bar_iac0;
        m_bar_fkac(kk,:) = m_bar_fkac0;
        m_bar_fcfac(kk,:) = m_bar_fcfac0;
        
        p = normpdf(log(m), log(muk), log(sigk));
        if isinf(sigk)
            p = ones(size(p));
        end
        pfl = abs((log(m) - log(muk)) ./ log(sigk)) > 3;

        m_bar_t(kk,:) = ac.true([], K0_nn, m, p);
        
        t0 = sum(p' .* K0_nn, 2);  % sum(x_bar .* A_bar{1}, 2);
        x_flag = t0 < (0.05 .* max(t0));
        m_bar_t(kk, x_flag) = NaN;
        
        m_bar_plac(kk,:) = m_bar_plac0;
        m_bar_intac(kk,:) = m_bar_intac0;
    end

end




%== START figures ========================================================%
figure(11);
plot(m_star, m_bar_t, 'k');
hold on;
plot(m_star, m_bar_intac);
plot(m_star, m_star, 'k--');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');

figure(10);
plot(m_star, (m_bar_intac - m_bar_t) ./ m_bar_t);
ylim([-0.25, 0.25]);
xlim([1e-3, 1e0]);
set(gca, 'XScale', 'log');

