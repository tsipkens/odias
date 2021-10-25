
function [b, Lb, x0] = gen_data(A, d, mu, s, w, d_star)

% Weighting for multuple modes. 
if ~exist('w', 'var'); w = []; end
if isempty(w); w = ones(size(mu)); end

if ~exist('d_star', 'var'); d_star = []; end
if isempty(d_star); f_plot = 0;
else; f_plot = 1; end


% Unimodal.
% x0 = normpdf(log10(d), log10(mu_d), log10(s_d));

% Bimodal.
% x0 = normpdf(log10(d), log10(mu), log10(s)) + ...
%     0.5 .* normpdf(log10(d), log10(mu / 3), log10(s / 1.15));

x0 = zeros(size(d));
for ii=1:length(mu)
    x0 = x0 + w(ii) .* ...
        normpdf(log10(d), log10(mu(ii)), log10(s(ii)));
end

b0 = A * x0;

N = 1e3;
[b, Lb] = tools.get_noise(b0, N, 1e-6);


% Plot data and distribution.
if f_plot
    d_star2 = exp(log(d_star) - (log(d_star(2)) - log(d_star(1))) ./ 2);  % shifted by 1/2 setpoint
    
    figure(gcf);
    plot(d_star, b0, '--', 'Color', 0.6.*[1,1,1]);
    hold on;
    % plot(d_star, b, 'b.', 'MarkerSize', 10);
    stairs(d_star2, b, 'b');  % more representative of TSI display

    hold off;
    set(gca, 'XScale', 'log');
    title('Real distribution and data');

    limy = ylim();
    hold on;
    eb = errorbar(d_star, b, 2 ./ (diag(Lb)), '.', 'Color', 0.85.*[1,1,1]);
    hold off;
    ylim([0, limy(2)]);
    uistack(eb,'bottom')
end

end
