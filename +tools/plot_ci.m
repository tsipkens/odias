
% PLOT_CI  Plot MAP estimate along with credible intervals.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-09

function h = plot_ci(d, x, Gpo, x0, cm)

if ~exist('x0', 'var'); x0 = []; end

% Plot color.
if ~exist('cm', 'var'); cm = []; end
% if isempty(cm); cm = [0.36, 0.79, 0.98]; end  % pastel blue
if isempty(cm); cm = [1, 0.37, 0.54]; end  % pastel red

% Credible interval limits for plotting.
x_high = x + 2 .* sqrt(diag(Gpo));
x_low = max(x - 2 .* sqrt(diag(Gpo)),0);

% Plot shaded region.
reg = [x_low; flipud(x_high)];
fill([d; flipud(d)], reg, cm, 'EdgeColor', 'none');
alpha(0.12);
set(gca, 'XScale', 'log');
hold on;

% Plot estimate.
h = plot(d, x, 'Color', cm, 'LineWidth', 1.5);

% Plot credible intervals.
plot(d, x_high, '--', 'Color', [cm, 0.5]);
plot(d, x_low, '--', 'Color', [cm, 0.5]);

% Plot true solution.
if ~isempty(x0)
    plot(d, x0, 'k--');
end

hold off;

if nargout==0; clear h; end

end

