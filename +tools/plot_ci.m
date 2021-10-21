
% PLOT_CI  Plot MAP estimate along with credible intervals.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-09

function h = plot_ci(d, x, Gpo, x0, cm)

if ~exist('x0', 'var'); x0 = []; end

% Plot color.
if ~exist('cm', 'var'); cm = []; end
% if isempty(cm); cm = [0.36, 0.79, 0.98]; end  % pastel blue
if isempty(cm); cm = [1, 0.37, 0.54]; end  % pastel red

% Parse error bound input.
% If not supplied, do not plot.
if ~exist('Gpo', 'var'); Gpo = []; end
if isempty(Gpo), Gpo = zeros(length(x)); end


% Credible interval limits for plotting.
x_high2 = x + 2 .* sqrt(diag(Gpo));
x_low2 = max(x - 2 .* sqrt(diag(Gpo)),0);
x_high1 = x + sqrt(diag(Gpo));
x_low1 = max(x - sqrt(diag(Gpo)),0);

% Plot shaded region.
reg = [x_low2; flipud(x_high2)];
fill([d; flipud(d)], reg, cm, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.15);

set(gca, 'XScale', 'log');
hold on;

reg = [x_low1; flipud(x_high1)];
fill([d; flipud(d)], reg, cm, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.22);


% Plot estimate.
h = plot(d, x, 'Color', cm, 'LineWidth', 1.5);


% Plot credible intervals.
% plot(d, x_high2, '--', 'Color', [cm, 0.5]);
% plot(d, x_low2, '--', 'Color', [cm, 0.5]);


% Plot true solution.
if ~isempty(x0)
    plot(d, x0, 'k--', 'LineWidth', 1);
end

hold off;
axis square;

if nargout==0; clear h; end

end

