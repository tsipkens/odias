
function [resid, xnorm, lvec, xvec] = lcurve(A, b)

% lvec = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5];
% n = length(lvec);

n = 50;
lvec = logspace(-2, 5, n);

xvec = zeros(size(A, 2), n);
xvec2 = xvec;

tools.textbar([0, n]);
for ii=1:n
    [xvec(:,ii), ~, ~, Gpo_inv] = invert.tikhonov(A, b, lvec(ii), 2);
    xvec2(:,ii) = chol(inv(Gpo_inv)) * xvec(:,ii);
    
    tools.textbar([ii, n]);
end

resid = sum((A * xvec - b) .^ 2);
xnorm = sum(xvec2 .^ 2);

%-{
% Diagnostic figure.
figure(gcf);
plot(resid, xnorm, '.-');
set(gca, 'XScale', 'log', 'YScale', 'log');
text(resid, xnorm, num2str(log10(lvec)'));
%}

%{
% Alternative plot.
plot(xvec);
%}

end
