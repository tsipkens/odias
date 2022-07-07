
% INVERT_INTAC  A quick inversion routine using the INTAC method.
%  Acts to shift the x-axis, with little other correction.

function y = invert_intac(x_star, b, x, varargin)

tools.textheader('Inversion using INTAC');
disp(' ');
disp('NOTE: Result is only qualitative.');
disp(' ');

x_bar = ac.intac(x_star, varargin{:});

y = interp1(x_bar, full(b), x, 'linear', 0);

% Correct for changing bin width.
dx_bar = log(x_bar(2:end)) - log(x_bar(1:(end-1)));  dx_bar = [dx_bar(1); dx_bar];
dx = log(x(2:end)) - log(x(1:(end-1)));  dx = [dx(1); dx];
dy = interp1(x_bar, dx_bar, x, 'linear', 'extrap') ./ dx;
y = y ./ dy;

tools.textheader();

end
