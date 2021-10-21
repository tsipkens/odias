
function xi = get_init(A, b, d, d_star)

% Interpolate from b.
xi = interp1(d_star, ...
    full(b) ./ (A * ones(length(d),1)), ...
    d);

% Filter out unphysical outputs.
xi(isnan(xi)) = 0;
xi(isinf(xi)) = 0;
xi = max(0,xi);

end