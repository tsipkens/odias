
function p = logdiff(d, x)

p = (log(x(2:end)) - log(x(1:end-1))) ./ ...
    (log(d(2:end)) - log(d(1:end-1)));

end
